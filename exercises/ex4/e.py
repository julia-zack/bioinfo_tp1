"""4e) Detector de ORFs. Escriba un codigo que:

  i)   levante una secuencia de ADN
  ii)  obtenga los 6 marcos de lectura posibles y determine los ORFs de cada uno
  iii) para cada ORF determine a) la longitud y b) la distribucion de aminoacidos
  iv)  en base a lo analizado determine, a partir del largo y la distribucion,
       una probabilidad de corresponder (o no) a una region codificante
  v)   para probar el programa disene: a) un control positivo, b) un control
       negativo; luego obtenga de una base de datos un gen eucariota completo y
       corra su programa. Compare con lo esperado.

--------------------------------------------------------------------------------
METHOD (how the probability of (iv) is computed)

Naive Bayes. Every ORF carries two pieces of evidence, its LENGTH and its amino
acid COMPOSITION, and each one is asked "does this look more like a real gene or
like random DNA?". Each answer is a log-odds (the log of a ratio of
probabilities): positive leans coding, negative leans chance. Taking them as
independent, the two log-odds ADD, and a sigmoid turns the sum into a
probability between 0 and 1.

  1) LENGTH SIGNAL.
     A random ORF runs until a stop appears; with P(codon = stop) = q ~ 3/64 its
     length follows a geometric of mean 1/q ~ 21 codons, so random ORFs are
     short. Real genes code for proteins of hundreds of residues, and those
     lengths are roughly log-normal (ln of the length is a bell curve). So:
         P(length | coding) = log-normal fitted to the real proteins
         P(length | chance) = geometric with q
     length_log_odds = log( P(length|coding) / P(length|chance) )
     It crosses zero near 90 aa: longer means more evidence of coding.

  2) COMPOSITION SIGNAL.
     Each amino acid gets a weight = log( f_natural / f_random ), where f_natural
     is its frequency in real proteins and f_random the one random DNA produces
     through the degeneracy of the genetic code. The weight is positive when the
     residue is more typical of a real protein (E, K, D...) and negative when it
     is more typical of chance (R, C...). An ORF's composition log-odds is the
     SUM of the weights of its residues, so it accumulates evidence with length.

  3) COMBINE.
     total_log_odds = length_log_odds + composition_log_odds + log(prior_odds)
     P(coding) = sigmoid(total_log_odds) = 1 / (1 + e^-total_log_odds)

Note: mixing a continuous density (log-normal) with a discrete one (geometric)
in the length ratio is an approximation, and the two signals are not perfectly
independent (longer ORFs also tend to look more natural). It is good enough to
rank and decide, not as a fully calibrated probability.

Both references (the log-normal of the length and the composition weights) are
learnt from the same Swiss-Prot proteins 4a uses, cached in data/.
"""

import math
import statistics

from sequences import (
    AA,
    CODON_TABLE,
    expected_aa_frequencies,
    find_orfs,
    generate_random_nt_sequence,
    stop_codon_probability,
)
from stats import freqs
from ncbi import fetch_genbank_record, get_real_sequences
from exercises.ex4.a import REAL_ORGANISMS

# How likely an ORF is to be coding before any evidence is in. 0.5 is neutral
# and leaves the decision to the two signals.
PRIOR_CODING = 0.5

# Transcript used when run() is called without arguments: a real human mRNA
# whose annotated CDS is the ground truth for (v).
DEMO_ACCESSION = "NM_001317077.2"


# ===========================================================================
# Basic maths
# ===========================================================================

def sigmoid(x):
    """Turn a log-odds (-inf to +inf) into a probability (0 to 1)."""
    if x < -700:            # e^-x overflows below this
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def lognormal_density(x, mu, sigma):
    """Log-normal density at x, with mu and sigma taken over ln(x)."""
    return (1.0 / (x * sigma * math.sqrt(2 * math.pi))) * \
           math.exp(-((math.log(x) - mu) ** 2) / (2 * sigma ** 2))


def geometric_density(k, q):
    """P(a random ORF is exactly k codons long): k-1 non-stops, then a stop."""
    return (1 - q) ** (k - 1) * q


# ===========================================================================
# The detector: learns from real proteins and scores an ORF
# ===========================================================================

class CodingScorer:
    """Trained detector. It learns two things from the real proteins:

      - how long real coding regions are (log-normal), against the length of a
        random ORF (geometric);
      - which amino acids are typical of a real protein and which of chance
        (the weights).

    Between them they turn any ORF into a probability that it codes.
    """

    def __init__(self, length_mu, length_sigma, stop_prob, aa_weight, prior):
        self.length_mu = length_mu         # mean of ln(length) of real proteins
        self.length_sigma = length_sigma   # standard deviation of ln(length)
        self.stop_prob = stop_prob         # q = P(codon is a stop)
        self.aa_weight = aa_weight         # {amino acid: log(f_natural/f_random)}
        self.log_prior_odds = math.log(prior / (1 - prior))

    def length_log_odds(self, length):
        """(signal 1) log( P(length|coding) / P(length|chance) )."""
        coding = lognormal_density(length, self.length_mu, self.length_sigma)
        chance = geometric_density(length, self.stop_prob)
        return math.log(coding / chance)

    def composition_log_odds(self, protein):
        """(signal 2) The weights of the ORF's residues, added up."""
        return sum(self.aa_weight.get(residue, 0.0) for residue in protein)

    def coding_log_odds(self, orf):
        """(iv) Both signals plus the prior."""
        return (self.length_log_odds(orf['length'])
                + self.composition_log_odds(orf['protein'])
                + self.log_prior_odds)

    def probability(self, orf):
        """(iv) The ORF's final probability of being coding."""
        return sigmoid(self.coding_log_odds(orf))


# ---------------------------------------------------------------------------
# Training: build the detector from real proteins
# ---------------------------------------------------------------------------

def load_real_proteins():
    """The real proteins of each organism (E. coli, yeast, human)."""
    return {name: get_real_sequences(organism, cache)
            for name, (organism, cache) in REAL_ORGANISMS.items()}


def coding_lengths(real_seqs):
    """Every real protein's length.

    A protein is what a coding region produces, so these lengths are what a real
    coding region looks like.
    """
    lengths = []
    for proteins in real_seqs.values():
        for protein in proteins:
            lengths.append(len(protein))
    return lengths


def natural_frequencies(real_seqs):
    """Frequency of each amino acid in real proteins.

    The three organisms are averaged rather than pooled, so human, which brings
    far more proteins, does not dominate the profile.
    """
    per_organism = [freqs("".join(proteins)) for proteins in real_seqs.values()]
    return {aa: sum(freq[aa] for freq in per_organism) / len(per_organism)
            for aa in AA}


def amino_acid_weights(real_seqs):
    """Weight of each amino acid, log( f_natural / f_random ).

    A positive weight votes coding, a negative one votes chance.
    """
    f_natural = natural_frequencies(real_seqs)
    f_random = expected_aa_frequencies()
    return {aa: math.log(f_natural[aa] / f_random[aa]) for aa in AA}


def train_scorer(real_seqs, prior=PRIOR_CODING):
    """Fit both references and return a detector ready to score ORFs."""
    log_lengths = [math.log(length) for length in coding_lengths(real_seqs)]
    length_mu = statistics.mean(log_lengths)
    length_sigma = statistics.pstdev(log_lengths)
    return CodingScorer(
        length_mu=length_mu,
        length_sigma=length_sigma,
        stop_prob=stop_codon_probability(),
        aa_weight=amino_acid_weights(real_seqs),
        prior=prior,
    )


def describe_scorer(scorer):
    """Print what the detector learnt, so it is not a black box."""
    typical = math.exp(scorer.length_mu)
    print("\nTrained detector:")
    print(f"  length: log-normal  mu={scorer.length_mu:.2f}  sigma={scorer.length_sigma:.2f}"
          f"  (typical length ~{typical:.0f} aa)")
    print(f"  null:   geometric   q={scorer.stop_prob:.4f}"
          f"  (random ORF ~{1/scorer.stop_prob:.0f} codons)")
    ordered = sorted(AA, key=lambda a: scorer.aa_weight[a], reverse=True)
    votes = lambda residues: "  ".join(f"{a}{scorer.aa_weight[a]:+.2f}" for a in residues)
    print(f"  vote coding: {votes(ordered[:5])}")
    print(f"  vote chance: {votes(ordered[-5:])}")


# ===========================================================================
# (i-iv) Apply the detector to a DNA sequence
# ===========================================================================

def score_sequence(dna, scorer):
    """(ii-iv) Every ORF of the six frames, each with its probability.

    find_orfs() gives (ii) the ORFs of the six frames and (iii-a) their length;
    freqs() gives (iii-b) their composition; the scorer gives (iv) the
    probability.
    """
    orfs = find_orfs(dna)
    for orf in orfs:
        orf['composition'] = freqs(orf['protein'])        # (iii-b)
        orf['probability'] = scorer.probability(orf)      # (iv)
    orfs.sort(key=lambda orf: orf['probability'], reverse=True)
    return orfs


def print_orf_table(orfs, top=8):
    """(iv) The ORFs most likely to be coding, best first."""
    print(f"\n  (iv) {len(orfs)} ORFs found. Most likely to be coding:")
    print(f"  {'frame':>6} {'start':>7} {'end':>7} {'length':>7} {'log-odds':>9} {'P(cod)':>8}")
    for orf in orfs[:top]:
        print(f"  {orf['frame']:>6} {orf['start']:>7} {orf['end']:>7} "
              f"{orf['length']:>7} {scorer_log_odds(orf):>9.1f} {orf['probability']:>8.3f}")


def scorer_log_odds(orf):
    """Recover the log-odds from the probability, for display."""
    p = min(max(orf['probability'], 1e-12), 1 - 1e-12)
    return math.log(p / (1 - p))


# ===========================================================================
# (v) Controls and a real gene
# ===========================================================================

# One codon per amino acid, enough to back-translate a protein into DNA for
# the positive control.
CODON_FOR = {}
for _codon, _residue in CODON_TABLE.items():
    if _residue != '*':
        CODON_FOR.setdefault(_residue, _codon)


def coding_dna_from_protein(protein):
    """Back-translate a protein into DNA, with a stop codon at the end."""
    return "".join(CODON_FOR[residue] for residue in protein) + "TAA"


# Real proteins from organisms outside the training set: not E. coli, not yeast,
# not human. Scoring one of those three organisms' own proteins would be
# circular, since they defined the references. If a jellyfish, a plant and an
# archaeon still get a high P, the detector generalises instead of memorising.
# Sequences from UniProt/Swiss-Prot.
FOREIGN_CONTROLS = {
    "Jellyfish - GFP (P42212)":
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFS"
        "YGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDF"
        "KEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADYQQNTPIGDGPVLLP"
        "DNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "Plant - Arabidopsis thaliana":
        "MPPKRNFRKRSFEEEEEDNDVNKAAISEEEEKRRLALEEVKFLQKLRERKLGIPALSSTAAQSSI"
        "GKVKPVEKTETEGEKEELVLQDTFAQETAVLIEDPNMVKYIEQELAKKRGRNIDDAEEVENELKR"
        "VEDELYKIPDHLKVKKRSSEESSTQWTTGIAEVQLPIEYKLKNIEETEAAKKLLQERRLMGRPKS"
        "EFSIPSSYSADYFQRGKDYAEKLRREHPELYKDRGGPQADGEAAKPSTSSSTNNNADSGKSRQAA"
        "TDQIMLERFRKRERNRVMRR",
    "Archaeon - Methanocaldococcus jannaschii":
        "MQLRLSSGNVLNEKVHKVGIIALGSFLENHGAVLPIDTDIKIASYIALKASILTGAKFLGVVIPS"
        "TEYEYVKHGIHNKPEEVYSYMRFLINEGKKIGVEKFLIVNCHGGNILVESFLKDLEYEFDIKVEM"
        "INITFTHASTEEVSVGYIIGIAKADEETLKEHNNFEKYPEVGMVGLKEARENNKAIDKEAKVVKR"
        "FGVKLDKKLGEKILNNAIEKVVEKIKEMIR",
}


def positive_controls(scorer):
    """(v-a) POSITIVE control: genes from organisms the detector never saw.

    Each real protein goes back to DNA and through the detector. A high P across
    all three shows it reaches beyond the proteomes it was trained on.
    """
    print("\n(v-a) POSITIVE control (organisms outside the training set):")
    print(f"  {'organism':<42} {'length':>7} {'P(cod)':>8}")
    for name, protein in FOREIGN_CONTROLS.items():
        orf = score_sequence(coding_dna_from_protein(protein), scorer)[0]
        print(f"  {name:<42} {orf['length']:>7} {orf['probability']:>8.3f}")


def negative_control(length, scorer):
    """(v-b) NEGATIVE control: random DNA should only hold short, low-P ORFs.

    The sequence is not seeded, so it changes every run. Some runs turn up a long
    ORF that scores high, which is worth watching: length on its own can be
    fooled, and composition is what pulls it back down.
    """
    dna = generate_random_nt_sequence(length)
    print(f"\n(v-b) NEGATIVE control: {length} nt of random DNA")
    print_orf_table(score_sequence(dna, scorer), top=3)


def annotated_cds(record):
    """The gene's annotated coding region, start and end on the + strand.

    This is the ground truth the prediction gets checked against.
    """
    for feature in record.features:
        if feature.type == "CDS" and feature.location.strand == 1:
            return int(feature.location.start), int(feature.location.end)
    return None


def real_gene(accession, scorer):
    """(v) Fetch a real eukaryotic gene, score it, and check the most probable
    ORF against the annotated coding region."""
    record = fetch_genbank_record(accession)
    dna = str(record.seq)
    print(f"\n(v) Real gene: {record.id}  {len(dna)} nt")
    print(f"    {record.description}")

    orfs = score_sequence(dna, scorer)
    print_orf_table(orfs)

    cds = annotated_cds(record)
    if cds is None:
        print("    (the record carries no CDS annotated on the + strand)")
        return
    best = orfs[0]
    match = (best['start'], best['end']) == cds
    print(f"\n    annotated CDS (truth):      start={cds[0]}  end={cds[1]}")
    print(f"    most probable ORF (call):   start={best['start']}  end={best['end']}"
          f"  P={best['probability']:.3f}")
    print(f"    -> {'HIT: the most probable ORF is the real CDS.' if match else 'no exact match.'}")


# ===========================================================================
# run()
# ===========================================================================

def run(source=DEMO_ACCESSION):
    real_seqs = load_real_proteins()
    scorer = train_scorer(real_seqs)
    describe_scorer(scorer)

    # (v) the three tests the exercise asks for
    positive_controls(scorer)
    negative_control(900, scorer)
    real_gene(source, scorer)

    print("\nPossible improvements (v, last question):")
    print("  - use a real codon usage to back-translate the positive control;")
    print("  - model the length per organism (prokaryote vs eukaryote) instead of one fit;")
    print("  - calibrate the probability against annotated CDS instead of assuming a prior;")
    print("  - add signals: codon usage, longest ORF per frame, Kozak sequence.")
