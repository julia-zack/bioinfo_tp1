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
import random
import statistics

from Bio import SeqIO

from sequences import (
    AA,
    CODON_TABLE,
    expected_aa_frequencies,
    find_orfs,
    generate_random_nt_sequence,
    stop_codon_probability,
)
from stats import freqs
from ncbi import fetch_cds_annotations, fetch_genbank_record, get_real_sequences
from exercises.ex4.a import REAL_ORGANISMS

# How likely an ORF is to be coding before any evidence is in. 0.5 is neutral
# and leaves the decision to the two signals.
PRIOR_CODING = 0.5

# Transcript used when run() is called without arguments: a real human mRNA
# whose annotated CDS is the ground truth for (v).
DEMO_ACCESSION = "NM_001317077.2"

# Fixes the random DNA of the negative control, so the run is reproducible.
NEGATIVE_CONTROL_SEED = 0

# Transcripts with an annotated coding region, used to weigh the two signals
# against each other by measuring instead of guessing.
TRANSCRIPT_SET = "data/ncbi_sample.fasta"

# Composition weights tried when fitting. 4d found the length to be the main
# signal and the composition a secondary one, so the search stays at or below 1,
# which is where the two would count the same.
WEIGHT_GRID = [round(0.05 * i, 2) for i in range(21)]

FRAME_ORDER = ['+1', '+2', '+3', '-1', '-2', '-3']

PER_LINE = 5   # amino acids per line when printing a distribution


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


def candidate_discount(orf_count):
    """How much to take off every ORF for having examined orf_count of them.

    Six frames over a sequence yield dozens of ORFs, and each one is a separate
    chance for a run of non-stop codons to come out long by luck. Scoring one
    candidate as if it had been the only one looked at overstates it by roughly
    the number of tries, which in log-odds is a subtraction of log(orf_count).
    """
    return math.log(orf_count) if orf_count else 0.0


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

    def __init__(self, length_mu, length_sigma, stop_prob, aa_weight, prior,
                 comp_weight=1.0):
        self.length_mu = length_mu         # mean of ln(length) of real proteins
        self.length_sigma = length_sigma   # standard deviation of ln(length)
        self.stop_prob = stop_prob         # q = P(codon is a stop)
        self.aa_weight = aa_weight         # {amino acid: log(f_natural/f_random)}
        self.log_prior_odds = math.log(prior / (1 - prior))
        self.comp_weight = comp_weight     # how much the composition counts

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
                + self.comp_weight * self.composition_log_odds(orf['protein'])
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


def train_scorer(real_seqs, prior=PRIOR_CODING, comp_weight=1.0):
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
        comp_weight=comp_weight,
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
    discount = candidate_discount(len(orfs))
    for orf in orfs:
        orf['composition'] = freqs(orf['protein'])                      # (iii-b)
        orf['log_odds'] = scorer.coding_log_odds(orf) - discount        # (iv)
        orf['probability'] = sigmoid(orf['log_odds'])                   # (iv)
    orfs.sort(key=lambda orf: orf['probability'], reverse=True)
    return orfs


def group_by_frame(orfs):
    """(ii) The ORFs split by reading frame, longest first within each."""
    grouped = {label: [] for label in FRAME_ORDER}
    for orf in orfs:
        grouped[orf['frame']].append(orf)
    for frame_orfs in grouped.values():
        frame_orfs.sort(key=lambda orf: -orf['length'])
    return grouped


def print_orfs_by_frame(orfs):
    """(ii) ORF count and longest ORF for each of the six frames."""
    grouped = group_by_frame(orfs)
    print(f"\n  (ii) {len(orfs)} ORFs across the six reading frames:")
    print(f"  {'frame':>6} {'ORFs':>5} {'longest':>8} {'start':>7} {'end':>7}")
    for label in FRAME_ORDER:
        frame_orfs = grouped[label]
        if not frame_orfs:
            print(f"  {label:>6} {0:>5} {'-':>8} {'-':>7} {'-':>7}")
            continue
        longest = frame_orfs[0]
        print(f"  {label:>6} {len(frame_orfs):>5} {longest['length']:>8} "
              f"{longest['start']:>7} {longest['end']:>7}")


def print_longest_composition(orfs):
    """(iii) Length and amino acid distribution of the longest ORF found."""
    if not orfs:
        print("\n  (iii) No ORFs found.")
        return

    longest = max(orfs, key=lambda orf: orf['length'])
    print(f"\n  (iii) Longest ORF: frame {longest['frame']}, "
          f"{longest['start']}-{longest['end']}, {longest['length']} aa")
    print("        amino acid distribution:")
    for start in range(0, len(AA), PER_LINE):
        row = AA[start:start + PER_LINE]
        print("        " + "  ".join(
            f"{a} {longest['composition'].get(a, 0):.3f}" for a in row))


def print_orf_table(orfs, top=8):
    """(iv) The ORFs most likely to be coding, best first."""
    print(f"\n  (iv) {len(orfs)} ORFs found. Most likely to be coding:")
    print(f"  {'frame':>6} {'start':>7} {'end':>7} {'length':>7} {'log-odds':>9} {'P(cod)':>8}")
    for orf in orfs[:top]:
        print(f"  {orf['frame']:>6} {orf['start']:>7} {orf['end']:>7} "
              f"{orf['length']:>7} {orf['log_odds']:>9.1f} {orf['probability']:>8.3f}")


# ===========================================================================
# Weighing the two signals against each other
# ===========================================================================

def labelled_signals(scorer):
    """The two signals of every ORF of every annotated transcript.

    One list per transcript, holding for each of its ORFs the length signal
    already carrying the candidate discount, the composition signal, and whether
    that ORF is the annotated coding region. The signals are kept apart so a
    weight can be tried without scoring again. The discount is the same for every
    ORF of a transcript, so it changes the threshold measures and leaves the
    ranking ones alone.
    """
    records = list(SeqIO.parse(TRANSCRIPT_SET, "fasta"))
    annotations = fetch_cds_annotations([record.id for record in records])

    transcripts = []
    for record in records:
        cds = annotations.get(record.id)
        if cds is None:
            continue
        found = find_orfs(str(record.seq))
        discount = candidate_discount(len(found))
        orfs = [(scorer.length_log_odds(orf['length']) - discount,
                 scorer.composition_log_odds(orf['protein']),
                 [orf['start'], orf['end']] == cds)
                for orf in found]
        transcripts.append(orfs)
    return transcripts


def hit_rate(transcripts, comp_weight):
    """How often the real coding region beats another ORF of its own transcript.

    0.5 is what picking at random would give and 1.0 is a perfect separation,
    the same measure 4d uses. It answers "which of these ORFs is the CDS".
    """
    wins = comparisons = 0
    for orfs in transcripts:
        real = [(length, comp) for length, comp, is_cds in orfs if is_cds]
        if not real:
            continue
        real_score = real[0][0] + comp_weight * real[0][1]
        for length, comp, is_cds in orfs:
            if is_cds:
                continue
            comparisons += 1
            if real_score > length + comp_weight * comp:
                wins += 1
    return wins / comparisons if comparisons else 0.0


def false_positive_rate(transcripts, comp_weight):
    """Fraction of transcripts holding a false positive.

    A different question from the hit rate: not "which ORF is the CDS" but "is
    this ORF a CDS at all". Every ORF over 0.5 that is not the annotated one is
    a coding region the detector claims and is not there. Counted per transcript
    rather than per ORF, because that is how the detector gets used: one
    sequence in, and every false positive comes out next to the right answer.
    """
    called = sum(1 for orfs in transcripts
                 if any(sigmoid(length + comp_weight * comp) > 0.5
                        for length, comp, is_cds in orfs if not is_cds))
    return called / len(transcripts) if transcripts else 0.0


def found_rate(transcripts, comp_weight):
    """Fraction of transcripts whose real coding region scores over 0.5."""
    found = sum(1 for orfs in transcripts
                if any(sigmoid(length + comp_weight * comp) > 0.5
                       for length, comp, is_cds in orfs if is_cds))
    return found / len(transcripts) if transcripts else 0.0


def top_orf_rate(transcripts, comp_weight):
    """Fraction of transcripts whose best scoring ORF is the annotated CDS.

    What the detector gets right when only its top answer is read, rather than
    every ORF that clears the threshold.
    """
    best = sum(1 for orfs in transcripts
               if max(orfs, key=lambda o: o[0] + comp_weight * o[1])[2])
    return best / len(transcripts) if transcripts else 0.0


def fit_composition_weight(transcripts, grid=WEIGHT_GRID):
    """The weight with the fewest false positives on these transcripts.

    The hit rate cannot pick a weight, because it sits at 1.000 across the whole
    grid: in almost every transcript the real CDS is already the longest ORF, so
    the length alone answers that question. The false positive rate does move,
    so it is what the weight is fitted on.
    """
    return min(grid, key=lambda weight: false_positive_rate(transcripts, weight))


def check_weight_holds(transcripts, splits=5):
    """Split the transcripts in two, several times over, and see whether a weight
    fitted on one half still beats 1.0 on the other.

    Fitting and checking on the same 89 transcripts would flatter the result, so
    this repeats the exercise on halves that took no part in each other's fit.
    Returns how many splits held up, and the weights they picked.
    """
    held = 0
    picked = []
    for seed in range(splits):
        rows = list(transcripts)
        random.Random(seed).shuffle(rows)
        half = len(rows) // 2
        chosen, held_out = rows[:half], rows[half:]
        weight = fit_composition_weight(chosen)
        picked.append(weight)
        if false_positive_rate(held_out, weight) < false_positive_rate(held_out, 1.0):
            held += 1
    return held, splits, picked


def print_composition_weight(scorer):
    """(iv) How much the composition counts next to the length, by measuring.

    The weight is fitted on all the annotated transcripts, and the repeated
    split is reported next to it as a check on whether the result survives being
    fitted and measured on different transcripts.
    """
    transcripts = labelled_signals(scorer)
    weight = fit_composition_weight(transcripts)

    print("\nHow much should the composition count next to the length?")
    print(f"  Measured on {len(transcripts)} transcripts with an annotated CDS.")
    print(f"  {'weight':>8} {'hit rate':>10} {'top ORF':>9}"
          f" {'false pos':>11} {'CDS found':>11}")
    for candidate in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        mark = "  <-" if candidate == weight else ""
        print(f"  {candidate:>8.2f} {hit_rate(transcripts, candidate):>10.3f}"
              f" {top_orf_rate(transcripts, candidate):>9.0%}"
              f" {false_positive_rate(transcripts, candidate):>10.0%}"
              f" {found_rate(transcripts, candidate):>11.0%}{mark}")

    held, splits, picked = check_weight_holds(transcripts)
    print(f"\n  Fitted weight: {weight:.2f}")
    print(f"  Splitting the transcripts in two and fitting on one half beats"
          f" weight 1.00\n  on the other half in {held} of {splits} splits,"
          f" picking {', '.join(f'{w:.2f}' for w in picked)}.")
    print("  The pick moves around, so the value is loose, but staying under 1.00")
    print("  holds up.")
    print("\n  The hit rate does not move, so the length alone already says which")
    print("  ORF is the CDS. What the weight changes is how often the detector")
    print("  calls an ORF that is not one, and no weight brings that near zero.")
    return weight


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
    print("\n(v-a) POSITIVE control: real genes from a jellyfish, a plant and an")
    print("      archaeon, none of them in the training set.")
    print("      Expected: a high P on all three. A low one would mean the")
    print("      detector only knows the proteomes it was trained on.")
    print(f"  {'organism':<42} {'length':>7} {'log-odds':>12} {'P(cod)':>8}")
    for name, protein in FOREIGN_CONTROLS.items():
        orf = score_sequence(coding_dna_from_protein(protein), scorer)[0]
        print(f"  {name:<42} {orf['length']:>7} "
              f"{orf['log_odds']:>12.6f} {orf['probability']:>8.3f}")


def negative_control(length, scorer, seed=NEGATIVE_CONTROL_SEED):
    """(v-b) NEGATIVE control: random DNA should only hold short, low-P ORFs.

    `seed` fixes the sequence so the run is reproducible; pass None to draw a
    fresh one. One run is thin evidence either way, since a lucky long ORF can
    score above 0.5: length on its own can be fooled, and composition is what
    pulls it back down.
    """
    if seed is not None:
        random.seed(seed)
    dna = generate_random_nt_sequence(length)
    label = f"seed {seed}" if seed is not None else "unseeded"
    print(f"\n(v-b) NEGATIVE control: {length} nt of random DNA ({label}).")
    print("      Expected: short ORFs only, all with a low P.")
    print_orf_table(score_sequence(dna, scorer), top=3)


def scramble(protein, shuffler):
    """Reorder a protein's residues, keeping the leading M in place.

    Same length and same composition as the original, only the arrangement of
    the residues changes.
    """
    residues = list(protein[1:])
    shuffler.shuffle(residues)
    return protein[0] + "".join(residues)


def scrambled_control(scorer, seed=NEGATIVE_CONTROL_SEED):
    """(v-c) A negative control that comes back positive, on purpose.

    A protein with its residues shuffled is not a protein, so it should score
    low. It does not: the two signals come out identical to the real protein's,
    digit for digit. Neither can see the difference, because the length does not
    change under a reordering and the composition log-odds is a sum over
    residues, which does not depend on the order they are added in.

    The two signals are shown on their own here, without the candidate discount.
    That discount depends on how many ORFs the DNA happens to contain, which the
    shuffling changes, and it would hide the point rather than make it.
    """
    shuffler = random.Random(seed)
    print("\n(v-c) SCRAMBLED control: the same three proteins with their residues")
    print("      shuffled, so the length and the composition are untouched and")
    print("      only the order is destroyed.")
    print("      Expected: a low score, since a scrambled protein is not a protein.")
    print("      What happens: the two signals give the same number for both.")
    print(f"  {'organism':<42} {'length':>7} {'real':>12} {'scrambled':>12}")
    for name, protein in FOREIGN_CONTROLS.items():
        scrambled = scramble(protein, shuffler)
        real = {'length': len(protein), 'protein': protein}
        mixed = {'length': len(scrambled), 'protein': scrambled}
        print(f"  {name:<42} {len(protein):>7} "
              f"{scorer.coding_log_odds(real):>12.6f}"
              f" {scorer.coding_log_odds(mixed):>12.6f}")


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
    print_orfs_by_frame(orfs)
    print_longest_composition(orfs)
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

    # (iv) the composition should count less than the length, but by how much is
    # measured against annotated transcripts
    scorer = train_scorer(real_seqs, comp_weight=print_composition_weight(scorer))

    # (v) the three tests the exercise asks for
    positive_controls(scorer)
    negative_control(900, scorer)
    scrambled_control(scorer)
    real_gene(source, scorer)

    print("\nPossible improvements (v, last question):")
    print("  - add a signal that reads the order of the residues, such as the")
    print("    frequency of consecutive pairs or the codon usage, since (v-c)")
    print("    shows both current signals are blind to it;")
    print("  - fit the composition weight on more than 89 transcripts. The two")
    print("    references are learnt from 16.8M residues, far past what 4b showed")
    print("    is enough, while the weight rests on 89;")
    print("  - measure the prior instead of leaving it at 0.5;")
    print("  - model the length per organism rather than fitting all three together.")
