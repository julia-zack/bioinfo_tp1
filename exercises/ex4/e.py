"""4e) Detector de ORFs. Escriba un codigo que:

  i)   levante una secuencia de ADN
  ii)  obtenga los 6 marcos de lectura posibles y determine los ORFs de cada uno
  iii) para cada ORF determine a) la longitud y b) la distribucion de aminoacidos
  iv)  en base a lo analizado determine, a partir del largo y la distribucion,
       una probabilidad de corresponder (o no) a una region codificante
  v)   para probar el programa disene: a) un control positivo, b) un control
       negativo; luego obtenga de alguna base de datos la secuencia de un gen
       eucariota completo y corra su programa. Compare con lo esperado.
       Que modificaciones podria agregar para mejorar su capacidad predictiva?

Piezas ya disponibles:
  sequences.get_six_frames_from_nt_seq  -> (ii)
  sequences.get_orf_sizes               -> (ii), (iii-a)
  sequences.calculate_aa_frequencies    -> (iii-b)
  stats.distributions_distance          -> insumo para (iv)
  ncbi.fetch_random_records             -> (v)
"""


import os

from Bio import SeqIO

from ncbi import fetch_cds_annotations, fetch_genbank_record
from sequences import AA, find_orfs, stop_codon_probability
from stats import freqs

# Transcript used when run() is called without arguments.
DEMO_ACCESSION = "NM_001317077.2"

FRAME_ORDER = ['+1', '+2', '+3', '-1', '-2', '-3']


def file_format(path):
    return "genbank" if path.endswith((".gb", ".gbk")) else "fasta"


def load_dna(source):
    """(i) Load one DNA sequence, from a local file or an NCBI accession.

    Returns a SeqRecord, so whatever annotations came with it stay attached.
    """
    if not os.path.exists(source):
        return fetch_genbank_record(source)

    records = list(SeqIO.parse(source, file_format(source)))
    if len(records) != 1:
        raise ValueError(
            f"{source} has {len(records)} records; load_dna reads one. "
            f"Use load_dna_set() for a file with several.")
    return records[0]


def load_dna_set(path):
    """(i) Load every DNA sequence in a file, for the evaluation sets."""
    return list(SeqIO.parse(path, file_format(path)))


def annotated_cds(record):
    """The CDS coordinates annotated in the record, if it has exactly one.

    0-based and half-open. Returns None when there is no single CDS to check
    a prediction against.
    """
    cds = [f for f in record.features if f.type == "CDS"]
    if len(cds) != 1:
        return None
    location = cds[0].location
    return int(location.start), int(location.end), location.strand


def orfs_by_frame(nt_seq):
    """(ii) The ORFs of the sequence, grouped by reading frame.

    Each ORF also carries its amino acid distribution (iii-b); its length
    (iii-a) comes from find_orfs().
    """
    grouped = {label: [] for label in FRAME_ORDER}
    for orf in find_orfs(nt_seq):
        orf['composition'] = freqs(orf['protein'])
        grouped[orf['frame']].append(orf)
    for orfs in grouped.values():
        orfs.sort(key=lambda o: -o['length'])
    return grouped


def print_orfs_by_frame(grouped):
    total = sum(len(orfs) for orfs in grouped.values())
    print(f"\n{total} ORFs across the six reading frames:\n")
    print(f"{'frame':>6} {'ORFs':>5} {'longest':>8} {'start':>7} {'end':>7}")
    for label in FRAME_ORDER:
        orfs = grouped[label]
        if not orfs:
            print(f"{label:>6} {0:>5} {'-':>8} {'-':>7} {'-':>7}")
            continue
        longest = orfs[0]
        print(f"{label:>6} {len(orfs):>5} {longest['length']:>8} "
              f"{longest['start']:>7} {longest['end']:>7}")


PER_LINE = 5   # amino acids per line when printing a distribution


def print_longest_composition(grouped):
    """(iii) Length and amino acid distribution of the longest ORF found."""
    orfs = [orf for frame_orfs in grouped.values() for orf in frame_orfs]
    if not orfs:
        print("\nNo ORFs found.")
        return

    longest = max(orfs, key=lambda o: o['length'])
    print(f"\nLongest ORF: frame {longest['frame']}, "
          f"{longest['start']}-{longest['end']}, {longest['length']} aa")
    print("  amino acid distribution:")
    for start in range(0, len(AA), PER_LINE):
        row = AA[start:start + PER_LINE]
        print("   " + "  ".join(f"{a} {longest['composition'][a]:.3f}" for a in row))

def expected_by_chance(orf_length, seq_length):
    """(iv) How many ORFs this long random DNA would produce on its own.


    TODO: reescribir esto
    A codon is a stop with probability q, so a run of `orf_length` codons
    escapes one with probability (1-q)**orf_length. Six frames over a sequence
    of `seq_length` give about 2*seq_length places for such a run to start, and
    the two multiplied are the number expected by chance.

    Below 1 the ORF is unlikely to be an accident; well above 1 it is the kind
    of thing random DNA produces routinely.
    """
    q = stop_codon_probability()
    positions = max(2 * seq_length, 1)
    return positions * (1 - q) ** orf_length


def score_orfs(grouped, seq_length, table):
    """(iv) Attach to each ORF its chance score and its probability of coding."""
    for frame_orfs in grouped.values():
        for orf in frame_orfs:
            orf['expected'] = expected_by_chance(orf['length'], seq_length)
            orf['probability'] = probability_of_coding(orf['expected'], table)
    return grouped


def print_scores(grouped, top=6):
    """(iv) The most likely coding ORFs, best first."""
    orfs = sorted((o for f in grouped.values() for o in f),
                  key=lambda o: o['expected'])
    print(f"\nMost likely to be coding:\n")
    print(f"{'frame':>6} {'start':>7} {'end':>7} {'aa':>5} "
          f"{'expected by chance':>19} {'P(coding)':>11}")
    for orf in orfs[:top]:
        print(f"{orf['frame']:>6} {orf['start']:>7} {orf['end']:>7} "
              f"{orf['length']:>5} {orf['expected']:>19.4g} "
              f"{orf['probability']:>11.3f}")

# Buckets of "expected by chance", from surprising to routine. The probability
# of each is measured, not assumed: see calibrate().
CALIBRATION_EDGES = [0.01, 0.1, 1, 10, 100]

TRANSCRIPT_SET = "data/ncbi_sample.fasta"


def labelled_orfs(records, annotations):
    """Every ORF of every annotated transcript, tagged as CDS or not.

    A transcript with a known CDS labels all of its own ORFs: one is the real
    coding region and the rest are spurious, which is exactly the mix a
    detector has to sort out.
    """
    labelled = []
    for record in records:
        cds = annotations.get(record.id)
        if cds is None:
            continue
        seq = str(record.seq)
        for frame_orfs in orfs_by_frame(seq).values():
            for orf in frame_orfs:
                labelled.append((expected_by_chance(orf['length'], len(seq)),
                                 orf['start'] == cds[0] and orf['end'] == cds[1]))
    return labelled


def calibrate(labelled):
    """(iv) Probability that an ORF is coding, per bucket of its score.

    Counted from labelled ORFs: within each bucket, the fraction that turned
    out to be the annotated CDS.
    """
    edges = CALIBRATION_EDGES + [float('inf')]
    table = []
    low = 0.0
    for high in edges:
        bucket = [is_cds for score, is_cds in labelled if low <= score < high]
        probability = sum(bucket) / len(bucket) if bucket else 0.0
        table.append((high, probability, len(bucket)))
        low = high
    return table


def probability_of_coding(expected, table):
    """(iv) Look up an ORF's score in the calibration table."""
    for high, probability, _ in table:
        if expected < high:
            return probability
    return 0.0


def print_calibration(table):
    print("\nProbability of being the annotated CDS, by score:\n")
    print(f"{'expected by chance':>20} {'ORFs':>7} {'P(coding)':>11}")
    low = 0.0
    for high, probability, count in table:
        label = f"{low:g} to {high:g}" if high != float('inf') else f"over {low:g}"
        print(f"{label:>20} {count:>7} {probability:>11.3f}")
        low = high

def run(source=DEMO_ACCESSION):
    record = load_dna(source)                       # (i)
    print(f"{record.id}  {len(record.seq)} nt")
    print(f"  {record.description}")

    grouped = orfs_by_frame(str(record.seq))        # (ii) y (iii)
    print_orfs_by_frame(grouped)

    print_longest_composition(grouped)

    transcripts = load_dna_set(TRANSCRIPT_SET)      # (iv)
    annotations = fetch_cds_annotations([r.id for r in transcripts])
    table = calibrate(labelled_orfs(transcripts, annotations))
    print_calibration(table)

    score_orfs(grouped, len(record.seq), table)
    print_scores(grouped)
