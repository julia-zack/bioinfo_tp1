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
METODO (como se calcula la probabilidad del punto iv)

La idea es un naive-Bayes: cada ORF trae dos evidencias independientes -- su
LARGO y su COMPOSICION de aminoacidos -- y a cada una le preguntamos "esto se
parece mas a un gen real o a ADN al azar?". La respuesta de cada evidencia es
un log-odds (el log de un cociente de probabilidades): positivo tira a
codificante, negativo a azar. Como son independientes, los dos log-odds SE
SUMAN, y una sigmoide convierte esa suma en una probabilidad entre 0 y 1.

  1) SEÑAL DE LARGO.
     Un ORF al azar dura hasta que aparece un STOP; como P(codon=stop) = q ~ 3/64,
     su largo sigue una geometrica de media 1/q ~ 21 codones: los ORFs al azar
     son cortos. Los genes reales codifican proteinas de cientos de aminoacidos,
     y esos largos son aproximadamente log-normales (ln del largo es una
     campana). Entonces:
         P(largo | codificante) = log-normal ajustada a las proteinas reales
         P(largo | azar)        = geometrica con q
     log_odds_largo = log( P(largo|codificante) / P(largo|azar) )
     Cruza cero cerca de ~90 aa: mas largo -> mas evidencia de codificante.

  2) SEÑAL DE COMPOSICION.
     Cada aminoacido tiene un peso = log( f_natural / f_random ), donde f_natural
     es su frecuencia en proteinas reales y f_random la que produce el ADN al
     azar (degeneracion del codigo). El peso es positivo si el aminoacido es mas
     tipico de proteina real (E, K, D...) y negativo si es mas tipico del azar
     (R, C...). El log-odds de composicion de un ORF es la SUMA de los pesos de
     sus residuos: acumula evidencia a medida que el ORF es mas largo.

  3) COMBINAR.
     log_odds_total = log_odds_largo + log_odds_composicion + log(prior_odds)
     P(codificante) = sigmoide(log_odds_total) = 1 / (1 + e^-log_odds_total)

Nota honesta: mezclar una densidad continua (log-normal) con una discreta
(geometrica) en el cociente del largo es una aproximacion; y las dos señales no
son perfectamente independientes (los ORFs largos tienden a verse mas naturales
tambien). Sirve para rankear y decidir, no como probabilidad calibrada al 100%.

Las referencias (log-normal del largo y pesos de composicion) se aprenden de las
mismas proteinas reales de Swiss-Prot que usa 4a, cacheadas en data/.
"""

import math
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
from ncbi import fetch_genbank_record, get_real_sequences
from exercises.ex4.a import REAL_ORGANISMS

# A priori, antes de mirar nada, cuan probable es que un ORF cualquiera sea
# codificante. 0.5 = neutro (no inclinamos la balanza de entrada); se podria
# estimar de datos, pero neutro deja que las dos señales hablen solas.
PRIOR_CODING = 0.5

# Transcripto usado cuando run() se llama sin argumentos: un mRNA humano real,
# con su region codificante (CDS) anotada, que sirve de "verdad" para el punto v.
DEMO_ACCESSION = "NM_001317077.2"


# ===========================================================================
# Funciones matematicas basicas
# ===========================================================================

def sigmoid(x):
    """Convierte un log-odds (de -inf a +inf) en una probabilidad (0 a 1)."""
    if x < -700:            # evita overflow de e^-x para log-odds muy negativos
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def lognormal_density(x, mu, sigma):
    """Densidad log-normal en x, con parametros mu y sigma sobre ln(x)."""
    return (1.0 / (x * sigma * math.sqrt(2 * math.pi))) * \
           math.exp(-((math.log(x) - mu) ** 2) / (2 * sigma ** 2))


def geometric_density(k, q):
    """P(un ORF al azar mida exactamente k codones): k-1 codones no-stop y luego
    un stop."""
    return (1 - q) ** (k - 1) * q


# ===========================================================================
# El detector: aprende de las proteinas reales y puntua un ORF
# ===========================================================================

class CodingScorer:
    """Detector entrenado. Aprende dos cosas de las proteinas reales:

      - que tan largas son las regiones codificantes de verdad (log-normal),
        contra el largo de un ORF al azar (geometrica);
      - que aminoacidos son tipicos de proteina real y cuales del azar (los pesos).

    Con eso convierte cualquier ORF en una probabilidad de ser codificante.
    """

    def __init__(self, length_mu, length_sigma, stop_prob, aa_weight, prior):
        self.length_mu = length_mu        # media de ln(largo) de proteinas reales
        self.length_sigma = length_sigma  # desvio de ln(largo)
        self.stop_prob = stop_prob         # q = P(codon es stop), para el nulo
        self.aa_weight = aa_weight          # {aminoacido: log(f_natural/f_random)}
        self.log_prior_odds = math.log(prior / (1 - prior))

    def length_log_odds(self, length):
        """(señal 1) log( P(largo|codificante) / P(largo|azar) )."""
        coding = lognormal_density(length, self.length_mu, self.length_sigma)
        chance = geometric_density(length, self.stop_prob)
        return math.log(coding / chance)

    def composition_log_odds(self, protein):
        """(señal 2) suma de los pesos de cada residuo del ORF."""
        return sum(self.aa_weight.get(residue, 0.0) for residue in protein)

    def coding_log_odds(self, orf):
        """(iv) las dos señales sumadas, mas el prior."""
        return (self.length_log_odds(orf['length'])
                + self.composition_log_odds(orf['protein'])
                + self.log_prior_odds)

    def probability(self, orf):
        """(iv) la probabilidad final de que el ORF sea codificante."""
        return sigmoid(self.coding_log_odds(orf))


# ---------------------------------------------------------------------------
# Entrenamiento: construir el detector a partir de proteinas reales
# ---------------------------------------------------------------------------

def load_real_proteins():
    """Las proteinas reales de cada organismo (E. coli, levadura, humano)."""
    return {name: get_real_sequences(organism, cache)
            for name, (organism, cache) in REAL_ORGANISMS.items()}


def coding_lengths(real_seqs):
    """(opcion A) El largo de cada proteina real ES el largo de una region
    codificante. Junta todos: es la referencia de "largos codificantes reales"."""
    lengths = []
    for proteins in real_seqs.values():
        for protein in proteins:
            lengths.append(len(protein))
    return lengths


def natural_frequencies(real_seqs):
    """Frecuencia de cada aminoacido en proteinas reales, en "partes iguales":
    promedia los tres organismos, para que el humano (que tiene muchas mas
    proteinas) no domine el perfil."""
    per_organism = [freqs("".join(proteins)) for proteins in real_seqs.values()]
    return {aa: sum(freq[aa] for freq in per_organism) / len(per_organism)
            for aa in AA}


def amino_acid_weights(real_seqs):
    """El peso de cada aminoacido = log( f_natural / f_random ). Positivo: vota
    codificante; negativo: vota azar."""
    f_natural = natural_frequencies(real_seqs)
    f_random = expected_aa_frequencies()
    return {aa: math.log(f_natural[aa] / f_random[aa]) for aa in AA}


def train_scorer(real_seqs, prior=PRIOR_CODING):
    """Ajusta las dos referencias y devuelve el detector listo para usar."""
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
    """Muestra que aprendio el detector, para que no sea una caja negra."""
    typical = math.exp(scorer.length_mu)
    print("\nDetector entrenado:")
    print(f"  largo: log-normal  mu={scorer.length_mu:.2f}  sigma={scorer.length_sigma:.2f}"
          f"  (largo tipico ~{typical:.0f} aa)")
    print(f"  nulo:  geometrica  q={scorer.stop_prob:.4f}  (ORF al azar ~{1/scorer.stop_prob:.0f} codones)")
    ordered = sorted(AA, key=lambda a: scorer.aa_weight[a], reverse=True)
    votes = lambda residues: "  ".join(f"{a}{scorer.aa_weight[a]:+.2f}" for a in residues)
    print(f"  votan codificante: {votes(ordered[:5])}")
    print(f"  votan azar:        {votes(ordered[-5:])}")


# ===========================================================================
# (i-iv) Aplicar el detector a una secuencia de ADN
# ===========================================================================

def score_sequence(dna, scorer):
    """(ii-iv) Todos los ORFs de los 6 marcos, cada uno con su probabilidad.

    find_orfs() ya da (ii) los ORFs de los seis marcos y (iii-a) su largo;
    freqs() da (iii-b) su composicion; el scorer da (iv) la probabilidad.
    """
    orfs = find_orfs(dna)
    for orf in orfs:
        orf['composition'] = freqs(orf['protein'])       # (iii-b)
        orf['probability'] = scorer.probability(orf)      # (iv)
    orfs.sort(key=lambda orf: orf['probability'], reverse=True)
    return orfs


def print_orf_table(orfs, top=8):
    """Los ORFs mas probables primero."""
    print(f"\n  {len(orfs)} ORFs encontrados. Mas probables de ser codificantes:")
    print(f"  {'frame':>6} {'start':>7} {'end':>7} {'largo':>6} {'log-odds':>9} {'P(cod)':>8}")
    for orf in orfs[:top]:
        print(f"  {orf['frame']:>6} {orf['start']:>7} {orf['end']:>7} "
              f"{orf['length']:>6} {scorer_log_odds(orf):>9.1f} {orf['probability']:>8.3f}")


def scorer_log_odds(orf):
    """El log-odds guardado, reconstruido de la probabilidad, solo para imprimir."""
    p = min(max(orf['probability'], 1e-12), 1 - 1e-12)
    return math.log(p / (1 - p))


# ===========================================================================
# (v) Controles y gen real
# ===========================================================================

# Un codon representativo por aminoacido, para "back-translation": convertir una
# proteina real en el ADN que la codifica (control positivo).
CODON_FOR = {}
for _codon, _residue in CODON_TABLE.items():
    if _residue != '*':
        CODON_FOR.setdefault(_residue, _codon)


def coding_dna_from_protein(protein):
    """El ADN que codifica una proteina, con un stop al final: convierte una
    proteina real en el gen que la produce (back-translation)."""
    return "".join(CODON_FOR[residue] for residue in protein) + "TAA"


# Proteinas reales de organismos AJENOS al entrenamiento: no son E. coli, ni
# levadura, ni humano. Son el control positivo HONESTO: como el detector se
# entreno con esos tres organismos, testear con una de sus propias proteinas
# seria circular (la proteina de test ayudo a definir la referencia). Si en
# cambio el detector le da P alta a genes de una medusa, una planta y una
# arqueobacteria, mostramos que GENERALIZA y no memoriza. Secuencias reales de
# UniProt/Swiss-Prot.
FOREIGN_CONTROLS = {
    "Medusa - GFP (P42212)":
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFS"
        "YGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDF"
        "KEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADYQQNTPIGDGPVLLP"
        "DNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "Planta - Arabidopsis thaliana":
        "MPPKRNFRKRSFEEEEEDNDVNKAAISEEEEKRRLALEEVKFLQKLRERKLGIPALSSTAAQSSI"
        "GKVKPVEKTETEGEKEELVLQDTFAQETAVLIEDPNMVKYIEQELAKKRGRNIDDAEEVENELKR"
        "VEDELYKIPDHLKVKKRSSEESSTQWTTGIAEVQLPIEYKLKNIEETEAAKKLLQERRLMGRPKS"
        "EFSIPSSYSADYFQRGKDYAEKLRREHPELYKDRGGPQADGEAAKPSTSSSTNNNADSGKSRQAA"
        "TDQIMLERFRKRERNRVMRR",
    "Arquea - Methanocaldococcus jannaschii":
        "MQLRLSSGNVLNEKVHKVGIIALGSFLENHGAVLPIDTDIKIASYIALKASILTGAKFLGVVIPS"
        "TEYEYVKHGIHNKPEEVYSYMRFLINEGKKIGVEKFLIVNCHGGNILVESFLKDLEYEFDIKVEM"
        "INITFTHASTEEVSVGYIIGIAKADEETLKEHNNFEKYPEVGMVGLKEARENNKAIDKEAKVVKR"
        "FGVKLDKKLGEKILNNAIEKVVEKIKEMIR",
}


def positive_controls(scorer):
    """(v-a) Control POSITIVO: genes de organismos que el detector nunca vio.
    Cada proteina real se convierte en su gen y se puntua; todos deberian dar
    P alta, mostrando que el detector generaliza mas alla del entrenamiento."""
    print("\n(v-a) Control POSITIVO (organismos ajenos al entrenamiento):")
    print(f"  {'organismo':<40} {'largo':>6} {'P(cod)':>8}")
    for name, protein in FOREIGN_CONTROLS.items():
        orf = score_sequence(coding_dna_from_protein(protein), scorer)[0]
        print(f"  {name:<40} {orf['length']:>6} {orf['probability']:>8.3f}")


def negative_control(length, scorer):
    """(v-b) Control NEGATIVO: ADN al azar. Solo deberia tener ORFs cortos y con
    P baja.

    Nota para el informe: el control negativo no tiene semilla fija, asi que el
    ADN al azar cambia en cada corrida. En alguna corrida, por azar, puede
    aparecer un ORF de ~120 aa que llega a P~0.26: sigue por debajo de 0.5 (o
    sea, bien clasificado como no codificante), pero es un lindo ejemplo de que
    el largo solo a veces se deja enganar y la composicion lo baja. Si se
    quieren numeros reproducibles, alcanza con ponerle una semilla al azar
    (una linea: random.seed(...) antes de generar la secuencia).
    """
    dna = generate_random_nt_sequence(length)
    print(f"\n(v-b) Control NEGATIVO: {length} nt de ADN al azar")
    print_orf_table(score_sequence(dna, scorer), top=3)


def annotated_cds(record):
    """La region codificante anotada del gen (start, end en el marco +), que es
    la "verdad" contra la que comparamos la prediccion."""
    for feature in record.features:
        if feature.type == "CDS" and feature.location.strand == 1:
            return int(feature.location.start), int(feature.location.end)
    return None


def real_gene(accession, scorer):
    """(v) Un gen eucariota real: lo levanta, lo puntua, y compara el ORF mas
    probable con la region codificante anotada."""
    record = fetch_genbank_record(accession)
    dna = str(record.seq)
    print(f"\n(v) Gen real: {record.id}  {len(dna)} nt")
    print(f"    {record.description}")

    orfs = score_sequence(dna, scorer)
    print_orf_table(orfs)

    cds = annotated_cds(record)
    if cds is None:
        print("    (el registro no trae un CDS anotado en el marco +)")
        return
    best = orfs[0]
    match = (best['start'], best['end']) == cds
    print(f"\n    CDS anotado (la verdad):     start={cds[0]}  end={cds[1]}")
    print(f"    ORF mas probable (prediccion): start={best['start']}  end={best['end']}"
          f"  P={best['probability']:.3f}")
    print(f"    -> {'ACIERTO: el ORF mas probable es el CDS real.' if match else 'no coincide exactamente.'}")


# ===========================================================================
# run()
# ===========================================================================

def run(source=DEMO_ACCESSION):
    real_seqs = load_real_proteins()
    scorer = train_scorer(real_seqs)
    describe_scorer(scorer)

    # (v) las tres pruebas de la consigna
    positive_controls(scorer)
    negative_control(900, scorer)
    real_gene(source, scorer)

    print("\nMejoras posibles (v, ultima pregunta):")
    print("  - usar un codon-usage real para el back-translation del control positivo;")
    print("  - modelar el largo por organismo (procariota vs eucariota) en vez de uno solo;")
    print("  - calibrar la probabilidad contra CDS anotados (medir, no asumir el prior);")
    print("  - sumar señales: uso de codones, marco con ORF mas largo, señal de Kozak.")
