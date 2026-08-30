"""4c) Elija una metrica que permita comparar dos distribuciones. (p.ej., RMSE).

Metrica elegida: distancia de variacion total, implementada en
stats.distributions_distance() y usada por 4a, 4b y 4d.
"""

from stats import distributions_distance


def run():
    """4c is a design decision, not a computation: it names the metric."""
    print(__doc__.strip())
    print("\n    TV(p, q) = 0.5 * sum |p(a) - q(a)|\n")
    print("Bounded between 0 and 1, reads as the fraction of probability mass")
    print("that has to move to turn one distribution into the other, and can be")
    print("restricted to a subset of amino acids by summing fewer terms.")
    print("\nThe justification is written up in informe.md.")

    # A worked example, so the number is not just asserted.
    p = {'A': 0.5, 'C': 0.5}
    q = {'A': 0.4, 'C': 0.6}
    print(f"\nExample: TV({p}, {q}) = {distributions_distance(p, q):.3f}")
