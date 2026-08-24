"""4c) Elija una metrica que permita comparar dos distribuciones (p.ej., RMSE).

Metrica elegida: distancia de variacion total (total variation distance),
implementada en stats.distributions_distance():

    TV(p, q) = 0.5 * sum_a |p(a) - q(a)|

Se eligio sobre RMSE porque:
  - Esta acotada entre 0 y 1, asi que el valor se interpreta directamente.
  - Se lee como "que fraccion de la masa de probabilidad hay que mover para
    convertir una distribucion en la otra".
  - No depende de la cantidad de bins, a diferencia del RMSE crudo.
"""

from stats import distributions_distance


def run():
    raise NotImplementedError(
        "4c es una decision de diseno: la metrica es stats.distributions_distance. "
        "Falta escribir la justificacion en el informe."
    )
