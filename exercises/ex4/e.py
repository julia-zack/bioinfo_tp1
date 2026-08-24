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


def run():
    raise NotImplementedError("4e todavia no esta implementado")
