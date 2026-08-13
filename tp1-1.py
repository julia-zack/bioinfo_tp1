# Valina (Val, V)
# Leucina (Leu, L)
# Treonina (Thr, T)
# Lisina (Lys, K)
# Triptófano (Trp, W)
# Histidina (His, H)
# Fenilalanina (Phe, F)
# Isoleucina (Ile, I)
# Arginina (Arg, R)
# Metionina (Met, M)	Alanina (Ala, A)
# Prolina (Pro, P)
# Glicina (Gly, G)
# Serina (Ser, S)
# Cisteína (Cys, C)
# Asparagina (Asn, N)
# Glutamina (Gln, Q)
# Tirosina (Tyr, Y)
# Ácido aspártico (Asp, D)
# Ácido glutámico (Glu, E)

import random
import numpy as np

aa = ['V','L','T','K','W','H','F','I','R','A','P','G','S','C','N','Q','Y','D','E']
aa_weights = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.04]

weighted_aa = np.random.choice(aa, size=len(aa), p=aa_weights)



def generate_sequence(length):
    return ''.join(random.choice(aa) for _ in range(length))


input_length = int(input("Enter the length of the amino acid sequence: "))
print (generate_sequence(input_length))