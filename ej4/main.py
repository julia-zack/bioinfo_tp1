
# 4a) Compare distribucion de aminoacidos de secuencia de proteinas al azar vs. secuencias de proteinas reales
# Diseño en dos pasos:
# Paso 1
# Elegimos las proteinas reales que vamos a usar: E. coli y humano
# Cuando se estabiliza la distribucion?
# -> definimos N como la cantidad de residuos a partir de la cual la estimacion de las frecuencias deja de cambiar
#
# Paso 2 Comparamos las tres distribuciones (azar, E. coli, humano)
# a tamaños de muestra << N, ~ N y >> N, para ver graficamente que por debajo de N la
# estimacion es poco confiable y por encima es estable.

# Correr desde esta carpeta:  python3 4-a.py
# Genera:  4a_convergencia.png  y  4a_regimenes.png
# Las secuencias reales ya vienen descargadas en data/. Si borras esos archivos,
# la primera corrida las vuelve a bajar de NCBI y las cachea de nuevo ahi.
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # para importar el paquete utils (dentro de ej4/)

from utils.utils import (
    AA,
    UMBRAL,
    freqs,
    random_residues,
    concatenate_til_n,
    curva_convergencia,
    encontrar_N,
    get_real_sequences,
)

# organismos reales que vamos a comparar (nombre -> (nombre NCBI, archivo cache))
REAL_ORGANISMS = {
    "E. coli": ("Escherichia coli", os.path.join(HERE, "utils", "data", "sequences_ecoli.json")),
    "Humano":  ("Homo sapiens",     os.path.join(HERE, "utils", "data", "sequences_human.json")),
}
COLOURS = {"Azar": "#c1666b", "E. coli": "#4a6fa5", "Humano": "#6a9955"}


def download_seqs():
    """Carga (o descarga y cachea) las secuencias reales de cada organismo."""
    real_seqs = {}
    for name, (org, cache) in REAL_ORGANISMS.items():
        real_seqs[name] = get_real_sequences(org, cache)
    return real_seqs


def calculate_data_length(real_seqs):
    """Cantidad de residuos de la fuente real MAS CHICA. Igualamos todas las
    fuentes a este total (TOT), asi la comparacion es justa: todas se miden con
    la misma cantidad de datos."""
    totales_reales = []
    for seqs in real_seqs.values():
        total_residuos = 0
        for proteina in seqs:
            total_residuos += len(proteina)
        totales_reales.append(total_residuos)
    return min(totales_reales)


def build_sources(real_seqs, length):
    """Arma, para cada fuente, una cadena de EXACTAMENTE `length` residuos:
    una secuencia al azar, y las reales recortadas a `length`."""
    sources = {"Azar": random_residues(length)}
    for name, seqs in real_seqs.items():
        secuencia_completa = "".join(seqs)
        sources[name] = secuencia_completa[:length]
    return sources


def analyze_convergence(sources):
    """Para cada fuente calcula: su curva de convergencia, su N, y su
    distribucion final (usando todos los residuos)."""
    curves = {}
    N_by_source = {}
    final_distributions = {}
    for name, res in sources.items():
        grid, ds = curva_convergencia(res)
        curves[name] = (grid, ds)
        N_by_source[name] = encontrar_N(grid, ds)
        final_distributions[name] = freqs(res)
    return curves, N_by_source, final_distributions


def print_N_by_source(N_by_source):
    print("\nTamano de estabilizacion N (residuos) por fuente:")
    for name, N in N_by_source.items():
        print(f"  {name}: N = {N}")


def representative_N(curves):
    """N donde se estabiliza el PROMEDIO de las curvas de todas las fuentes.
    Las tres convergen juntas, asi que un solo N resume "cuando alcanza para
    todas" (y evita el ruido de que el N por fuente salte entre corridas)."""
    todas_las_curvas = []
    for grid, ds in curves.values():
        todas_las_curvas.append(ds)
    curva_media = list(np.mean(todas_las_curvas, axis=0))
    grid_comun = curves["Azar"][0]   # todas comparten la misma grilla
    return encontrar_N(grid_comun, curva_media)


def mark_N_on_xaxis(ax, N):
    """Agrega N como un tick mas del eje X, con el mismo formato que los 10^k
    pero resaltado (negrita + cuarto color) y una sola vez."""
    # 1) armamos los ticks "decada" (10^2, 10^3, ...) que entran en el rango
    lo, hi = ax.get_xlim()
    exponente_min = int(np.ceil(np.log10(lo)))
    exponente_max = int(np.floor(np.log10(hi)))
    decadas = []
    for k in range(exponente_min, exponente_max + 1):
        decadas.append(10 ** k)
    # 2) sumamos N a la lista de ticks (sin duplicar si cae en una decada)
    ticks = sorted(set(decadas + [N]))
    ax.set_xticks(ticks)
    # 3) etiqueta: N va como numero comun; las decadas como 10^k
    etiquetas = []
    for t in ticks:
        if t == N:
            etiquetas.append(str(int(N)))
        else:
            exponente = int(round(np.log10(t)))
            etiquetas.append(r'$10^{%d}$' % exponente)
    ax.set_xticklabels(etiquetas)
    # 4) resaltamos la etiqueta de N (negrita + color)
    for lbl, t in zip(ax.get_xticklabels(), ticks):
        if t == N:
            lbl.set_color('#d1462f')       # cuarto color, para que resalte
            lbl.set_fontweight('bold')


def plot_convergence(curves, N):
    """Figura 1: una curva de convergencia por fuente, con N marcado."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, (grid, ds) in curves.items():
        ax.plot(grid, ds, marker='o', color=COLOURS[name], label=name)
    ax.axhline(UMBRAL, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0,
               label=f'umbral = {UMBRAL}')
    ax.axvline(N, color='#d1462f', linestyle=':', linewidth=1.3)
    ax.set_xscale('log')
    ax.set_xlabel('Residuos analizados')
    ax.set_ylabel('Distancia a la distribucion final')
    ax.set_title('4a - Paso 1: ¿cuando se estabiliza la distribucion de aa?')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)

    mark_N_on_xaxis(ax, N)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, '4a_convergencia.png'), dpi=150)


def choose_regimes(N, data_length):
    """Tres tamaños de muestra para comparar: mucho menor que N, cerca de N,
    y toda la data disponible (mucho mayor que N)."""
    R_small = 1000                # mucho menor que N: todavia ruidoso (probar 500 para mas ruido)
    R_med = min(N, data_length)   # cerca de N
    R_large = data_length         # toda la data disponible (>> N)
    return [
        (f'<< N  (R = {R_small})', R_small),
        (f'~ N   (R = {R_med})',   R_med),
        (f'>> N  (R = {R_large})', R_large),
    ]


def plot_regimes(real_seqs, final_distributions, N, data_length):
    """Figura 2: distribucion de aa de las tres fuentes a tres tamaños de
    muestra (<< N, ~ N, >> N)."""
    regimenes = choose_regimes(N, data_length)

    # ordenamos los aa por su frecuencia en E. coli, para que el grafico se lea mejor
    orden = sorted(AA, key=lambda a: final_distributions["E. coli"].get(a, 0), reverse=True)
    posiciones_base = np.arange(len(orden))   # una posicion en x por aminoacido
    ancho_barra = 0.27

    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for ax, (etiqueta_regimen, R) in zip(axes, regimenes):
        # distribucion de aa de cada fuente usando R residuos
        muestras = {
            "Azar":    freqs(random_residues(R)),
            "E. coli": freqs(concatenate_til_n(real_seqs["E. coli"], R)),
            "Humano":  freqs(concatenate_til_n(real_seqs["Humano"], R)),
        }
        # una tanda de barras por fuente, corridas a izquierda/centro/derecha
        for j, nombre in enumerate(["Azar", "E. coli", "Humano"]):
            distribucion = muestras[nombre]
            alturas = []
            for a in orden:
                alturas.append(distribucion.get(a, 0))
            posiciones = posiciones_base + (j - 1) * ancho_barra
            ax.bar(posiciones, alturas, ancho_barra, label=nombre, color=COLOURS[nombre])

        ax.axhline(1 / len(AA), color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
        ax.set_title(f'Muestra {etiqueta_regimen} residuos')
        ax.set_ylabel('Frecuencia')
        ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    axes[0].legend(frameon=False, ncol=3)
    axes[-1].set_xticks(posiciones_base)
    axes[-1].set_xticklabels(orden)
    axes[-1].set_xlabel('Aminoacido (ordenados por frecuencia en E. coli)')
    fig.suptitle('4a - Paso 2: distribucion de aa a distintos tamanos de muestra')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, '4a_regimenes.png'), dpi=150)


def main():
    real_seqs = download_seqs()
    data_length = calculate_data_length(real_seqs)
    sources = build_sources(real_seqs, data_length)

    # Paso 1: cuando se estabiliza la distribucion de aa
    curves, N_by_source, final_distributions = analyze_convergence(sources)
    print_N_by_source(N_by_source)
    N = representative_N(curves)
    print(f"\nN representativo (promedio de las tres fuentes): {N}")
    plot_convergence(curves, N)

    # Paso 2: comparar las distribuciones a << N, ~ N, >> N
    plot_regimes(real_seqs, final_distributions, N, data_length)

    plt.show()
    print("\nGuardados: 4a_convergencia.png y 4a_regimenes.png")


if __name__ == "__main__":
    main()
