# TP 1 — Análisis de secuencias

---

## 4a

> **4a)** Compare Distribución de aminoácidos de secuencia de proteínas al azar vs. secuencias de proteínas reales.

### Método

Se comparan cuatro fuentes con **exactamente 79.453 residuos** cada una (el total de la fuente real más chica, *E. coli*), para que ninguna tenga ventaja por cantidad de datos, más un perfil combinado que sirve de referencia:

| Fuente | Qué es |
|---|---|
| **Random DNA** | Control al azar: ADN aleatorio traducido, descartando los stops. Los residuos no tienen información, pero respetan el código genético. En la figura 1 es una muestra de 79.453 residuos; en la figura 2, donde funciona como referencia, se usan sus frecuencias **exactas**: la probabilidad de cada codón bajo ADN uniforme, agrupada por aminoácido, que con GC = 0,5 es (número de codones)/61. |
| **E. coli**, **Yeast**, **Human** | Proteínas reales (300 por organismo, RefSeq, vía NCBI). |
| **Natural** | Los tres organismos combinados en partes iguales (79.453 residuos de cada uno, 238.359 en total). Es la distribución de referencia que usa 4e. Al ser una referencia y no una fuente más en la comparación, usa todos los residuos disponibles. |

Dos decisiones que conviene explicitar:

- **El control es ADN al azar traducido, no aminoácidos uniformes.** Un control uniforme atribuiría a la biología diferencias que en realidad explica la degeneración del código (Leu tiene 6 codones, Met tiene 1). El caso uniforme igual está en el gráfico: es la línea punteada en 1/20 = 0,05.
- **Las listas de proteínas se mezclan antes de recortar** (`SAMPLE_SEED`). Sin mezclar, quedarse con los primeros 79.453 residuos significaba quedarse con las primeras ~116 proteínas que devolvió NCBI, y eso da una distribución 0,065 lejos de la real, cuando una muestra al azar del mismo tamaño queda a 0,008.

### Resultados

![Distribución de aminoácidos](exercises/ex4/aa_distribution.png)

![Perfil natural vs control](exercises/ex4/aa_natural_vs_random.png)

Las distancias entre distribuciones (distancia de variación total, ver 4c):

| | vs Natural |
|---|---|
| Human | 0,043 |
| E. coli | 0,062 |
| Yeast | 0,069 |
| **Random DNA** | **0,125** |

### Discusión

**Las secuencias reales no son uniformes, pero el azar tampoco.** El control de ADN al azar ya sale muy lejos de 1/20: Leu 0,098 y Met 0,016, porque la frecuencia esperada de cada aminoácido es (número de codones)/61. Comparar contra un control uniforme habría exagerado el efecto de la biología.

**La diferencia real se separa en dos partes.** Tomando el control como referencia:

| aa | Codones | Random DNA | Natural | Rango organismos |
|---|---|---|---|---|
| L | 6 | 0,098 | 0,096 | 0,092–0,098 |
| S | 6 | 0,099 | 0,079 | 0,064–0,088 |
| R | 6 | 0,100 | 0,057 | 0,048–0,063 |
| E | 2 | 0,032 | 0,064 | 0,060–0,069 |
| K | 2 | 0,033 | 0,057 | 0,042–0,073 |
| C | 2 | 0,033 | 0,015 | 0,011–0,019 |

- **Leucina**: el código explica todo. Es el aminoácido más abundante tanto al azar como en proteínas reales, simplemente porque tiene 6 codones.
- **Arginina**: el código predice 0,100 pero los organismos usan la mitad. Acá sí hay selección.
- **Cisteína**: predicha en 0,033, observada en 0,015. Es reactiva y forma puentes disulfuro, así que su uso está restringido.
- **Glutámico y lisina**: al revés, los organismos los usan casi el doble de lo que predice el código.

**Los organismos se diferencian entre sí casi tanto como del azar.** *E. coli*–Yeast da 0,119 y Yeast–Human 0,100, contra 0,125 de Natural al control. Es decir: "natural" no es un blanco angosto sino una nube ancha. El perfil combinado queda a 0,043–0,069 de cada organismo y a 0,125 del control, un margen de apenas ~1,8×.

Esto es una limitación para 4e: la distribución de aminoácidos sola no alcanza para separar codificante de no codificante con confianza, y el largo del ORF va a tener que aportar la mayor parte del poder de discriminación.

**Qué aminoácidos sirven para discriminar.** En `aa_natural_vs_random.png` las barras verticales muestran el rango entre los tres organismos. Cuando ese rango es chico y la separación al control es grande, el residuo discrimina; cuando el rango se superpone con el control, no aporta nada.

- Sirven: **E, N, K, R, C** — separación grande y consistente en los tres organismos.
- No sirven: **A, G, V, I, D, P, Y, W** — la diferencia entre organismos es mayor que la diferencia contra el azar.

### Cómo reproducirlo

```
python3 main.py 4a
```

Código en `exercises/ex4/a.py`. Las secuencias quedan cacheadas en `data/`; si se borran, la primera corrida las vuelve a bajar de NCBI.

---

## 4c

> **4c)** Elija una métrica que permita comparar dos distribuciones. (p.ej., RMSE).

Se eligió la **distancia de variación total** (TV), implementada en `stats.distributions_distance()`:

```
TV(p, q) = 0,5 · Σ |p(a) − q(a)|
```

Se prefirió sobre RMSE porque:

- Está acotada entre 0 y 1, así que el número se interpreta solo.
- Se lee directamente como "qué fracción de la masa de probabilidad hay que mover para convertir una distribución en la otra".
- No depende de la cantidad de bins, a diferencia del RMSE crudo.
- Permite restringirla a un subconjunto de aminoácidos sumando menos términos, que es lo que se explota en 4d.

---

## 4d

> **4d)** Escriba un código que dadas dos distribuciones las compare y obtenga la métrica correspondiente. Utilícela para: i) ver cómo evolucionan las distribuciones obtenidas en 4a al aumentar el tamaño de la muestra (evaluar cuándo la distribución deja de cambiar). ii) estudiar de manera sistemática las diferencias entre las distribuciones obtenidas en 4a para diferentes muestras (comparar distribuciones de secuencias naturales y al azar).

La parte **(i)** es el análisis de convergencia de 4b: ahí la métrica se aplica a cada fuente contra su propia distribución final. Esta sección es la parte **(ii)**.

### Distancias entre todas las fuentes

```
              Random DNA     E. coli       Yeast       Human     Natural
  Random DNA       0.000       0.135       0.166       0.110       0.125
     E. coli       0.135       0.000       0.119       0.086       0.062
       Yeast       0.166       0.119       0.000       0.100       0.069
       Human       0.110       0.086       0.100       0.000       0.043
     Natural       0.125       0.062       0.069       0.043       0.000
```

Lo importante no son los valores sino su escala relativa: **la distancia entre organismos (0,086–0,119) es del mismo orden que la distancia de cada organismo al azar (0,110–0,166)**. *Yeast* está a 0,119 de *E. coli* y a 0,166 del control: apenas 1,4× más lejos del azar que de otro organismo real.

Es decir, "natural" no es un punto sino una nube ancha, y cualquier detector que use un único perfil de referencia hereda esa limitación.

### Qué tan bien separa la métrica una secuencia individual

La pregunta que importa para 4e no es cuánto se diferencian dos corpus grandes, sino si la métrica puede decidir sobre **una** secuencia. Se midió con AUC: se toman fragmentos de proteínas reales y controles al azar del mismo largo, se calcula la distancia de cada uno al perfil Natural, y se mide la probabilidad de que el fragmento real quede más cerca. 0,5 es azar, 1,0 es separación perfecta.

![Discriminación según el largo](exercises/ex4/orf_discrimination.png)

```
  length      all 20   E R K N C         R C
      16       0.458       0.542       0.621
      30       0.470       0.569       0.617
      60       0.500       0.633       0.699
     120       0.509       0.714       0.773
     250       0.546       0.794       0.845
     400       0.612       0.872       0.909
```

AUC de cada aminoácido por separado, con fragmentos de 250 residuos:

```
  R: 0.835   E: 0.816   C: 0.783   D: 0.689   K: 0.594   N: 0.589   S: 0.551   H: 0.526
```

### Discusión

**Usar los 20 aminoácidos es casi inútil.** Se queda en 0,46–0,55 hasta los 250 residuos: prácticamente indistinguible de tirar una moneda. La razón es de ruido, no de biología: en un fragmento de 250 residuos cada frecuencia tiene un error de muestreo de ±0,013, y sumar veinte términos de error tapa las cinco o seis diferencias que sí son reales.

**Un subconjunto chico funciona mucho mejor.** Con sólo R y C la AUC sube a 0,845 en 250 residuos y 0,909 en 400. R y C son justamente los casos donde el código genético predice frecuencia alta y la biología entrega baja (R: 0,100 esperado contra 0,057 observado; C: 0,033 contra 0,015), así que la discrepancia es grande y consistente en los tres organismos.

**La composición no sirve en ORFs cortos.** A 16 residuos —la mediana de los ORFs medidos en 3c— hasta el mejor subconjunto da 0,62, y los 20 juntos dan 0,458, *peor que azar*. La señal recién se vuelve usable arriba de ~120 residuos.

### Consecuencia para 4e

1. **El largo es la señal principal.** Es lo único que discrimina en ORFs cortos, y es donde 3c ya mostró la diferencia más fuerte (máximo real 1121 contra 231 del control).
2. **La composición es señal secundaria**, y sólo aporta arriba de ~100 residuos.
3. **Conviene usar pocos aminoácidos, no los 20.** Incluir los que no discriminan degrada activamente el resultado.

Estos tres puntos justifican la forma del score de 4e: probabilidad base por largo, corregida por composición sólo cuando el ORF es lo bastante largo.

### Advertencia metodológica

El subconjunto R/C se eligió mirando los mismos datos con los que después se lo evaluó, así que 0,909 es optimista. La versión rigurosa elige los aminoácidos usando dos organismos y evalúa sobre el tercero, dejado afuera.

### Cómo reproducirlo

```
python3 main.py 4d
```

Código en `exercises/ex4/d.py`.
