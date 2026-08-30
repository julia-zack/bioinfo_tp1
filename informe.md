# TP 1 — Análisis de secuencias

---

## Datos: de dónde salen las secuencias reales

Todo el objetivo 4 compara secuencias al azar contra secuencias reales, así que **qué se entiende por "secuencia real" define los resultados**. Las proteínas se bajan de **Swiss-Prot**, no de RefSeq, y la diferencia no es cosmética.

### Qué es Swiss-Prot

UniProt tiene dos mitades:

| | Swiss-Prot | TrEMBL |
|---|---|---|
| Curado por | personas | software |
| Tamaño | ~570.000 entradas (todas las especies) | ~250 millones |
| Se lo llama | *reviewed* | *unreviewed* |

**Swiss-Prot es la mitad revisada a mano**: un curador lee la bibliografía y escribe la entrada. La regla que más importa acá es que hay **una entrada por gen por organismo**; las isoformas viven como anotaciones adentro de la entrada, no como registros separados.

Se accede con el mismo cliente Entrez que ya usaba el código, porque NCBI espeja Swiss-Prot: alcanza con cambiar el filtro de la query (`swissprot[filter]` en vez de `refseq[filter]`). Es además el formato que pide el objetivo 3 del TP, el que parsean `SeqIO.read(handle, "swiss")` y `ExPASy.get_sprot_raw()`.

### Por qué no RefSeq

RefSeq es **genómico**: un pipeline de anotación automática emite un registro de proteína por cada gen predicho de cada genoma secuenciado. Eso trae tres problemas para un estudio de composición de aminoácidos, los tres verificados sobre las queries que usaba el código:

1. **Redundancia por isoformas.** Los primeros hits de humano en RefSeq eran cuatro isoformas de *maestro heat-like repeat-containing protein family member 1* y cinco de *serine/threonine-protein phosphatase 2A activator*. Las 400 "proteínas" humanas eran más bien ~100 genes con duplicados. Secuencias casi idénticas no aportan información nueva.
2. **Entradas no verificadas.** Todo ORF que llame el predictor se convierte en registro. La muestra de *E. coli* estaba dominada por `hypothetical protein`, dominios `DUF####` de función desconocida, proteínas de fago y secuencias marcadas `partial`.
3. **Duplicación por cepa.** `"Escherichia coli"[Organism] AND refseq[filter]` matchea **6.666.151** proteínas — una por cada cepa anotada — contra 23.300 en Swiss-Prot.

El contraste en las tres queries:

| Organismo | RefSeq | Swiss-Prot |
|---|---|---|
| *E. coli* (todas las cepas) | 6.666.151 | 23.300 |
| *E. coli* K-12 | — | **6.074** |
| *S. cerevisiae* S288C | 6.029 | 7.923 |
| *H. sapiens* | 197.929 | 20.616 |

**Para *E. coli* se usa K-12, no la especie entera.** `"Escherichia coli"[Organism]` matchea todas las cepas secuenciadas, y aun dentro de Swiss-Prot eso reintroduce por la ventana la redundancia por cepa que se acababa de sacar por la puerta. K-12 es el proteoma de *un* organismo. (Swiss-Prot asigna las entradas a nivel de cepa y no de subcepa, así que `"...str. K-12 substr. MG1655"[Organism]` no matchea nada.) Levadura y humano ya eran de un solo organismo.

Las 20.616 entradas humanas de Swiss-Prot son la comprobación de que efectivamente hay una por gen: el proteoma humano revisado ronda los 20.000 genes. Las 197.929 de RefSeq son ese mismo proteoma multiplicado por variantes de transcripto y predicciones.

### Cómo se elige la muestra

Cambiar de base de datos arregla la redundancia y las entradas basura, pero es un problema **independiente** de cómo se elige qué bajar, y ese segundo problema era real.

`esearch` devuelve los ids en el orden por defecto de NCBI, que es por recencia de depósito. Tomar los primeros *n* es entonces quedarse con lo último que se depositó, no con una muestra del organismo. Se nota mirando la cabeza de las listas: los primeros hits de levadura son `Uncharacterized protein YNL155C-A` y parecidos, y los de humano son péptidos chicos tipo `Small humanin-like peptide 6`.

El sesgo se midió. Bajando del mismo pool de 6.074 entradas de *E. coli* K-12, los primeros 400 ids contra 400 sorteados al azar:

| Muestra | Largo medio |
|---|---|
| Primeros 400 (cabeza de la lista) | 256 |
| 400 al azar | **304** |

Tomar la cabeza sesga hacia proteínas 19% más cortas, y es la muestra al azar la que se acerca al largo medio conocido del proteoma de K-12 (~316). El sesgo además apuntaba en direcciones distintas según el organismo — humano daba 746 por cabeza contra 513 al azar, o sea *al revés* —, lo cual es peor que un sesgo constante: contamina las comparaciones **entre** organismos, no sólo los valores absolutos.

### La solución: bajar el proteoma completo

Para estos tres organismos el conjunto revisado no es una muestra del proteoma, **es el proteoma**. Así que no se muestrea: se baja todo.

| Organismo | Proteínas | Residuos | Largo medio | Largo medio conocido |
|---|---|---|---|---|
| *E. coli* K-12 | 6.062 | 1.848.964 | 305 | ~316 |
| *S. cerevisiae* | 7.888 | 3.544.698 | 449 | ~450 |
| *H. sapiens* | 20.591 | 11.454.685 | 556 | ~560 |

Los tres largos medios caen sobre el valor conocido. Antes ninguno lo hacía. Como no hay selección, no hay sesgo de selección posible — la pregunta de cómo muestrear desaparece en lugar de resolverse.

(Los totales quedan un poco abajo del número de entradas de la query porque `clean_sequences()` descarta las que tienen caracteres no estándar: selenocisteína, posiciones ambiguas.)

Esto además es lo que hace **medible** el ejercicio 4b. Preguntar "¿cuántas proteínas alcanzan?" contra una referencia construida con las mismas 400 proteínas contesta sobre todo "bajaste 400"; contra el proteoma entero, una muestra de 400 es el 2% de la referencia y la respuesta pasa a ser sobre el organismo. El detalle está en 4b.

### Dos detalles de implementación

- **La paginación de `esearch` no es una optimización.** `esearch` devuelve como máximo 10.000 ids por request y el proteoma humano revisado tiene 20.616 entradas, así que una sola llamada devolvería la mitad más reciente. `fetch_id_pool()` pagina hasta juntar la lista completa. Sin eso, incluso un sorteo al azar sería un sorteo dentro de un pool ya sesgado.
- **`fetch_proteins()` acepta `n` para tomar una muestra al azar** en vez de todo, con `seed` fijo para que sea reproducible. No se usa en 4a, pero es lo que permite los experimentos de 4b que varían el tamaño del corpus.

Código en `ncbi.py`; las secuencias quedan cacheadas en `data/swissprot_*.json`. El nombre del archivo lleva la base de datos justamente para que un cache viejo de RefSeq no se reutilice sin querer.

---

## 4a

> **4a)** Compare Distribución de aminoácidos de secuencia de proteínas al azar vs. secuencias de proteínas reales.

### Método

Se comparan cuatro fuentes con **exactamente 1.848.964 residuos** cada una (el total de la fuente real más chica, el proteoma de *E. coli* K-12), para que ninguna tenga ventaja por cantidad de datos, más un perfil combinado que sirve de referencia:

| Fuente | Qué es |
|---|---|
| **Random DNA** | Control al azar: ADN aleatorio traducido, descartando los stops. Los residuos no tienen información, pero respetan el código genético. En la figura 1 es una muestra de 1.848.964 residuos; en la figura 2, donde funciona como referencia, se usan sus frecuencias **exactas**: la probabilidad de cada codón bajo ADN uniforme, agrupada por aminoácido, que con GC = 0,5 es (número de codones)/61. |
| **E. coli**, **Yeast**, **Human** | Proteomas revisados completos de Swiss-Prot (6.062 / 7.888 / 20.591 proteínas; ver la sección *Datos*). |
| **Natural** | Los tres organismos combinados en partes iguales (1.848.964 residuos de cada uno, 5.546.892 en total). Es la distribución de referencia que usa 4e. |

Dos decisiones que conviene explicitar:

- **El control es ADN al azar traducido, no aminoácidos uniformes.** Un control uniforme atribuiría a la biología diferencias que en realidad explica la degeneración del código (Leu tiene 6 codones, Met tiene 1). El caso uniforme igual está en el gráfico: es la línea punteada en 1/20 = 0,05.
- **Las listas de proteínas se mezclan antes de recortar** (`SAMPLE_SEED`). Levadura y humano tienen más residuos que *E. coli*, así que se los recorta, y sin mezclar ese recorte se quedaría con las proteínas que NCBI devuelve primero, que son las depositadas más recientemente. Medido contra la distribución del proteoma completo: humano sin mezclar queda a 0,019 y mezclado a 0,004. *E. coli* da 0,000 en los dos casos porque es la fuente más chica y no se recorta nada.

### Resultados

![Distribución de aminoácidos](exercises/ex4/aa_distribution.png)

![Perfil natural vs control](exercises/ex4/aa_natural_vs_random.png)

Las distancias entre distribuciones (distancia de variación total, ver 4c):

| | vs Natural |
|---|---|
| Human | 0,045 |
| E. coli | 0,063 |
| Yeast | 0,072 |
| **Random DNA** | **0,131** |

### Discusión

**Las secuencias reales no son uniformes, pero el azar tampoco.** El control de ADN al azar ya sale muy lejos de 1/20: Leu 0,098 y Met 0,016, porque la frecuencia esperada de cada aminoácido es (número de codones)/61. Comparar contra un control uniforme habría exagerado el efecto de la biología.

**La diferencia real se separa en dos partes.** Tomando el control como referencia:

| aa | Codones | Random DNA | Natural | Rango organismos |
|---|---|---|---|---|
| L | 6 | 0,098 | 0,100 | 0,095–0,106 |
| S | 6 | 0,098 | 0,077 | 0,057–0,092 |
| R | 6 | 0,098 | 0,052 | 0,045–0,057 |
| E | 2 | 0,033 | 0,065 | 0,059–0,071 |
| K | 2 | 0,033 | 0,058 | 0,044–0,073 |
| D | 2 | 0,033 | 0,052 | 0,047–0,058 |
| C | 2 | 0,033 | 0,016 | 0,012–0,022 |

- **Leucina**: el código explica todo. Es el aminoácido más abundante tanto al azar (0,098) como en proteínas reales (0,100), simplemente porque tiene 6 codones.
- **Arginina**: el código predice 0,098 pero los organismos usan casi la mitad, 0,052. Acá sí hay selección, y es la diferencia más grande de las veinte.
- **Cisteína**: predicha en 0,033, observada en 0,016. Es reactiva y forma puentes disulfuro, así que su uso está restringido.
- **Glutámico, aspártico y lisina**: al revés, los organismos los usan cerca del doble de lo que predice el código. Los tres están cargados.

**Los organismos se diferencian entre sí casi tanto como del azar.** *E. coli*–Yeast da 0,122 y Yeast–Human 0,101, contra 0,131 de Natural al control. Es decir: "natural" no es un blanco angosto sino una nube ancha. El perfil combinado queda a 0,045–0,072 de cada organismo y a 0,131 del control, un margen de apenas ~1,8×.

Esto es una limitación para 4e: la distribución de aminoácidos sola no alcanza para separar codificante de no codificante con confianza, y el largo del ORF va a tener que aportar la mayor parte del poder de discriminación.

**Qué aminoácidos sirven para discriminar.** En `aa_natural_vs_random.png` las barras verticales muestran el rango entre los tres organismos. El criterio es comparar dos cantidades por aminoácido: la **separación** contra el control (|natural − azar|) y el **rango** entre organismos. Sirve el que tiene separación grande y rango chico; el que varía más entre organismos que contra el azar no aporta nada.

| aa | Separación | Rango | ¿Sirve? |
|---|---|---|---|
| R | 0,046 | 0,012 | sí, con margen 4× |
| E | 0,032 | 0,012 | sí |
| D | 0,020 | 0,011 | sí |
| C | 0,017 | 0,011 | sí |
| K | 0,025 | 0,029 | no: varía más entre organismos |
| S | 0,021 | 0,035 | no |
| A | 0,008 | 0,041 | no |

- Sirven: **R, E, D, C** — separación grande y consistente en los tres organismos. (Q, T y H también cumplen el criterio, pero con separaciones de 0,009–0,011, demasiado chicas para aportar.)
- No sirven: **A, G, V, I, L, M, W, Y, F, P** por separación nula, y **K, S, N** porque, aunque su separación es grande, la diferencia entre organismos lo es más.

El caso de **K** es el más instructivo: 0,033 al azar contra 0,058 natural parece una señal fuerte, pero los organismos van de 0,044 a 0,073 entre sí. Un ORF rico en lisina puede ser perfectamente humano o perfectamente levadura; no dice si es real. Este ranking se confirma después de manera independiente en 4d, midiendo AUC por aminoácido: R 0,905, E 0,722, C 0,691, y recién después K 0,629.

### Cómo reproducirlo

```
python3 main.py 4a
```

Código en `exercises/ex4/a.py`. Las secuencias quedan cacheadas en `data/swissprot_*.json`; si se borran, la primera corrida vuelve a bajar los tres proteomas (unos dos minutos).

---

## 4b

> **4b)** Analice cómo cambian las distribuciones al aumentar el tamaño de la secuencia "al azar" analizada, y al incrementar el número de secuencias reales analizadas. ¿Cuándo es suficiente?

### Método

La consigna nombra **dos ejes distintos**, y no son el mismo: se puede crecer en cantidad de **residuos** o en cantidad de **proteínas enteras**. Se miden por separado.

En ambos casos la pregunta "¿cuándo es suficiente?" se responde igual: se compara la distribución obtenida con una muestra parcial contra la distribución final de esa misma fuente (usando *todos* sus datos), usando una métrica (la que luego usaremos en 4c), y se busca el punto donde la distancia cae por debajo de un umbral de 0,01.

- **Eje 1 — residuos.** Aplica a las cuatro fuentes, incluido el control al azar. Llamamos **N** al número de residuos donde se cruza el umbral.
- **Eje 2 — número de secuencias.** Aplica sólo a los organismos reales: el control al azar se genera como una única cadena, así que "número de secuencias" no está definido para él. Llamamos **K** al número de proteínas enteras donde se cruza el umbral.

Cada curva se promedia sobre 5 barajadas del orden (`REPS`), para que no dependa de qué proteínas vinieron primero.

### Resultados

![Estabilización de la distribución](exercises/ex4/aa_stabilized.png)

**Eje 1 (residuos).** Todas las fuentes se estabilizan en el mismo orden de magnitud, ~3·10⁴ residuos:

| Fuente | N (residuos), mediana de 6 corridas |
|---|---|
| Random DNA | 31.400 |
| E. coli | 33.500 |
| Yeast | 31.600 |
| Human | 31.300 |

**Eje 2 (proteínas enteras).** Alcanza con unas **240–425 proteínas**, pero eso equivale a muchos más residuos que el N del eje 1:

| Organismo | K (proteínas) | de un proteoma de | Largo medio | K en residuos |
|---|---|---|---|---|
| E. coli | ~240 | 6.062 | 305 | ~73.000 |
| Yeast | ~252 | 7.888 | 449 | ~113.000 |
| Human | ~425 | 20.591 | 556 | ~236.000 |

Una vez pasado ese punto, la figura de regímenes muestra que la distribución efectivamente deja de moverse:

![Distribuciones a distintos tamaños de muestra](exercises/ex4/aa_sample_sizes.png)

Con R = 1000 (<< N) las barras son ruido: levadura aparece con K en 0,096 cuando su valor real es 0,073, *E. coli* con L en 0,119 contra 0,106, humano con P en 0,075 contra 0,064. Ninguno se sostiene. Con R ≈ N y R >> N los dos paneles son prácticamente idénticos — comparar L, S, R, E o K entre el segundo y el tercero: las diferencias quedan en el tercer decimal. Ahí está la respuesta a "cuándo es suficiente": el panel del medio ya contiene toda la información que aporta el de abajo, usando el **1,7%** de los datos (31.878 residuos contra 1.848.964).

### Discusión

**Una proteína entera vale mucho menos que su cantidad de residuos.** Es el resultado más interesante de los dos ejes juntos. Con residuos sueltos alcanza con ~31.000; sumando proteínas enteras hacen falta entre 2 y 8 veces más residuos para llegar al mismo umbral. La razón es que los residuos dentro de una proteína **no son independientes**: una proteína de membrana es rica en hidrofóbicos, una ribosomal es rica en K y R. Cada proteína aporta una composición propia, así que agregar proteínas de a bloques agrega menos información que agregar la misma cantidad de residuos al azar.

**El efecto crece con el largo medio.** La penalización por correlación se lee dividiendo los residuos que hacen falta por cada eje:

| Organismo | Largo medio | K en residuos | N (residuos) | Penalización |
|---|---|---|---|---|
| E. coli | 305 | ~73.000 | ~33.500 | 2,2× |
| Yeast | 449 | ~113.000 | ~31.600 | 3,6× |
| Human | 556 | ~236.000 | ~31.300 | 7,5× |

La penalización crece monótonamente con el largo medio de las proteínas del organismo. Cuanto más larga la proteína, más correlacionados los residuos que aporta y peor rinde cada uno: en humano hacen falta 7,5 veces más residuos si vienen empaquetados en proteínas que si vinieran sueltos.

**Consecuencia práctica:** los proteomas completos que usa 4a (6.062 / 7.888 / 20.591 proteínas) están entre 14 y 48 veces por encima de K. El perfil de referencia que hereda 4e no tiene problema de tamaño de muestra. Pero un trabajo que bajara "unas 100 proteínas por organismo" estaría por debajo del umbral en los tres casos, y en humano por un factor de 4.

### Advertencias metodológicas

- **N y K son ruidosos.** Como el orden se baraja sin semilla fija, cambian de corrida en corrida: sobre 6 corridas, N va de 17.800 a 35.400 y K de 206 a 666. Los valores de las tablas son medianas; lo que se sostiene es el orden de magnitud (N ~ 3·10⁴ residuos, K ~ 2·10² a 4·10² proteínas), no la cifra exacta. El ruido viene de que la curva es casi plana donde se la corta: cerca del umbral cae 0,006 a lo largo de una celda entera de la grilla, así que un desvío vertical de 0,002 corre el cruce por casi toda la celda.
- **La referencia sigue conteniendo a la muestra**, aunque ahora casi no importe. Cada curva se compara contra la distribución final de su propia fuente. Con el proteoma completo, una muestra de K = 425 proteínas es el 2% de las 20.591 de la referencia, así que la contaminación es despreciable. Con el corpus de 400 proteínas que se usaba antes era del 50%, y eso invalidaba la medición — ver abajo.

### Por qué estos números son creíbles, y antes no lo eran

Hay un test que distingue una medición real de un artefacto: **variar la cantidad de datos disponibles**. Si K mide una propiedad del organismo, no debería depender de cuántas proteínas bajamos.

Con un corpus de 400 proteínas, K salía siempre alrededor de la mitad del corpus, cualquiera fuera el corpus:

| Corpus | K | K/corpus |
|---|---|---|
| 25 | 22 | 0,88 |
| 50 | 41 | 0,82 |
| 100 | 77 | 0,77 |
| 200 | 128 | 0,64 |
| 398 | 207 | 0,52 |

O sea que "K = 200" no decía *"este organismo necesita 200 proteínas"*, decía *"bajaste 400 proteínas"*. La medición era circular: la referencia estaba construida con las mismas proteínas que se estaban evaluando.

Con el proteoma humano completo el test se aprueba:

| Corpus | K | K/corpus |
|---|---|---|
| 200 | 130 | 0,65 |
| 400 | 203 | 0,51 |
| 1.000 | 304 | 0,30 |
| 2.500 | 462 | 0,18 |
| 6.000 | 592 | 0,10 |
| 20.591 | 481 | 0,02 |

A partir de ~2.500 proteínas K deja de seguir al corpus y se estanca en 460–590. Ese plateau es el número que significa algo: **humano necesita del orden de 500 proteínas**, no las ~200 que reportaba la versión anterior.

Como control externo independiente: para una fuente uniforme, la distancia esperada por puro error de muestreo es 1,74/√R, que cruza el umbral de 0,01 en R ≈ 30.300. El N medido para el control al azar es ~31.400. Antes, con el corpus chico, la curva daba ~21.000 — sesgado hacia abajo, exactamente en la dirección que predice la contaminación. Ahora teoría y medición coinciden.

### Cómo reproducirlo

```
python3 main.py 4b
```

Código en `exercises/ex4/b.py`; las funciones de convergencia están en `stats.py` (`convergence_curve` para el eje 1, `count_convergence_curve` para el eje 2).

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
  Random DNA       0.000       0.155       0.165       0.108       0.131
     E. coli       0.155       0.000       0.122       0.090       0.063
       Yeast       0.165       0.122       0.000       0.101       0.072
       Human       0.108       0.090       0.101       0.000       0.045
     Natural       0.131       0.063       0.072       0.045       0.000
```

Lo importante no son los valores sino su escala relativa: **la distancia entre organismos (0,090–0,122) es del mismo orden que la distancia de cada organismo al azar (0,108–0,165)**. *Yeast* está a 0,122 de *E. coli* y a 0,165 del control: apenas 1,4× más lejos del azar que de otro organismo real.

Es decir, "natural" no es un punto sino una nube ancha, y cualquier detector que use un único perfil de referencia hereda esa limitación.

### Qué tan bien separa la métrica una secuencia individual

La pregunta que importa para 4e no es cuánto se diferencian dos corpus grandes, sino si la métrica puede decidir sobre **una** secuencia. Se midió con AUC: se toman fragmentos de proteínas reales y controles al azar del mismo largo, se calcula la distancia de cada uno al perfil Natural, y se mide la probabilidad de que el fragmento real quede más cerca. 0,5 es azar, 1,0 es separación perfecta.

![Discriminación según el largo](exercises/ex4/orf_discrimination.png)

```
  length      all 20   E R K N C         R C
      16       0.490       0.537       0.614
      30       0.445       0.566       0.624
      60       0.473       0.631       0.703
     120       0.490       0.712       0.800
     250       0.568       0.820       0.892
     400       0.609       0.879       0.936
```

AUC de cada aminoácido por separado, con fragmentos de 250 residuos:

```
  R: 0.905   E: 0.722   C: 0.691   K: 0.629   D: 0.623   S: 0.568   H: 0.534   M: 0.497
```

### Discusión

**Usar los 20 aminoácidos es casi inútil.** Se queda en 0,44–0,57 hasta los 250 residuos: prácticamente indistinguible de tirar una moneda. La razón es de ruido, no de biología: en un fragmento de 250 residuos cada frecuencia tiene un error de muestreo de ±0,013, y sumar veinte términos de error tapa las cinco o seis diferencias que sí son reales.

**Un subconjunto chico funciona mucho mejor.** Con sólo R y C la AUC sube a 0,892 en 250 residuos y 0,936 en 400. R y C son justamente los casos donde el código genético predice frecuencia alta y la biología entrega baja (R: 0,098 esperado contra 0,052 observado; C: 0,033 contra 0,016), así que la discrepancia es grande y consistente en los tres organismos.

**La composición no sirve en ORFs cortos.** A 16 residuos —la mediana de los ORFs medidos en 3c— hasta el mejor subconjunto da 0,61, y los 20 juntos dan 0,49, indistinguible del azar. La señal recién se vuelve usable arriba de ~120 residuos.

### Consecuencia para 4e

1. **El largo es la señal principal.** Es lo único que discrimina en ORFs cortos, y es donde 3c ya mostró la diferencia más fuerte (máximo real 1121 contra 231 del control).
2. **La composición es señal secundaria**, y sólo aporta arriba de ~100 residuos.
3. **Conviene usar pocos aminoácidos, no los 20.** Incluir los que no discriminan degrada activamente el resultado.

Estos tres puntos justifican la forma del score de 4e: probabilidad base por largo, corregida por composición sólo cuando el ORF es lo bastante largo.

### Advertencia metodológica

El subconjunto R/C se eligió mirando los mismos datos con los que después se lo evaluó, así que 0,936 es optimista. La versión rigurosa elige los aminoácidos usando dos organismos y evalúa sobre el tercero, dejado afuera.

### Cómo reproducirlo

```
python3 main.py 4d
```

Código en `exercises/ex4/d.py`.
