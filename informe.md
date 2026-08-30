# TP 1: análisis de secuencias

Cada ejercicio se corre por separado con `python3 main.py <ejercicio>`, por ejemplo `python3 main.py 4a`. El código está en `exercises/`, las funciones compartidas en `sequences.py`, `stats.py`, `ncbi.py` y `plots.py`. Las secuencias quedan cacheadas en `data/`; si se borran, la primera corrida vuelve a bajar los tres proteomas.

---

## Datos

Las proteínas reales se bajan de **Swiss-Prot**, la mitad revisada a mano de UniProt, donde hay una entrada por gen por organismo y las isoformas figuran como anotaciones adentro de la entrada. Se accede con el mismo cliente Entrez que usa el resto del TP, cambiando el filtro de la query a `swissprot[filter]`.

La primera versión usaba RefSeq y la cambiamos. RefSeq es genómico: un pipeline de anotación emite un registro por cada gen predicho de cada genoma secuenciado, así que la misma proteína aparece repetida por isoforma y por cepa, y conviven con entradas sin verificar. La query `"Escherichia coli"[Organism] AND refseq[filter]` devuelve 6.666.151 proteínas contra 23.300 en Swiss-Prot, y entre los primeros resultados humanos había cuatro isoformas de una misma proteína seguidas de cinco de otra. Para *E. coli* usamos además la cepa K-12 y no la especie entera, porque `"Escherichia coli"[Organism]` matchea todas las cepas secuenciadas y eso reintroduce la redundancia por cepa aun dentro de Swiss-Prot.

De cada organismo se baja el proteoma revisado completo, no una muestra:

| Organismo | Proteínas | Residuos | Largo medio | Largo medio conocido |
|---|---|---|---|---|
| *E. coli* K-12 | 6.062 | 1.848.964 | 305 | ~316 |
| *S. cerevisiae* | 7.888 | 3.544.698 | 449 | ~450 |
| *H. sapiens* | 20.591 | 11.454.685 | 556 | ~560 |

Los tres largos medios coinciden con el valor conocido de cada proteoma, que es la comprobación de que la descarga no quedó sesgada. Bajar todo evita tener que elegir una muestra, y eso importa porque `esearch` devuelve los ids ordenados por fecha de depósito: tomar los primeros *n* hubiera significado quedarse con lo último depositado. Medido sobre K-12, los primeros 400 ids dan un largo medio de 256 residuos contra 304 de 400 sorteados al azar. El pool de ids se pagina, porque `esearch` devuelve como máximo 10.000 por request y el proteoma humano tiene 20.616 entradas.

Los totales quedan algo por debajo del número de entradas de cada query porque se descartan las secuencias con caracteres no estándar, como selenocisteína o posiciones ambiguas.

---

## 4a

> **4a)** Compare Distribución de aminoácidos de secuencia de proteínas al azar vs. secuencias de proteínas reales.

Se comparan cuatro fuentes con la misma cantidad de residuos (1.848.964, el tamaño del proteoma de *E. coli* K-12, que es el más chico de los tres) para que ninguna tenga ventaja por volumen de datos: una fuente al azar y los proteomas revisados de *E. coli*, levadura y humano. A los tres organismos combinados en partes iguales los llamamos *Natural*; esa es la distribución de referencia que después usa 4e.

La fuente al azar no se arma sorteando aminoácidos con igual probabilidad sino traduciendo ADN aleatorio y descartando los codones stop. La diferencia importa porque el código genético no reparte los codones de manera pareja (leucina tiene seis y metionina uno solo): sortear aminoácidos uniformemente daría una fuente que no se parece a ninguna secuencia obtenida de ADN, y la comparación terminaría midiendo esa desigualdad del código antes que lo que distingue a una proteína real. El caso uniforme igual aparece en las figuras, como la línea punteada en 0,05. Levadura y humano tienen más residuos que *E. coli*, así que se los recorta, y las listas se mezclan antes de hacerlo para que el recorte no termine quedándose con las proteínas depositadas más recientemente.

![Distribución de aminoácidos](exercises/ex4/aa_distribution.png)

![Perfil natural vs secuencias al azar](exercises/ex4/aa_natural_vs_random.png)

Medidas con la métrica de 4c, las distancias al perfil Natural son 0,045 para humano, 0,063 para *E. coli*, 0,072 para levadura y 0,131 para la fuente al azar.

Leucina aparece con la misma frecuencia en las dos fuentes, 0,098 al azar y 0,100 en los organismos, así que la cantidad de leucina de una secuencia no dice nada sobre su origen. Arginina es el caso opuesto: 0,098 al azar contra 0,052 en los organismos, una diferencia lo bastante grande como para distinguirlas. Lo mismo ocurre con cisteína (0,033 contra 0,016) y, en sentido inverso, con glutámico, aspártico y lisina, que los organismos usan cerca del doble. Los aminoácidos donde las dos fuentes coinciden no aportan información; los que se separan son los candidatos para el detector de 4e.

El resultado que más condiciona lo que sigue es que los organismos se diferencian entre sí casi tanto como del azar: *E. coli* y levadura están a 0,122, contra los 0,131 que separan a Natural del azar. Dicho de otro modo, dos proteínas reales de organismos distintos pueden estar tan lejos entre sí como una real de una al azar, así que usar un único perfil de referencia deja un margen de apenas 1,8×. Para 4e esto significa que la composición de aminoácidos por sí sola no va a alcanzar, y que habrá que explorar otros criterios.


---

## 4b

> **4b)** Analice cómo cambian las distribuciones al aumentar el tamaño de la secuencia "al azar" analizada, y al incrementar el número de secuencias reales analizadas. ¿Cuándo es suficiente?

La consigna nombra dos ejes que no son el mismo: se puede crecer en cantidad de residuos o en cantidad de proteínas enteras. Los medimos por separado, pero contra el mismo blanco, que es lo que hace comparables los dos resultados. Ese blanco es el proteoma completo del organismo, que es lo que significa "la distribución de este organismo"; para la fuente al azar es la distribución exacta de su modelo nulo, calculada a partir de la tabla de codones. En ambos casos se toma una muestra parcial, se mide su distancia al blanco con la métrica de 4c, y se busca dónde cae por debajo de 0,01. Llamamos **N** al número de residuos donde se cruza el umbral y **K** al número de proteínas enteras. El eje de proteínas aplica sólo a los organismos, porque la fuente al azar se genera como una única cadena y no tiene "secuencias". Cada curva se promedia sobre 5 barajadas del orden, para que no dependa de qué proteínas vinieron primero.

![Estabilización de la distribución](exercises/ex4/aa_stabilized.png)

Las cuatro fuentes se estabilizan en el mismo orden de magnitud, alrededor de 3·10⁴ residuos: 34.200 la fuente al azar, 34.000 *E. coli*, 28.400 levadura y 32.400 humano. En proteínas enteras alcanza con unas 252 para *E. coli*, 220 para levadura y 406 para humano.

Los dos ejes no dan el mismo número. Las 406 proteínas humanas suman unos 226.000 residuos, siete veces el N del eje 1, así que hacen falta más residuos para estabilizar la distribución cuando llegan agrupados en proteínas enteras que cuando se cuentan sueltos. La brecha además es mayor cuanto más largas son las proteínas del organismo: 2,3 en *E. coli*, que promedia 305 residuos por proteína, 3,5 en levadura con 449 y 7,0 en humano con 556.

![Distribuciones a distintos tamaños de muestra](exercises/ex4/aa_sample_sizes.png)

Los tres gráficos de la figura de arriba muestran lo mismo de otra forma, con la distribución completa a tres tamaños de muestra. Con R = 1000, muy por debajo de N, las barras son ruido: levadura da K en 0,096 cuando su valor real es 0,073, y *E. coli* da L en 0,119 contra 0,106. Con R ≈ N y con R >> N los dos gráficos de abajo son casi idénticos, con diferencias en el tercer decimal: el del medio ya contiene toda la información que aporta el último, usando el 1,7% de los datos.

Conviene aclarar que N y K son ruidosos. Cada corrida mezcla los datos en un orden distinto y da un resultado un poco diferente, y como la curva llega al umbral casi horizontal, una variación mínima en su altura mueve bastante el punto donde lo cruza: sobre 6 corridas de la fuente al azar, N va de 21.900 a 34.400. Lo que se sostiene es el orden de magnitud y la relación entre los dos ejes, no la cifra exacta de una corrida.

---

## 4c

> **4c)** Elija una métrica que permita comparar dos distribuciones. (p.ej., RMSE).

Usamos la distancia de variación total, implementada en `stats.distributions_distance()`:

```
TV(p, q) = 0,5 · Σ |p(a) − q(a)|
```

La elegimos sobre RMSE porque está acotada entre 0 y 1, así que el valor se interpreta sin referencia externa, y porque se lee directamente como qué fracción de la masa de probabilidad hay que mover para convertir una distribución en la otra. Tampoco depende de la cantidad de bins, y permite restringirla a un subconjunto de aminoácidos sumando menos términos, que es lo que se usa en 4d.

---

## 4d

> **4d)** Escriba un código que dadas dos distribuciones las compare y obtenga la métrica correspondiente. Utilícela para: i) ver cómo evolucionan las distribuciones obtenidas en 4a al aumentar el tamaño de la muestra (evaluar cuándo la distribución deja de cambiar). ii) estudiar de manera sistemática las diferencias entre las distribuciones obtenidas en 4a para diferentes muestras (comparar distribuciones de secuencias naturales y al azar).

La parte (i) es el análisis de convergencia de 4b, donde la métrica se aplica a cada fuente contra su propia distribución final. Esta sección es la parte (ii).

```
              Random DNA     E. coli       Yeast       Human     Natural
  Random DNA       0.000       0.155       0.165       0.108       0.131
     E. coli       0.155       0.000       0.122       0.090       0.063
       Yeast       0.165       0.122       0.000       0.101       0.072
       Human       0.108       0.090       0.101       0.000       0.045
     Natural       0.131       0.063       0.072       0.045       0.000
```

La matriz confirma con todos los pares lo que en 4a se veía en uno: las distancias entre organismos (0,090 a 0,122) están en el mismo orden que las distancias al azar (0,108 a 0,165). Levadura está a 0,122 de *E. coli* y a 0,165 del azar, apenas 1,4 veces más lejos.

Para 4e, sin embargo, la pregunta no es cuánto se diferencian dos conjuntos grandes sino si la métrica alcanza para decidir sobre una secuencia sola. Para medirlo tomamos fragmentos de proteínas reales y secuencias al azar del mismo largo, calculamos la distancia de cada uno al perfil Natural, y contamos con qué frecuencia el fragmento real queda más cerca que el aleatorio. Llamamos a eso la tasa de aciertos: 0,5 es lo que daría elegir al azar y 1,0 es separación perfecta.

![Discriminación según el largo](exercises/ex4/orf_discrimination.png)

```
  largo       20 aa   E R K N C         R C
      16       0.489       0.535       0.618
      30       0.474       0.570       0.636
      60       0.454       0.638       0.711
     120       0.467       0.694       0.777
     250       0.574       0.828       0.895
     400       0.607       0.877       0.935
```

Usar los veinte aminoácidos es prácticamente inútil: se queda entre 0,45 y 0,57 hasta los 250 residuos. El motivo es ruido, no composición: en un fragmento de 250 residuos cada frecuencia tiene un error de muestreo, y sumar veinte términos de error tapa las cinco o seis diferencias que sí existen. Restringiendo la métrica a los aminoácidos que en 4a se separaban entre las dos fuentes, la tasa sube a 0,895 con R y C solos, y midiendo aminoácido por aminoácido el orden es R 0,900, E 0,733, C 0,699, K 0,636 y D 0,605.

El límite aparece con los fragmentos cortos. A 16 residuos, que es la mediana de los ORFs medidos en 3c, el mejor subconjunto da 0,62 y los veinte juntos dan 0,49, indistinguible del azar. La señal recién se vuelve usable arriba de unos 120 residuos.

De ahí salen tres cosas para 4e. El largo es la señal principal, porque es la única que discrimina en ORFs cortos y porque es donde 3c mostró la diferencia más marcada (máximo real de 1121 residuos contra 231 en las secuencias al azar). La composición es una señal secundaria que sólo aporta valor arriba de unos 100 residuos. Y conviene usar pocos aminoácidos en lugar de los veinte, porque incluir los que no discriminan empeora el resultado.

Queda una objeción posible: R y C se eligieron porque eran los que mejor separaban en estos datos, y después se midió con ellos sobre esos mismos datos, así que parte del 0,935 podría ser casualidad de la muestra contada como acierto. Para descartarlo repetimos el experimento dejando un organismo afuera. Se busca el mejor par entre los 190 posibles usando sólo los otros dos, se arma el perfil de referencia también sólo con esos dos, y recién entonces se mide sobre el organismo que quedó afuera, que no participó ni de la elección ni de la referencia.

| Organismo afuera | Par elegido | Tasa en los dos de entrenamiento | Tasa en el que quedó afuera |
|---|---|---|---|
| *E. coli* | E R | 0,904 | 0,936 |
| Levadura | E R | 0,850 | 0,882 |
| Humano | C R | 0,959 | 0,886 |

Arginina aparece en los tres pares y el segundo aminoácido alterna entre glutámico y cisteína, así que la elección no depende del organismo que se mire. Las tasas sobre el organismo dejado afuera van de 0,882 a 0,936, el mismo rango que las de entrenamiento, de modo que el 0,935 de la tabla anterior no era producto de haber elegido mirando los mismos datos. Que en dos de los tres casos la tasa sea mayor sobre el organismo dejado afuera no es raro: los organismos difieren en cuán fácil es distinguirlos del azar, y eso pesa más que la diferencia entre elegir y evaluar.

