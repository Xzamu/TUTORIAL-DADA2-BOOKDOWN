---
output:
  bookdown::gitbook:
    css: style.css
---

<style>
.rmdnote {
  padding: 10px;
  margin-bottom: 15px;
  border: 1px solid #abcdef;
  background-color: #f3f3f3;
  border-radius: 5px;
}
</style>

# Introducción a DADA2

Repositorio de GitHub [aquí](https://github.com/benjjneb/dada2)

Cita:

> Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., y Holmes, S. P. (2016). 
> DADA2: High-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583.

Descripción:
DADA2 es un paquete de software de código abierto utilizado para modelar y corregir errores en secuencias de amplicones  secuenciados con el protocolo Illumina. Infiere secuencias de una muestra de manera exacta y resuelve diferencias con una resolución de un nucleótido.


Los datos utilizados pueden ser consultados en el [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) bajo el nombre de proyecto [PRJNA428495](https://www.ebi.ac.uk/ena/browser/view/prjna428495).

Descargamos los archivos FASTQ en la sección "Generated FASTQ files: FTP" y los guardamos en una carpeta.

Los binarios del paquete DADA2 están disponibles a través de 
[Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)

## Instalación de DADA2 en R (versión 4.3)
Este bloque instala el paquete DADA2 a través de Bioconductor. Es importante tener instalado BiocManager para poder acceder a los paquetes de Bioconductor.


```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dada2")
library(dada2)
```

## Configuramos el directorio en el que se ubican las reads

```r
path <- "/Users/" # Cambiamos al directorio que contiene las reads
list.files(path)
```

## Ordenamos los archivos y extraemos el nombre de las muestras
Este bloque ordena los archivos FASTQ y extrae el nombre de las muestras. Asegúrate de que los patrones de búsqueda (pattern="_1.fastq.gz" y pattern="_2.fastq.gz") coincidan con la nomenclatura de tus archivos.


```r
fnFs <- sort(list.files(path, pattern="_1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz"))
sample.names <- sapply(strsplit(fnFs, "_"), '[', 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
```

## Establecer ruta para archivos filtrados
Define una ruta basada en la variable ``sample.names`` donde guardar los archivos filtrados. El código genera nombres de archivo con un sufijo ("_F_filt.fastq.gz" o "_R_filt.fastq.gz") para diferenciarlos.


```r
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
print(fnFs)
print(fnRs)
print(filtFs)
print(filtRs)
```

## Filtrar y recortar lecturas
La función filterAndTrim se utiliza para el filtrado y recorte de las lecturas, es decir, se eliminan las lecturas de baja calidad y las recorta a una longitud específica.
Se especifican varios parámetros, como las longitudes de truncamiento, la calidad máxima permitida, la eliminación de secuencias de fagos (``rm.phix``), la compresión y el uso de múltiples hilos (``multithread``). El resultado se almacena en la variable ``out``.
Windows 10 permite el 'multi-theading', excepto en el comando ``filterAndTrim``.


```r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```

## Estimación de la tasa de error
Se estima la tasa de error para las lecturas adelante (errF) y las lecturas reversas (errR) utilizando la función ``learnErrors``. 
La opción ``multithread=TRUE`` se utiliza para acelerar el proceso.


```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

## Dereplicar lecturas
Se realiza la eliminación de duplicados de lecturas (dereplicación) para las lecturas adelante (derepFs) y las lecturas reversas (derepRs) utilizando la función ``derepFastq``. Se asignan nombres a las muestras con ``sample.names``.


```r
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

## Inferencia de secuencias
Se utilizan las tasas de error estimadas en el paso anterior para la inferencia de secuencias únicas utilizando la función ``dada`` para las lecturas adelante (dadaFs) y las lecturas reversas (dadaRs).


```r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

## Fusionar lecturas emparejadas
Aquí se fusiona las secuencias emparejadas utilizando la función ``mergePairs``. 
Las secuencias únicas y las lecturas originales se utilizan para la fusión, y los resultados se almacenan en la variable mergers.


```r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```

## Crear tabla de secuencias
Se crea una tabla de secuencias utilizando la función ``makeSequenceTable``. 
Esta tabla contendrá información sobre las secuencias y su abundancia.


```r
seqtab <- makeSequenceTable(mergers)
```

## Eliminar quimeras
Las quimeras son secuencias que resultan de la combinación de dos o más secuencias parentales diferentes durante el proceso de amplificación por PCR. Estas formaciones ocurren en ciclos posteriores de PCR cuando hay una alta concentración de cebadores parcialmente extendidos que compiten con los cebadores originales.

Se realiza la eliminación de quimeras utilizando la función ``removeBimeraDenovo``. 
Se especifica el método de eliminación como "consensus".


```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```


Las quimeras puede introducir artefactos¿¿'¿'¿ en los análisis de datos, dando lugar a interpretaciones erróneas.

En el análisis de datos usando DADA2, las quimeras se eliminan utilizando el código removeBimeraDenovo. Esta función es una interfaz de conveniencia para la eliminación de quimeras. DADA2 ofrece varios métodos para identificar quimeras, como la identificación a partir de secuencias agrupadas y la identificación por consenso entre muestras. Al utilizar el método de consenso, por ejemplo, las muestras en una tabla de secuencias se verifican independientemente en busca de quimeras, y se toma una decisión de consenso sobre cada variante de secuencia.

DADA2 es muy meticuloso al tratar con estas quimeras para asegurar que los datos analizados sean precisos. Por ejemplo, al utilizar removeBimeraDenovo con el método "pooled", todas las muestras en la tabla de secuencias se agrupan para la identificación de quimeras. Si se usa el método "consensus", las muestras en una tabla de secuencias se verifican independientemente en busca de quimeras, y se toma una decisión de consenso sobre cada variante de secuencia. Esto es vital ya que las quimeras pueden variar entre las muestras y es esencial asegurarse de que no afecten los resultados.

En conclusión, eliminar quimeras o bimeras es esencial para obtener una representación precisa de las comunidades microbianas en los análisis de datos. Si no se eliminan, podrían llevar a interpretaciones erróneas de la estructura y diversidad de la comunidad. DADA2 proporciona herramientas efectivas para identificar y eliminar estas quimeras, asegurando así la calidad y precisión de los resultados obtenidos.

## Asignar taxonomía
Se asigna la taxonomía a las secuencias utilizando la función ``assignTaxonomy``. 
Se proporciona un archivo de referencia para la asignación taxonómica, y se utiliza el ``multithread`` para acelerar el proceso. La base de datos taxonómica utilizada será un archivo fasta de entrenamiento derivado de la versión 138.1 del [Proyecto Silva](https://www.arb-silva.de/) y formateado para su uso con DADA2. Este archivo puede ser descargado de [Zenodo](https://zenodo.org/records/4587955).


```r
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/UBBBC/Desktop/Samuel/PRJNA428495/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
```

## Impresión de resultados
Se imprime en la consola las primeras filas de la asignación taxonómica resultante.


```r
unname(head(taxa))
```

## Alinear secuencias
'getSequences' extrae las secuencias de la tabla sin quimeras.
'names(sequences)' asigna a cada secuencia su propio nombre para su fácil identificación.


```r
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences
```

'AlignSeqs' de la biblioteca 'DECIPHER' alinea las secuencias de ADN.
'DNAStringSet' convierte las secuencias en un formato adecuado para el alineamiento.
'anchor=NA' indica que no se usará un anclaje durante el alineamiento.


```r
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
```
## Crear matriz de distancia
'phyDat' convierte las alineaciones en un formato que puede ser utilizado por funciones de filogenética.
'type="DNA"' especifica que los datos son secuencias de ADN.
'dist.ml' calcula una matriz de distancias utilizando el método de máxima verosimilitud.


```r
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
```

## Realizar Neighbor Joining
'NJ' realiza un análisis filogenético utilizando el método Neighbor Joining.
La salida es un árbol filogenético basado en la matriz de distancias.


```r
treeNJ <- NJ(dm)  # NJ es Neighbor Joining
```




```r
fit <- pml(treeNJ, data=phang.align)
```

'pml' (de la biblioteca phangorn) ajusta un modelo de máxima verosimilitud al árbol filogenético.
'treeNJ' es el árbol generado por Neighbor Joining.
'data=phang.align' son los datos de alineación convertidos previamente con 'phyDat'.

## Actualizamos el modelo con nuevos parámetros

```r
fitGTR <- update(fit, k=4, inv=0.2)
```

``update``: Actualiza un modelo previamente ajustado (fit) con nuevos parámetros.
``k=4``: Define el número de categorías de tasa para el modelo gamma de tasas de sustitución de nucleótidos variables entre sitios.
``inv=0.2``: Establece la proporción de sitios invariables en 0.2.

"JUSTIFICACIÓN"

## Optimizamos el modelo GTR con parámetros de inversión y tasa de gamma


```r
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
```
optim.pml^[
Optimiza un modelo filogenético basado en máxima verosimilitud.
]


## Desacoplamos el paquete phangorn para evitar conflictos


```r
detach("package:phangorn", unload=TRUE)
```

## Leemos los metadatos de las muestras desde un archivo CSV


```r
samdf <- read.csv("metadata.csv", header=TRUE, row.names = 1)
```

## Verificamos que los nombres de las filas de la tabla de secuencias coincidan con los metadatos


```r
rownames(seqtabNoC) %in% rownames(samdf)
all(rownames(seqtabAll) %in% samdf$run)
```

## Construimos un objeto phyloseq con la tabla de OTU, los metadatos de las muestras, la tabla taxonómica y el árbol filogenético


```r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
```


Phyloseq `<span class="tippy html-widget html-fill-item-overflow-hidden html-fill-item" height="500" id="htmlwidget-e3fe14b1b2d1d892ca76" width="960"></span>
<script type="application/json" data-for="htmlwidget-e3fe14b1b2d1d892ca76">{"x":{"opts":{"content":"Crea un objeto phyloseq que integra la tabla de OTU, los datos de las muestras, la información taxonómica y el árbol filogenético."},"text":"¿Qué hace?"},"evals":[],"jsHooks":[]}</script>`{=html}.



## Eliminamos las muestras sintéticas del objeto phyloseq


```r
ps <- prune_samples(sample_names(ps) != "Mock", ps)
```

## Imprimimos el objeto phyloseq resultante
ps
