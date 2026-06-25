# MADOBIS_Project
Bioinformatic Identification of Novel Open Reading Frames in Differentially Expressed Unassigned Transcripts in Breast Cancer

Trabajo de Fin de Máster — Máster en Análisis de Datos Ómicos y Biología de Sistemas (MADOBIS)
Autor: Marcos Antonio Dorantes Benitez



## Descripción

Este repositorio contiene el código desarrollado para el Trabajo de Fin de Máster centrado en la identificación bioinformática de nuevos ORFs (Open Reading Frames) en transcritos no anotados que se expresan de forma diferencial entre tejido sano y tejido tumoral en cáncer de mama (TCGA-BRCA).

El proyecto combina dos grandes bloques de trabajo:


Análisis transcriptómico y multi-ómico (TFM_DorantesBenitez.Rmd): procesamiento de datos de expresión a nivel de isoforma, análisis exploratorio, expresión diferencial, anotación funcional, identificación de smORFs y modelado de redes de regulación.
Predicción estructural de péptidos candidatos (AlphaFoldGenerico_DorantesBenitez.ipynb): predicción de la estructura 3D de los péptidos codificados por los smORFs identificados como relevantes, mediante AlphaFold.



## Estructura del repositorio

MADOBIS_Project/

├── TFM_DorantesBenitez.Rmd                  # Pipeline completo en R (análisis ómico)

├── AlphaFoldGenerico_DorantesBenitez.ipynb  # Notebook de predicción estructural (AlphaFold)

└── README.md                                # Descripción del repositorio


Nota: los datos de entrada (matrices de expresión TCGA-BRCA, metadata de muestras, etc.) no se incluyen en el repositorio por su tamaño. Ver sección "Data and Code Availability Statement" en el manuscrito.



## TFM_DorantesBenitez.Rmd

Script de R Markdown que implementa el pipeline completo de análisis, dividido en las siguientes secciones:

### 1.1. Tratamiento de datos
Carga de librerías de Bioconductor/CRAN (NOISeq, biomaRt, clusterProfiler, FactoMineR, mixOmics, MORE, igraph, etc.).
Carga y limpieza de la matriz de expresión a nivel de isoforma de TCGA-BRCA (RNA-seq, RSEM normalizado) y de la metadata de muestras.
Identificación de viales por paciente a partir de los códigos de barras TCGA y emparejamiento de muestras (tumor primario vs. tejido normal adyacente del mismo paciente).

### 1.2. Análisis exploratorio (PCA)
PCA de réplicas técnicas (viales 01A/01B y 11A/11B) para descartar efecto lote.
PCA global (Normal / Tumor primario / Metástasis).
PCA de muestras pareadas Tumor vs. Normal — base para el análisis de expresión diferencial.

### 1.3. Expresión diferencial
*Test de Wilcoxon pareado* (una muestra, H₀: mediana de la diferencia = 0) sobre la matriz diferencial (Sano − Tumoral) por paciente, con filtrado de isoformas con exceso de ceros, cálculo del estadístico de contraste (Z) y corrección de p-valores por FDR.
*NOISeq como método complementario* (no paramétrico, basado en log-fold-change y probabilidad de diferencial expresión), incluyendo control de calidad exploratorio (RNA composition, distribución de counts, sensitivity plot).
Comparación e intersección de isoformas significativas entre ambos métodos (diagramas de Venn).

### 1.4. Anotación funcional y caracterización de transcritos
Consulta a Ensembl (biomaRt) y UCSC Genome Browser (RMySQL) para recuperar coordenadas, biotipo y anotación de los transcritos diferencialmente expresados.
Identificación de transcritos no anotados / lncRNA con potencial codificante no descrito.
Consulta a la API REST de sORFs.org para la detección de small ORFs (smORFs) dentro de dichos transcritos.
Análisis de enriquecimiento funcional (GO / KEGG) con clusterProfiler y enrichplot.

### 1.5. Modelado multi-ómico / redes de regulación
CCA regularizado (rCCA) con mixOmics para evaluar la correlación entre la expresión de smORFs y genes codificantes.
MORE (PLS1 con selección de variables por permutación) para modelar relaciones regulador–diana y construir la red smORF–gen (networkMORE).
Exportación de redes de correlación a formato .gml con igraph para su visualización en Cytoscape.
Intersección de nodos entre las redes MORE y rCCA mediante diagramas de Venn.

### Salida
El resultado de este script es una lista de smORFs candidatos (transcritos no anotados, diferencialmente expresados entre tejido tumoral y sano, con evidencia de codificación y relación regulatoria con genes conocidos), que constituyen el input del siguiente bloque.


## AlphaFoldGenerico_DorantesBenitez.ipynb

Notebook genérico para la predicción de la estructura tridimensional de los péptidos codificados por los smORFs candidatos seleccionados en el análisis anterior, mediante AlphaFold.

Toma como entrada las secuencias peptídicas derivadas de los smORFs de interés.
Genera las predicciones de estructura 3D y las métricas de confianza asociadas (p. ej. pLDDT).
Permite una primera caracterización estructural de péptidos sin homología conocida, como evidencia adicional de su potencial función biológica.

Pensado para ejecutarse en Google Colab (recomendado, por la necesidad de GPU) o en un entorno local con GPU compatible.



## Requisitos

### Para TFM_DorantesBenitez.Rmd

R (≥ 4.2 recomendado) y RStudio.

Paquetes de CRAN: dplyr, ggplot2, ggplotify, ggrepel, ggvenn, VennDiagram, FactoMineR, factoextra, patchwork, reshape2, stringr, httr, jsonlite, igraph, RMySQL.

Paquetes de Bioconductor: NOISeq, biomaRt, clusterProfiler, enrichplot, GO.db, org.Hs.eg.db, AnnotationDbi, mixOmics, MORE, AnnotationDbi.

Instalación rápida:

```
install.packages(c("dplyr", "ggplot2", "ggplotify", "ggrepel", "ggvenn", 
                   "VennDiagram", "FactoMineR", "factoextra", "patchwork", 
                   "reshape2", "stringr", "httr", "jsonlite", "igraph", "RMySQL"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("NOISeq", "biomaRt", "clusterProfiler", "enrichplot", 
                       "GO.db", "org.Hs.eg.db", "AnnotationDbi", "mixOmics", "MORE"))
```                   

### Para AlphaFoldGenerico_DorantesBenitez.ipynb

Python ≥ 3.8

Jupyter Notebook / Google Colab

AlphaFold (o ColabFold) y sus dependencias asociadas (biopython, numpy, etc. según la implementación utilizada en el notebook)



## Datos de entrada

El script de R asume la presencia de los siguientes archivos en una carpeta local (no incluidos en el repositorio por tamaño/licencia):

- `BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt` matriz de expresión de isoformas TCGA-BRCA (nivel 3, normalizada por RSEM).
- `sample_type.csv` metadata con el tipo de muestra (tumor primario / normal / metástasis) por código de barras TCGA.
- `Galaxy_multiBigwigSummary_bin_counts.tabular` archivo tabular generado a partir de multiBigwigSummary en GALAXY (Version 3.5.4+galaxy0) usando los archivos generados por el propio pipeline `smORFs_hg19.bed` y los archivos BigWig tracks for hg19 (`PhyloCSF+0.bw, +1.bw, +2.bw`) del Broad Institute. 
- `filtered_sorfs.csv` listado de ORFs de la base de datos https://sorfs.org/database, siguiendo las indicaciones del manucristo.



## Uso

Descargar los datos de entrada y colocarlos en el directorio de trabajo o en una carpeta data/, siguiendo las indicaciones de las rutas del propio pipelina.
Abrir TFM_DorantesBenitez.Rmd en RStudio y ejecutar los chunks secuencialmente (las secciones más largas permiten guardar/cargar objetos intermedios en .RData para evitar recomputar pasos costosos).
Exportar la lista de smORFs candidatos y sus secuencias peptídicas.
Cargar dichas secuencias en AlphaFoldGenerico_DorantesBenitez.ipynb (idealmente en Google Colab) para obtener las predicciones estructurales.



## Contexto académico

Trabajo de Fin de Máster del Máster en Análisi de Datos Ómicos y Biología de Sistemas (MADOBIS).



## Licencia

Este repositorio se proporciona con fines académicos. Si deseas reutilizar el código, por favor cita el Trabajo correspondiente.
