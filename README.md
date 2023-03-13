# **Proyecto Final WBDS** || Enero-Marzo 2023
#### > Aschettino Giovanna

## Objetivo del proyecto
El análisis de expresión diferencial de genes a partir de datos de RNA-seq se realiza casi exclusivamente en el lenguaje R. Con el objetivo de aplicar lo aprendido en el curso y tomando de ejemplo información contenida en **Gene Expression Omnibus** de NCBI se genera un script para el análisis de expresión diferencial de genes en lenguaje python. 

## Instalación
* Python versión >= 3.8
* Paquetes de python:
  * adjustText
  * datetime
  * bioinfokit
  * matplotlib
  * numpy
  * os
  * pandas
  * pickle
  * plotly
  * pydeseq2

## Solución de problemas

- Instalación de paquetes:

    Se pueden visualizar los paquetes instalados con el siguiente comando:

         pip list

    De ser necesario se lo instala con 

         pip install <nombre del paquete>

- Compatibilidad del dataframe inicial.

    A modo de ejemplo se usó una matriz con el conteo de reads por gen descargada de GSO (GSE152075). Este archivo inicial tiene un formato de matriz PacientesxGenes. De tener otro tipo de input se debería reajustar la entrada al formato.
    La matriz de clasificación tiene dos columnas, Paciente y clasificación. En el ejemplo sale del mismo archivo inicial, de tener otro dataframe debería acomodarse en df_clase
    Cuando se quiera usar otro dataframe se deberá indicar en la tercer sección el path al archivo nuevo
    
```
df_conteo=pd.read_csv("InputData/GSE152075_raw_counts_GEO.txt", sep=' ', header=0)
```
    Si se tiene un archivo de clasificación se deberá agregar de la misma manera que el anterior, reemplazando el bloque:
```
df_clase=pd.DataFrame(columns=['id', 'condicion'])
df_clase.id=df_conteo.columns
df_clase.condicion=df_clase.id.str.split('_').str[0].replace("POS", "Positivo").replace("NEG", "Negativo")
df_clase=df_clase.set_index('id').rename_axis(None)
```
por
```
df_clase=pd.read_csv("**Nuevo path**", sep=' ', header=0)
```

- Path de trabajo:
El programa va a correr en el home a menos que se indique donde se quiere correr como una lista en la variable. si está trabajando en más de una pc y se tienen diferentes paths se pueden agregar a la lista separandolos por coma, siempre respetando que al final se encuentre el directorio actual. 

```pathlist=["**path1**","**path2**", os.getcwd()]]```


## Input y Output data

#### Input
- Matriz de conteo de reads por gen
- Tabla con la clasificación de cada paciente

#### Output
- Histograma con el comportamiento del p-valor
- Tabla con los 30 genes upregulados ordenados por el p-valor ajustado y por el log fold change
- Tabla con los 30 genes downregulados ordenados por el p-valor ajustado y por el log fold change
- Volcano plot sin anotar
- Volcano plot con los genes más significativos anotados.

## Pasos del script

El acceso a la información puede venir de múltiples plataformas. 

En esta oportunidad se decidió buscar información en ncbi y a partir de ahí buscar la información en GEO _https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi_. 

Con el id del proyecto se puede realizar la búsqueda y nos lleva a una página similar como la siguiente: 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075

Aquí se puede descargar al final de la página el archivo *_raw_counts_GEO.txt.gz que va a contener la información necesaria.

Una vez realizada la descarga y la descompresión del archivo, se puede proceser a correr el script ```wbds_DEA.py```, teniendo en cuenta de hacer las modificaciones de cambio de nombre de archivo comentado en el paso anterior. Este script se puede dividir conceptualmente en 10 secciones:
1. Carga de paquetes
2. Seteo de path de trabajo
3. Lectura de archivos:
  El archivo inicial contiene el conteo de lecturas por gen, por paciente.
  Se genera un archivo con la clasificación de los pacientes en base al nombre de muestra. Es específico para ese dataframe. Es necesario contar con esa clasificación para poder continuar.
4. Exploración y filtrado de datos
  Se aplican criterios para el filtrado de pacientes y de genes de acuerdo a las características del protocolo. 
  En este caso se elije eliminar los pacientes sin clasificación (criterio usado con frecuencia) y eliminar aquellos genes sin conteo en ningún paciente, dependiendo el protocolo, este valor podría ser mayor.
5. Expresión diferencial de genes
6. Análisis estadístico de la expresión diferencial
7. Evaluación inicial de resultados.
  Se observa el pvalor y el pvalor ajustado. Se visualiza mediante histogramas.
8. Reajuste del log fold change para poder ser visualizado graficamente.
9. Generación de tablas con los genes más influyentes tanto para aquellos downregulados como los up regulados.
10. Visualización de información a partir de Volcano Plots. 

## Resultados

Los datos obtenidos de la base de datos GEO corresponden a un estudio realizado para evaluar la respuesta antiviral en pacientes con SARS-CoV-2.
En base al análisis obtenido a partir del script y la publicación incial se puede ver que los genes significativos estadísticamente se correlacionan con los publicados, siendo importantes los genes asociados a respuesta antiviral.

### Listado con los 30 genes down reculados más significativos

![Tabla_Genes_down_bypaj](https://user-images.githubusercontent.com/54379644/224775418-af4eb28a-9531-47db-8fc0-3094b5623c59.png)

### Listado con los 30 genes up reculados más significativos

![Tabla_Genes_Up_bypadj](https://user-images.githubusercontent.com/54379644/224775442-bad27b0c-7a8f-4221-a2b7-4a03fb39f884.png)


## Links a bibliografía consultada

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075
- https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000849
- https://www.sciencedirect.com/science/article/pii/S0168170223000151#bib0041
- https://www.sciencedirect.com/science/article/pii/S258900422030777X

##
Proyecto Final WBDS - Análisis de datos de expresión de genes con python by Giovanna Aschettino is marked with CC0 1.0. To view a copy of this license, visit http://creativecommons.org/publicdomain/zero/1.0

