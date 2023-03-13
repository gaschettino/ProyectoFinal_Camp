# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 17:54:05 2023

@author: Giovanna Aschettino

toy data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075
https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000849
https://www.sciencedirect.com/science/article/pii/S0168170223000151#bib0041
https://www.sciencedirect.com/science/article/pii/S258900422030777X

"""
# %% Paquetes a cargar
import os
import datetime
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pickle as pkl
import matplotlib.pylab as plt
from bioinfokit import visuz
import numpy as np
from adjustText import adjust_text

# from bioinfokit import analys
# import glob
# import seaborn as sns

# %% Inicializo la carpeta donde voy a trabajar.

print('Paso 1: Inicialización')
pathlist=["/home/usuario/OneDrive/RnaSeq", "C:/Users/Asus/Documents/ProyectoFinal_Camp", "/mnt/c/Users/Asus/Documents/ProyectoFinal_Camp"]
flag=False
for p in pathlist:
    if os.path.isdir(p) and not flag:
        flag=True
        os.chdir(p)
        path=p
        datestring = datetime.datetime.now().strftime("%y%m%d_%H%M")
        outfolder=path+'/'+datestring
        os.mkdir(outfolder)
        print('Directory change to ' + path)
        print('Outpath ' + outfolder)

if not flag:
    print('#### Error: Dirpath no existe')
    exit()
    
del flag, pathlist, p


# %% Lectura y conformación de dataframes

# Lectura de la count matrix con los pacientes positivos y negativos tomados de GSE152075.
df_conteo=pd.read_csv("GSE152075_raw_counts_GEO.txt", sep=' ', header=0)


# Generación de una matriz con las condiciones de COVID de los pacientes del estudio.
df_clase=pd.DataFrame(columns=['id', 'condicion'])
df_clase.id=df_conteo.columns
df_clase.condicion=df_clase.id.str.split('_').str[0].replace("POS", "Positivo").replace("NEG", "Negativo")
df_clase=df_clase.set_index('id').rename_axis(None)

# COmo necesito matriz paciente x gen se va a modificar el dataframe
df_conteo=df_conteo.T
# Con estos dos dataframes vamos a continuar trabajando.

# %% Exploración y filtrado de datos

# Aquellos genes que no tienen lecturas asignadas serán eliminados
# Aquellos que no tienen ninguna anotación respecto a la condición serán eliminados

#Si bien en el caso ejemplo el paso de eliminación por condición no es necesario, como es frecuente que ocurra y para futura utilidad del script se decide dejarlo.
npacientes=len(df_clase)
mtras = ~df_clase.condicion.isna()
print(f'De un total de {npacientes} pacientes se quedaron {len(mtras)}')
df_conteo = df_conteo.loc[mtras]
df_clase = df_clase.loc[mtras]


# Se borran aquellos genes que no tienen lecturas en ningún paciente
ngenes=len(df_conteo.columns)
genes_keep = df_conteo.columns[df_conteo.sum(axis=0) > 0]
print(f'De un total de {ngenes} genes se quedaron {len(genes_keep)}')
df_conteo=df_conteo[genes_keep]

# %% Expresión diferencial de genes

# A partir de los dos dataframes se realiza la expresión diferencial teniendo en cuenta la condición de los pacientes estudiados. En este caso: positivos y negativos.

# Primero lo que se debe hacer es pasar la información contenida a un tipo de archivo especial.DeseqDataSet, no se puede abrir con el gestor de archivos, se debe seguir trabajando sobre él. 
dds = DeseqDataSet(df_conteo, df_clase, design_factors="condicion", refit_cooks=True, n_cpus=2,)

# Se calculan algunos parámetros para la expresión diferencial. Se hacen varios pasos de fitting.
# Fitting dispersions, dispersion trend curve, MAP dispersions y LFCs
dds.deseq2()


# Se guarda el archivo intermedio con el único fin de poder reanudar el script desde este punto de ser necesario. 
print("Guardando el archivo dds.pkl")
with open(os.path.join(outfolder, "dds.pkl"), "wb") as f:
    pkl.dump(dds, f)

# %% Lectura de archivo dds.pkl 
# Fue generado en el punto anterior NO ES NECESARIO descomentarlo si el paso anterior no falló
file=outfolder+'/dds.pkl'
dds=pd.read_pickle(file, compression='infer')

# %% Análisis estadístico

# A partir del resultado anterior, lo que se calcula ahora son los p-valor y los p-valor ajustados para la expresion diferencial. Los resultados se van a guardar en el parámetro results_df del objeto.

estadistica= DeseqStats(dds, n_cpus=8)
estadistica.summary()

with open(os.path.join(outfolder, "estadistica.pkl"), "wb") as f:
    pkl.dump(estadistica, f)
    
print('res_estad.result_df')
df_summary=estadistica.results_df

filepath=outfolder+'/res_estad.result_df.txt'
df_summary.to_csv(filepath, sep='\t',index=False)
# %% Lectura de archivo estadistica.pkl 
# Fue generado en el punto anterior NO ES NECESARIO descomentarlo si el paso anterior no falló

# file=outfolder+"/estadistica.pkl"
# estadistica=pd.read_pickle(file, compression='infer')
# df_summary=estadistica.results_df

# file=outfolder+'/res_estad.result_df.txt'
# df_summary = pd.read_csv(file, sep="\t", header=0)
# df_summary=df_summary.set_index("Unnamed: 0").rename_axis(None)

# %% Evaluación de resultados.

# Para evaluar qué genes son importantes en este estudio, se hizo el análisis del pvalor y el pvalor ajustado. Aquellos significativamente relevantes serán los que tengan un pvalor ajustado menor a 0.01. 
# Esto va a ayudar tanto a tener resultados coherentes como para ayudar a reductir el ruido obtenido.

# Por otro lado, un parámetro que se calculó es logFoldChange. Se dice que un gen con un LFC que toma un valor de 2 está expresado 2 veces más en la condicion positiva respecto de la control

# Se van a eliminar aquellos genes en donde no tenga resultados de p valor ajustado.
ngenes=len(df_summary)
df_summary=df_summary.dropna(subset=['padj'])
print(f'Luego de la eliminación de valores NaN del padj, de un total de {ngenes} genes se quedaron {len(df_summary)}')

ngenes=len(df_summary)
df_summary=df_summary.dropna(subset=['pvalue'])
print(f'Luego de la eliminación de valores NaN del pvalor, de un total de {ngenes} genes se quedaron {len(df_summary)}')

# Este paso se dio probablemente porque el filtrado inicial fue conservador.

# %% Histograma de la distribución de los genes usando el pvalor ajustado y sin ajustar
## LUEGO PONERLO DENTRO DEL MISMO GRAFICO GRANDE
plotname=outfolder+'/histograma_pvalor.png'
fig = plt.figure()
plt.hist(df_summary.pvalue, 20, color="grey", histtype ='bar') 
plt.xlabel('pvalor')
plt.ylabel('Cantidad de apariciones')
plt.title('Distribución del valor del pvalor\n para el total de los genes analizados\n', fontweight ="bold")
plt.show()
fig.savefig(plotname, dpi=fig.dpi)

plotname=outfolder+'/histograma_pajustado.png'
fig = plt.figure()
plt.hist(df_summary.padj, 20, color="grey", histtype ='bar') 
plt.xlabel('Valor padj')
plt.ylabel('Cantidad de reads')
plt.title('Distribución del valor del padj \n para el total de los genes analizados\n', fontweight ="bold")
plt.show()
fig.savefig(plotname, dpi=fig.dpi)

# Se puede observar que la distribución es similar, pero no igual.

# %% Reajuste del LFC para que se pueda plotear.
# Lo que se hace es achivar el log2, es decir se trata de remover el ruido asociado con los bajos niveles de conteo de reads.
# El objetivo de esto es simplemente que se puedan visualizar mejor los datos.

# Los valores de LFC se guardan en el mismo lugar que los antereiores, son REEMPLAZADOS.
estadistica.lfc_shrink()

with open(os.path.join(outfolder, "shrunk_estadistica.pkl"), "wb") as f:
    pkl.dump(estadistica, f)

# %% Lectura de archivo dds.pkl 
# Fue generado en el punto anterior NO ES NECESARIO descomentarlo si el paso anterior no falló
file=outfolder+'/shrunk_estadistica.pkl'
estadistica=pd.read_pickle(file, compression='infer')

# %% Filtrado del archivo estadistica para poder hacer el volcano plot

df_shrink=estadistica.results_df

filepath=outfolder+'/datafinal.txt'
df_shrink.to_csv(filepath, sep='\t',index=False)

down = df_shrink[(df_shrink['log2FoldChange']<=-2)&(df_shrink['padj']<=0.01)].sort_values('padj')
up = df_shrink[(df_shrink['log2FoldChange']>=2)&(df_shrink['padj']<=0.01)].sort_values('padj')

# %% VolcanoPlot 1

plotname=outfolder+'/VolcanoPlot_1.png'
fig = plt.figure()
plt.scatter(x=df_shrink.log2FoldChange,y=df_shrink.padj.apply(lambda x:-np.log10(x)),s=1,label="Not significant",color="grey")

plt.scatter(x=down['log2FoldChange'],y=down['padj'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
plt.scatter(x=up['log2FoldChange'],y=up['padj'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

plt.xlabel("log2FoldChange")
plt.ylabel("-logFDR")
plt.axvline(-2,color="grey",linestyle="--")
plt.axvline(2,color="grey",linestyle="--")
plt.axhline(2,color="grey",linestyle="--")
plt.legend()
fig.savefig(plotname, dpi=fig.dpi)

# %% VolcanoPlot 2


plotname=outfolder+'/VolcanoPlot_2.png'
fig = plt.figure()
plt.scatter(x=df_shrink['log2FoldChange'],y=df_shrink['padj'].apply(lambda x:-np.log10(x)),s=1,label="Not significant",color="grey")
# highlight down- or up- regulated genes
plt.scatter(x=down['log2FoldChange'],y=down['padj'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
plt.scatter(x=up['log2FoldChange'],y=up['padj'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")
texts=[]
for i,r in up.iterrows():
    texts.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=i))
adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=0.1))
plt.xlabel("log2FoldChange")
plt.ylabel("-logFDR")
plt.axvline(-2,color="grey",linestyle="--")
plt.axvline(2,color="grey",linestyle="--")
plt.axhline(2,color="grey",linestyle="--")
plt.legend()
fig.savefig(plotname, dpi=fig.dpi)

# %% VolcanoPlot 3
df_plot=df_shrink[['index','log2FoldChange','padj']]
df_plot=df_plot.dropna(subset=['padj'])
df_plot.rename(columns={"index":'GeneNames', "log2FoldChange":'log2FC', "padj":'p-value'}, inplace=True)

filenames=outfolder+'/VolcanoPlot_3'
visuz.GeneExpression.volcano(df=df_plot, lfc='log2FC', pv='p-value', plotlegend=True, legendpos='upper right', legendanchor=(1.46,1), color=("#00239CFF", "grey", "#E10600FF"), valpha=0.5, geneid="GeneNames", gstyle=2, sign_line=True, figname=filenames, axtickfontname='Verdana', axlabelfontname='Verdana')


# %% VolcanoPlot 4
genesfiltrados=up['index'][0:10].tolist()+down['index'][0:10].tolist()

filenames=outfolder+'/VolcanoPlot_4'
visuz.GeneExpression.volcano(df=df_plot, lfc='log2FC', pv='p-value', plotlegend=True, legendpos='upper right', legendanchor=(1.46,1), color=("#00239CFF", "grey", "#E10600FF"), valpha=0.5, geneid="GeneNames",genenames=tuple(genesfiltrados), gfont=6, dotsize=4, sign_line=True, figname=filenames, axtickfontname='Verdana', axlabelfontname='Verdana', )

# %%

from sklearn.preprocessing import StandardScaler

# Separating out the features
x = df_conteo.columns.values

# Separating out the target
y = df.loc[:,['target']].values

# Standardizing the features
x = StandardScaler().fit_transform(x)
