#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np

############## Se produce una matriz para todos los clusters de la busqueda ##############################

#Se crea el archivo order_out_mci que contiene todos los genes en una sola linea
#este conserva el orden original de los clusters, de mayor a menor
#get_ipython().system('tr "\\n" "\\t" < data_out_mci > order_out_mci')
#get_ipython().system("sed -i 's/\\t$//' order_out_mci")

#Se leen los archivos usando pandas
genes1 = pd.read_csv('order_out_mci', sep='\t')
data1 = pd.read_csv('evalue_filter2', sep='\t', lineterminator='\n', header=None)

#data.describe()
#genes.describe()

#Se crea una matriz vacia que respeta el orden original de los genes
cols1=rows1=[x for x in genes1]
matrix1=pd.DataFrame(-1, columns=cols1, index=rows1, dtype='float')

#Se rellena la matriz con los valores de la tercer columna de data1
for i_index, i_row in data1.iterrows():
  valor_actual = matrix1[i_row[0]][i_row[1]] 
  if valor_actual == -1 or valor_actual > i_row[2]: 
    matrix1[i_row[0]][i_row[1]] = i_row[2]

#Se transforman los valores de la matris usando -log10(x + 1e-200)
def minus_log(x = float):
    if x != -1:
        return(-(np.log(x+10**(-200))/np.log(10)))
    else:
        return(0)
log_matrix1=matrix1.applymap(minus_log)

#########################################################################################################

import seaborn as sns
import matplotlib.pyplot as plt

######## Se crea el heatmap a partir de la matriz que contiene todolos clusters de la busqueda ##########
sns.set(rc={'figure.figsize':(100,100)})

heat1 = sns.heatmap(log_matrix1, annot=False)

fig = heat1.get_figure()
fig.savefig("All_clusters_heat_map.png") 
plt.close()

###################################### Se obtiene una lista de indices #########################################

index_list_temp=open('data_index').read().split('\n')
ind_list=[]
for i in index_list_temp:
    if i != '':
        ind_list.append(int(i))


###########################   Se realizan dendrogramas+heatmaps para cada cluster   ############################

#Obtengo el archivo 'gene_number' indicando la cantidad de genes por cluster:
#get_ipython().system("awk '{print NF}' data_out_mci > gene_number")
#get_ipython().system('tr "\\n" "\\t" < gene_number > A')
#get_ipython().system('mv A gene_number')

#Guardo la cantidad en la lista 'gene_num' la cantidad dde genes por cluster
gene_num=open('gene_number')
gene_num=list(gene_num)[0].split('\t')
if gene_num[len(gene_num)-1] == '':
    gene_num.pop(len(gene_num)-1)

#A partir del df log_matrix1, selecciono las submatrices correspondientes a
#cada cluster y elaboro un dendograma
prev_num=0
num=-1
index=-1
for number in gene_num:
    num+=int(number)
    index+=1
    #print(num, prev_num)
    if num +1-prev_num > 1:
        print("Creando drendrograma para el cluster: "+str(ind_list[index]))
        sub_log_matrix=log_matrix1.iloc[prev_num:num+1,prev_num:num+1]
        prev_num=num+1
        sns.clustermap(sub_log_matrix, xticklabels=False, yticklabels=True,
                       cbar_pos=(0.8, 0.8, 0.05, 0.18) ,square=True,
                       dendrogram_ratio=(0.00001,0.2))
        #plt.show()
        plt.savefig('Cluster_'+str(ind_list[index])+'_Dendrogram')
        plt.close()
    else:
        print("Salteandose dendrograma para el cluster pequeño Nº: "+ str(ind_list[index]))    

    #plt.imshow(sub_log_matrix, cmap='hot', interpolation='nearest')
    #plt.show()


### Cosas que quedan por hacer:
### -Detectar y descartar outliers (Opcionalmente)
### -Obtener una matriz simetrica
### -Marcar de ser posible la seleccion
### -Marcar con un mapa de colores cada cluster 



