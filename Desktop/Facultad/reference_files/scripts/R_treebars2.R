#!/usr/bin/env Rscript

install.packages("dplyr")
install.packages("dendextend")
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

library("DECIPHER")   
library('dendextend')
library('dplyr')

## Quitar al momento de integrarlo al código
setwd("/home/jpereira/Escritorio/respaldo.Inflation_45/refeed_Toxo_Neo/R_trees/") 

##### Dataframe y arbol con genes seleccionados (anotados SRS)
#df_cluster <- read.csv2( "df_filtered_selection_2.0", header = T, sep = "\t", dec = "." )
#tree_CL <- ReadDendrogram('mafft_all_selection.treefile')

df_cluster <- read.csv2( "df_filtered_whole_2.0", header = T, sep = "\t", dec = "." )
tree_CL <- ReadDendrogram('mafft_all_whole.treefile')

labels_CL <- labels(tree_CL)
labels_CL <- gsub(" ", "_", labels_CL)
df_CL <- df_cluster %>% slice(match(labels_CL, df_cluster$genID))

##### Se asocia a un intervalo de valores un gradiente de colores 
common_vals_gradient <- function( values_list, 
                                  cols = c( colors()[109], colors()[30]),
                                  steps = 20,
                                  outlier = 1/5) # Diferencia equivalente a 1/5 de la longitud 
  # del intervalo de valores comunes
{
  gradient <- colorRampPalette( c(cols[1], cols[2]) )(steps)              # Gradiente de colores
  common_vals <- quantile(values_list, probs = c(0.025,0.975)) %>% unname 
  range_val <- common_vals[2] - common_vals[1]                            # Rango de valores comunes
  
  out_cols <- c()
  for ( value in values_list ){
    # Valores menores al 2.5% de los datos
    if ( value < common_vals[1] ) {
      if ( (common_vals[1] - value) > range_val*outlier ){ ## Se determina si el valor ees un outlier
        index <- 1
      } else {
        index <- 1
      }
    } 
    # Valores mayores al 97.5% de los datos
    if ( value > common_vals[2] ) {
      #print(20)
      if ( (value - common_vals[2] ) > range_val*outlier ){ ## Se determina si el valor ees un outlier
        index <- steps 
      } else {
        index <- steps
      }
    }
    # Resto de los valores, correspondiente al 95% mas comúm   
    if ( (common_vals[1] < value) & (value < common_vals[2]) ) {
      index <- ((value - common_vals[1])/(range_val/steps)) %>% ceiling()  ## Rango entre %2.5 - %97.5
    }
    out_cols <- append(out_cols, gradient[index])
  }
  return(out_cols)
}
## Colores para todos los clusters. colclust_2 

colGC_prom <- common_vals_gradient(values_list = df_CL$prom_cg)
colGC_genes <- common_vals_gradient(values_list = df_CL$genes_cg)
colLongs_prom <- common_vals_gradient(values_list = df_CL$prom_longs, cols = c("#FFC685", "#9E3D22"))
colLongs_genes <- common_vals_gradient(values_list = df_CL$genes_longs, cols = c("#FFC685", "#9E3D22"))
colclust_2 <- c(); {
  row_index <- as.numeric(rownames(df_CL))
  all_clust <- unique(df_CL$clust_num)
  for ( row in row_index ){
    if (df_CL[row,]$clust_size >= 5){
      for ( col in 1:length(all_clust) )
        if ( all_clust[col] == df_CL[row,]$clust_num ) {
          colclust_2 <- append(colclust_2, col)
        } 
    } else {
      colclust_2 <- append(colclust_2, "black")  
    }
  }}
colorgs <- c(); {org_names <- df_CL$Org_name %>% unique
smol_cols <- c(rev(cols34)[9], rev(cols34)[10])
 {
  row_index <- as.numeric(rownames(df_CL))
  for ( row in row_index ) {
    for (indx in 1:length(org_names)){
      if ( org_names[indx] == df_CL[row,]$Org_name ) {
        colorgs <- append(colorgs, smol_cols[indx])
      }
    }
  }}}

##############3

summary(df_CL0$product)

##########3
# Sustituir con  NAs los "-----".
df_CL0 <- df_CL
df_CL0["product"][df_CL0["product"] == "-----"] <- NA 

df_CL0$product <- as.character(df_CL0$product) 
df_CL0$description <- as.character(df_CL0$description)

df_CL0 <- within(df_CL0, product <- ifelse(is.na(product), description, product))
df_CL0$product

df_CL0$product <- as.factor(df_CL0$product) 
df_CL0$description <- as.factor(df_CL0$description)

## Homogeneniso los factores del DF usando Regexp
df_CL0$product <- gsub('.*SAG.*', 'SRS', df_CL0$product)
df_CL0$product <- gsub('.*SRS.*', 'SRS', df_CL0$product)
df_CL0$product <- gsub('.*srs.*', 'SRS', df_CL0$product)
df_CL0$product <- gsub('hypothetical protein.*', 'hypothetical protein', df_CL0$product)

## Lista colores asociado a cada producto del gen 
colprod <- c(); {prod_names <- df_CL0$product %>% unique
#smol_cols <- c(rev(cols34)[9], rev(cols34)[10])
{
  row_index <- as.numeric(rownames(df_CL0))
  for ( row in row_index ) {
    for (indx in 1:length(prod_names)){
      if ( prod_names[indx] == df_CL0[row,]$product ) {
        colprod <- append(colprod, indx)
      }
    }
  }}}; 


col_bars <- cbind( colLongs_genes, colLongs_prom, colGC_genes, colGC_prom, colorgs,colprod,colclust_2)
col_labels <- c( "", "Size", "", "%GC" ,"Org", "product","Cluster")

eP <- list(lwd = 1, t.cex = 0.3, p.col = "plum")
nP <- list(pch = NA, lab.cex = 0.5) 

tree_CL2 <- tree_CL %>% 
  set( "labels", c(1:length(labels(tree_CL))) ) #%>%
#set("hang_leaves", 0)

## 34 distinct colors
cols34 <- c("darkgrey","#2f4f4f","#7f0000","#808000","#483d8b","#008000",
            "#3cb371","#4682b4","#000080","#9acd32","#20b2aa","#32cd32",
            "#8b008b","#b03060","#ff4500","#ff8c00","#7fff00","#00fa9a",
            "#8a2be2","#dc143c","#f4a460","#0000ff","#da70d6","#ff00ff",
            "#1e90ff","#f0e68c","#fa8072","#ffff54","#ff1493","#7b68ee",
            "#afeeee","#7fffd4","#ffe4c4","#ffb6c1")

palette(rev(cols34))

pdf(file="tree_AW_colbars_test.pdf", height=15, width=30)
plot(tree_CL2, nodePar = nP, edgePar = eP)#, xaxt="n") #, yaxt="n")
colored_bars(colors = col_bars,
             dend = tree_CL2,
             sort_by_labels_order = FALSE,
             y_scale = 0.375,
             y_shift = -0.1,
             rowLabels = col_labels )
dev.off()

############ HANG TREE ############
tree_CL2 <- tree_CL %>% 
  set( "labels", c(1:length(labels(tree_CL))) ) %>%
  set("hang_leaves", -1)
  
pdf(file="tree_AW_colbars_hang_test.pdf", height=15, width=30)
plot(tree_CL2, nodePar = nP, edgePar = eP)#, xaxt="n") #, yaxt="n")
colored_bars(colors = col_bars,
             dend = tree_CL2,
             sort_by_labels_order = FALSE,
             y_scale = 0.375,
             y_shift = -0.1,
             rowLabels = col_labels )
dev.off()
########### SUB DENDROGRAMS ###################

# Algoritmo de visualización
# El ancho de las barras en total debe ser 1/10 del total del gráfico
# Deben de estar ubicadas levemente por debajo del 0
# Tomar la altura de la hoja mas baja, bajar todos los nodos a dicha altura 
# en caso que se muestren solo los números

{
#Difucultades de automatización
# You will often needs to adjust the y_scale, y_shift and the text_shift
# parameters, in order to get the bars in the location you would want
# This can probably be done automatically, but will require more work.
# Since it has to do with the current mar settings, the number of groups,
# and each computer's specific graphic device. patches for smarter defaults 
# will be appreciated

## Grupos de genes de genes monofiléticos pertenecientes a un solo cluster del MCL
A1 <- c(3:7)
A2 <- c(8:89)
A3 <- c(93:110) #93:108
A4 <- c(111:117) 
A5 <- c(130:136) 
A6 <- c(137:141)
A7 <- c(142:171)
A8 <- c(172:180)
A9 <- c(193:199)
A10 <- c(202:206)
A11 <- c(212:218)
A12 <- c(219:223) ## Grafico 12 no se aprecian bien las barras
A13 <- c(258:309) #A13-1 238:303, #A13-2, 304:309
A14 <- c(315:338)

######################
A_list <- list(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14)
i=0
for (A in A_list) {
  i=i+1
  tree_A <- find_dendrogram(tree_CL2, A) 
  if (is.null(tree_A)) {
    print(paste("Null en árbol:", i))
    break
  } 
  
  ## Se obtiene la altura mínima de cada árbol y se la sustrae a todos los nodos 
  ## de forma de acercar el dednrograma al eje x
  h <- tree_A %>% get_nodes_attr("height") %>% min()
  tree_Ah <- tree_A %>% raise.dendrogram(-h)
  leaves_H <- function(n, h){ 
    if (is.leaf(n)) {
      attributes(n)$height  <- attributes(n)$height-h  }
    n
  } ## Se aplica dendrapply solamente a las hojas para evitar excesiva recursividad
  tree_Ah <- dendrapply(tree_Ah, leaves_H, h ) 
  
  col_bars_A <- cbind( colLongs_genes, colGC_genes, colprod,colorgs)
  col_labels_A <- c(  "Size",  "%GC" ,"Product", "Org")
  
  ## Ajustes 'a priori' del tamaño de la imagen, y elementos del dendrograma
  ## según el número total de hojas en el dendrograma
  fig_H = (length(A))**(1/2) + 3
  fig_W = fig_H*1.5
  pdf(file=paste("tree_A",i,"_colbars_test_product.pdf"), height= fig_H, width=fig_W)
  plot(tree_Ah, nodePar = nP, edgePar = eP, yaxt="n")#, xaxt="n") #, yaxt="n")
  colored_bars(colors = col_bars_A[min(A):max(A),],
               dend = tree_Ah,
               sort_by_labels_order = FALSE,
               y_scale = 1/fig_H*1.2-0.05,
               y_shift = -0.1,
               text_shift = 40/fig_H - 1.5,
               rowLabels = col_labels_A )
  dev.off()
  print(paste("Finalizando árbol: ", i))
}
}  


#######################################################################################
A <- A1
{tree_A <- find_dendrogram(tree_CL2, A) 
## Ajusto altura del arbol y las hojas 
h <- tree_A %>% get_nodes_attr("height") %>% min()
tree_Ah <- tree_A %>% raise.dendrogram(-h)
leaves_H <- function(n, h){ 
  if (is.leaf(n)) {
    attributes(n)$height  <- attributes(n)$height-h  }
  n
} ## Se aplica dendrapply solamente a las hojas para evitar excesiva recursividad
tree_Ah <- dendrapply(tree_Ah, leaves_H, h ) 

col_bars <- cbind( colLongs_genes, colGC_genes, colorgs)
col_labels <- c(  "Size",  "%GC" ,"Org")


fig_H = (length(A))**(1/2) + 3
fig_W = fig_H*1.5
pdf(file="tree_A2_colbars_test.pdf", height= fig_H, width=fig_W)
plot(tree_Ah, nodePar = nP, edgePar = eP, yaxt="n")#, xaxt="n") #, yaxt="n")
colored_bars(colors = col_bars[min(A):max(A),],
             dend = tree_Ah,
             sort_by_labels_order = FALSE,
             y_scale = 1/fig_H*1.2 - 0.05,
             y_shift = -0.1,
             text_shift = 40/fig_H -1.5,
             rowLabels = col_labels )
dev.off()
}

?colored_bars
tree_A %>% plot( nodePar = nP, edgePar = eP )

A <- list(c(1,2))
B <- list(c(3,4,5))
l <- list(A,B)

for (i in 1:length(l)){
  print("Bucle 1")
  print(l[[i]])
  for (j in 1:length(l[[i]])) {
    print("Bucle 2")
    print(j)
  }
}

#####################################################
x <- c(1,2,3,4)
y <- c(1,2,3,4)


df_CL[,c(4,6)]
?replace
df_CL[c(31,173,202),]
