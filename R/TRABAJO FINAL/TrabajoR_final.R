##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos<-read.table("archivo-entrada.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
head(datos)
tail(datos)
dim(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

control_media<-apply(datos[,2:4],1,mean) #a cada fila se hace la media de las columnas 2 a la 4 (control)
condi_media<-apply(datos[,5:7],1,mean) #a cada fila se hace la media de las columnas 5 a la 7 (condicion)
control_de<-apply(datos[,2:4],1,sd) #a cada fila se hace la desviacion estandar de las columnas 2 a la 4 (control)
condi_de<-apply(datos[,5:7],1,sd) #a cada fila se hace la desviacion estandar de las columnas 5 a la 7 (condicion)
datos<-cbind(datos, control_media, condi_media, control_de, condi_de)
dim(datos) #para comprobar que se me han añadido las 4 columnas.
names(datos) #para comprobar que las columnas tiene los nombres adecuados.

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

control_media2<-log2(control_media)  #pongo log2 para que los datos pequeños sean más fáciles de representar
control_de2<-log2(control_de)
condi_media2<-log2(condi_media)
condi_de2<-log2(condi_de)
png("Histograma control-condi.png")
par(mfrow = c(1, 2)) #en una imagen (fila), dos histogramas (columnas).
hist(control_media2, breaks = 15, col = "lightgreen",
     main = "Media control", xlab = "Media", ylab = "Frecuencia")
hist(condi_media2, breaks = 15, col = "purple",
     main = "Media condicion", xlab = "Media", ylab = "Frecuencia")
par(mfrow = c(1, 1))
dev.off()

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

png("Boxplot control-condi.png")
boxplot(control_media2, condi_media2,
        names = c("Control", "Condición"),
        col = c("lightgreen", "purple"),
        main = "Medias control-condicion",
        ylab = "Media")
legend("topleft", legend = c("Control", "Condición"), fill = c("lightgreen", "purple"))
dev.off()

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(unlist(datos[, 2:4])) #shapiro test para las columnas del control: p valor<0,05
shapiro.test(unlist(datos[, 5:7])) #shapiro test para las columnas de condicion: p valor<0,05
#como el pvalor es < 0,05, aplico wilcox.test porque los datos no son normales.
pvalues <- numeric(nrow(datos))  #crea una lista vacia con tantas filas como filas tiene el conjutno de datos (619)
for (i in 1:nrow(datos)) {  #i es cada valor/cada fila
  control <- as.numeric(datos[i, c("C2-TPM", "C1-TPM", "C3-TPM")])  #genera un vector con los 3 datos de control
  condicion <- as.numeric(datos[i, c("1W5-TPM", "1W6-TPM", "1W5-COV")])  #genera un vector con los 3 datos de condicion
  pvalues[i] <- wilcox.test(  #compara el vector generado por cada fila de control-condicion
    control,
    condicion,
    exact = FALSE  #fuerza una aproximacion normal para evitar que el pvalor se quede en 1 decimal
  )$p.value
}

datos$pvalue_wilcox <- pvalues
head(datos)
dim(datos)


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

exp_dif <- datos[!is.na(datos$pvalue_wilcox) & datos$pvalue_wilcox < 0.05, ]
nrow(exp_dif) #salen 0 genes con expresion diferencial con pvalue<0,05 asi que lo cambio a pvalor<0,1.
exp_dif <- datos[!is.na(datos$pvalue_wilcox) & datos$pvalue_wilcox < 0.1, ]
nrow(exp_dif)  #salen 468 genes diferencialmente expresados.
exp_dif[, 1]#listado de genes con expresion diferencial.
exp_dif
sobre <- exp_dif$condi_media>exp_dif$control_media  #la tabla exp_dif no tiene nombres de columnas, por lo que puedo usar directamente la formula para sacar la media de cada grupo.
sobre
sob<-exp_dif[which(sobre), 1] #cuenta los TRUE, que son los genes que cumplen la condicion de que la media de la condicion>media de control
length(sob)  #salen 69 genes sobreexpresados.

reprimidos <- exp_dif$condi_media<exp_dif$control_media
reprimidos
rep<-exp_dif[which(reprimidos), 1] #cuenta los TRUE, que son los genes que cumplen la condicion de que la media de la condicion<media de control
length(rep)  #salen 399 genes reprimidos.

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

mat <- as.matrix(exp_dif[, 2:7])  #paso exp_dif a una matriz para que pueda ser usada para hacer el heatmap.
class(mat)  #compruebo que se ha convertido en matriz.
rowv = as.dendrogram(hclust(as.dist(1-cor(t(mat)))))
#as.dendogram permite organizar los datos de forma jerarquizada para que se pueden representar como heatmap.
#1-cor correlaciona columnas cercanas entre sí, pero como queremos correlacionar filas (genes) hay que trasponer las filas por las columnas usando t(filtered).
#hclust agrupa las filas mas parecidas y as.dist permite que se genere un formato de distancias.
colv = as.dendrogram(hclust(as.dist(1-cor(mat)))) #aqui no usamos t(filtered) porque queremos agrupar las columnas.
png("Heatmap genes con expresion diferencial.png")
heatmap(mat, Rowv=rowv, Colv=colv, cexCol=0.7,
        col = rev(redblue(256)), scale = "row")
dev.off()
