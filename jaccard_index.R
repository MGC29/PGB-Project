setwd("~/Documents/pgb/project_neuroD2/chipseq_neuroD2/JACCARD index/all beds")
library(GenomicRanges)
install.packages(ComplexHeatmap)

library(ComplexHeatmap)
library(circlize)


load_bed <- function(filepath) {
  bed <- read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
  GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3))
}


# Lista de archivos BED
files <- c("ASCL1_.bed", "NEUROD1.bed", "NFYC.bed", 
           "MYOD1_mm10_.bed", "neurod2_peaks.coordinates.bed", 
           "TranscriptionFactors_bHLH_brain_mm.bed")

# Carga los archivos como GRanges
gr_list <- lapply(files, load_bed)
names(gr_list) <- gsub(".bed", "", files)  # Nombra las listas sin extensión

n <- length(gr_list)
jaccard_matrix <- matrix(0, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(gr_list)

for (i in 1:n) {
  for (j in 1:n) {
    intersection <- intersect(gr_list[[i]], gr_list[[j]])
    union <- union(gr_list[[i]], gr_list[[j]])
    jaccard_matrix[i, j] <- length(intersection) / length(union)
  }
}

heatmap(jaccard_matrix, 
        col = colorRampPalette(c("blue", "yellow", "red"))(100), 
        scale = "none", 
        main = "Jaccard Index Heatmap", 
        xlab = "TFs", 
        ylab = "TFs", 
        margins = c(10, 10))

# Renombrar los nombres de los factores
colnames(jaccard_matrix) <- c("ASCL1", "NeuroD1", "NeuroD2", "NFYC", "MYOD1", "TFs_Brain")
rownames(jaccard_matrix) <- c("ASCL1", "NeuroD1", "NeuroD2", "NFYC", "MYOD1", "TFs_Brain")

heatmap(jaccard_matrix, 
        col = colorRampPalette(c("green", "yellow", "red"))(100), 
        scale = "none", 
        main = "Jaccard Index Heatmap", 
        xlab = "TFs", 
        ylab = "TFs", 
        margins = c(10, 10))

heatmap(jaccard_matrix, 
        col = colorRampPalette(c("blue", "yellow", "red"))(100), 
        scale = "none", 
        main = "Jaccard Index Heatmap", 
        xlab = "TFs", 
        ylab = "TFs", 
        margins = c(10, 10), 
        cellnote = round(jaccard_matrix, 2))  # Valores redondeados a 2 decimales

# Identificar pares con overlap mayor al 1%
threshold <- 0.01
overlap_pairs <- which(jaccard_matrix > threshold, arr.ind = TRUE)

# Lista los pares con overlap significativo
significant_pairs <- data.frame(
  TF1 = rownames(jaccard_matrix)[overlap_pairs[, 1]],
  TF2 = colnames(jaccard_matrix)[overlap_pairs[, 2]],
  Jaccard_Index = jaccard_matrix[overlap_pairs]
)

print(significant_pairs)


# Filtrar la fila o columna correspondiente a NeuroD2
neurod2_overlap <- jaccard_matrix["NeuroD2", ]  # Extraer los índices de la fila de NeuroD2
print(neurod2_overlap)

# Definir umbral para el overlap
threshold <- 0.01  # 1%

# TFs con overlap significativo
significant_overlap <- neurod2_overlap[neurod2_overlap > threshold]
print(significant_overlap)

# Gráfico de barras del overlap con NeuroD2
barplot(
  neurod2_overlap,
  main = "Jaccard Index with NeuroD2",
  xlab = "TFs",
  ylab = "Jaccard Index",
  col = ifelse(neurod2_overlap > threshold, "red", "blue"),
  las = 2,  # Rotar nombres en el eje x
  cex.names = 0.8  # Reducir tamaño de las etiquetas
)
abline(h = threshold, col = "green", lty = 2)  # Línea para el umbral

write.csv(data.frame(TF = names(significant_overlap), Jaccard_Index = significant_overlap), 
          "NeuroD2_overlap.csv", 
          row.names = FALSE)

#########################################################
# Definir los TFs que deseas comparar
tf_subset <- c("ASCL1", "NeuroD1", "NFYC", "MYOD1")

# Filtrar los índices de Jaccard para esos TFs contra NeuroD2
neurod2_subset <- neurod2_overlap[tf_subset]
print(neurod2_subset)
#######################################
# Gráfico de barras para los TFs seleccionados
barplot(
  neurod2_subset,
  main = "Jaccard Index: NeuroD2 vs Selected TFs",
  xlab = "TFs",
  ylab = "Jaccard Index",
  col = ifelse(neurod2_subset > threshold, "red", "blue"),  # Colores según el umbral
  las = 2,  # Rotar nombres en el eje x
  cex.names = 0.8  # Reducir tamaño de las etiquetas
)
abline(h = threshold, col = "green", lty = 2)  # Línea para el umbral

# Cambiar los nombres del vector para que sean personalizados
names(neurod2_subset) <- c("ASCL1", "NeuroD1", "NFYC", "MYOD1")

# Vuelve a graficar con los nombres personalizados
barplot(
  neurod2_subset,
  main = "Jaccard Index: NeuroD2 vs Selected TFs",
  xlab = "TFs",
  ylab = "Jaccard Index",
  col = ifelse(neurod2_subset > threshold, "red", "blue"),
  las = 2,
  cex.names = 0.8
)
abline(h = threshold, col = "green", lty = 2)

write.csv(data.frame(TF = names(neurod2_subset), Jaccard_Index = neurod2_subset), 
          "NeuroD2_vs_Selected_TFs.csv", 
          row.names = FALSE)


# Define los TFs a comparar
tf_subset <- c("ASCL1", "NeuroD1", "NFYC", "MYOD1", "NeuroD2")  # Incluye NeuroD2

# Filtra la matriz de Jaccard para estos TFs
jaccard_filtered <- jaccard_matrix[tf_subset, tf_subset]
print(jaccard_filtered)


heatmap(
  jaccard_filtered,
  col = colorRampPalette(c("blue", "yellow", "red"))(100), 
  scale = "none", 
  main = "Jaccard Index Heatmap: NeuroD2 and Selected TFs",
  xlab = "TFs", 
  ylab = "TFs", 
  margins = c(10, 10)  # Ajustar márgenes para etiquetas largas
)


# Renombra las filas y columnas
rownames(jaccard_filtered) <- c("ASCL1", "NeuroD1", "NFYC", "MYOD1", "NeuroD2")
colnames(jaccard_filtered) <- c("ASCL1", "NeuroD1", "NFYC", "MYOD1", "NeuroD2")

# Vuelve a graficar el heatmap
heatmap(
  jaccard_filtered,
  col = colorRampPalette(c("blue", "yellow", "red"))(100), 
  scale = "none", 
  main = "Jaccard Index Heatmap: NeuroD2 and Selected TFs",
  xlab = "TFs", 
  ylab = "TFs", 
  margins = c(10, 10)
)

install.packages("pheatmap")
library(pheatmap)

# Crear el heatmap
pheatmap(
  jaccard_filtered,
  color = colorRampPalette(c("blue", "yellow", "red"))(100),
  main = "Jaccard Index Heatmap: NeuroD2 and Selected TFs",
  cluster_rows = FALSE,  # Sin clustering si no es necesario
  cluster_cols = FALSE
)





