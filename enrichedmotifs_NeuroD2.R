setwd("~/Documents/pgb/project_neuroD2/chipseq_neuroD2/dataset_neurod2/promotersND2/resultadopromoters")

library(dplyr)


#######################################################
#ANALISIS WITH FIMO MEME#
# Leer el archivo
fimo_results_nd2_promoters <- read.table("fimo.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Crear una clave para coordenadas únicas en el archivo de FIMO
fimo_results_nd2_promoters <- fimo_results_nd2_promoters %>%
  mutate(coord_key = paste(sequence_name, start, stop, sep = "_"))
########################################################


#######################################################
#ANALISIS WITH THE PROMOTERS,ANNOTATED GENES FROMO PEAKS NED2#
# Cargar los datos de picos con genes asociados
promotersnd2_genes <- read.table("promotersnd2_genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(promotersnd2_genes) <- c("chr", "start", "end", "peak_id", "intensity", "region", 
                     "gene_chr", "gene_start", "gene_end", "gene_info")
head(promotersnd2_genes)
# Extraer el nombre del gen de la columna "gene_info"
promotersnd2_genes$gene_name <- sub(".*gene_name=([^;]+);.*", "\\1", promotersnd2_genes$gene_info)

# Crear una clave para coordenadas únicas en el archivo de picos
promotersnd2_genes <- promotersnd2_genes %>%
  mutate(coord_key = paste(chr, start, end, sep = "_"))

########################################################


###################################################
#unir datos usando solapamientos
library(GenomicRanges)

# Convertir fimo_results_nd2_promoters a GRanges
fimo_gr <- makeGRangesFromDataFrame(fimo_results_nd2_promoters,
                                    seqnames.field = "sequence_name",
                                    start.field = "start",
                                    end.field = "stop")

# Convertir promotersnd2_genes a GRanges
promoters_gr <- makeGRangesFromDataFrame(promotersnd2_genes,
                                         seqnames.field = "chr",
                                         start.field = "start",
                                         end.field = "end")

# Buscar solapamientos
overlap_indices <- findOverlaps(fimo_gr, promoters_gr)

# Crear una tabla combinada con los solapamientos
combined_data_nd2 <- as.data.frame(cbind(
  fimo_results_nd2_promoters[queryHits(overlap_indices), ],
  promotersnd2_genes[subjectHits(overlap_indices), ]
))

# Ver los primeros resultados
head(combined_data_nd2)

colnames(combined_data_nd2)


library(tibble)

# Reparar nombres duplicados automáticamente
combined_data_nd2 <- as_tibble(combined_data_nd2, .name_repair = "unique")

# Ver los primeros resultados
head(combined_data_nd2)

colnames(combined_data_nd2)

##################################################################
#Hacer las gràficas 
library(ggplot2)
library(dplyr)

####################################################################
#########1_Motivos más Frecuentes y Genes Asociados
##Shows the most common motifs and to which genes they are associated.
####################################################################

# Contar la frecuencia de motivos por gen
motif_gene_counts <- combined_data_nd2 %>%
  group_by(motif_id, gene_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

# Gráfica de barras
ggplot(motif_gene_counts, aes(x = reorder(motif_id, -count), y = count, fill = gene_name)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Frequency of motifs found for NeuroD2 TF and associated genes.",
       x = "Motifs",
       y = "Frequency",
       fill = "Gene") +
  theme_minimal()

###########10Motifs mas frecuentes
# Filtrar los 10 motivos más frecuentes
top_motifs <- motif_gene_counts %>%
  top_n(2, count)

#######5motifs more frequents
# Paso 2: Filtrar los 5 motivos más frecuentes
top_motifs <- motif_gene_counts %>%
  group_by(motif_id) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 3) # Tomar los 5 motivos más frecuentes

# Paso 3: Seleccionar los 3 genes más frecuentes por motivo
top_motifs_genes <- motif_gene_counts %>%
  filter(motif_id %in% top_motifs$motif_id) %>%
  group_by(motif_id) %>%
  slice_max(order_by = count, n = 3) %>%
  ungroup()

# Paso 4: Crear la gráfica de barras
ggplot(top_motifs_genes, aes(x = reorder(motif_id, -count), y = count, fill = gene_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Frequency of motifs found for NeuroD2 TF and associated genes",
       x = "Motifs",
       y = "Frequency",
       fill = "Gene") +
  theme_minimal()
################################################################################################
#biological pathway analysis
gene_list <- unique(combined_data_nd2$gene_name)

gene_list <- unique(combined_data_nd2$gene_info) # Adjust based on your column

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db") # For mouse genes
library(clusterProfiler)
library(org.Mm.eg.db)

##########################################################


####################################################################
##2_Intensidad Promedio de Picos por Motivo
##Identifies which motifs are associated with stronger peaks.
####################################################################


library(ggplot2)
library(dplyr)

# Calcular intensidad promedio de picos por motivo
motif_peak_intensity <- combined_data_nd2 %>%
  group_by(motif_id) %>%
  summarise(avg_intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_intensity)) # Ordenar por intensidad promedio

#Filtrar los 10 motivos más relevantes
top_motifs <- motif_peak_intensity %>%
  slice_head(n = 3) # Tomar los primeros 10

# Crear la gráfica de barras
ggplot(top_motifs, aes(x = reorder(motif_id, -avg_intensity), y = avg_intensity)) +
  geom_bar(stat = "identity", fill = "steelblue") + # Barras con color personalizado
  geom_text(aes(label = round(avg_intensity, 2)), hjust = -0.150) + # Etiquetas con los valores
  coord_flip() + # Rotar los ejes para mayor legibilidad
  labs(title = "Motifs by Peak Average Intensity",
       x = "Motif",
       y = "Peak Average Intensity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Mejorar etiquetas
  scale_fill_brewer(palette = "Blues") # Colores personalizados opcionales


####################################################################
#### 3_Relación entre la Puntuación del Motivo y la Intensidad del Pico
##Explores whether high scoring motifs correspond to intense peaks
####################################################################

# Gráfica de dispersión
library(ggplot2)
library(dplyr)

# Seleccionar los 10 genes con más ocurrencias
top_genes <- combined_data_nd2 %>%
  count(gene_name, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(gene_name)

# Filtrar los datos para incluir solo esos genes
filtered_data <- combined_data_nd2 %>%
  filter(gene_name %in% top_genes)

# Crear la gráfica con los datos filtrados
ggplot(filtered_data, aes(x = score, y = intensity, color = gene_name)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(
    title = "Relationship between motif score and peak intensity (Top 10 Genes)",
    x = "Motif Score",
    y = "Peak Intensity",
    color = "Gen Asociado"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#######
####################################################################
#### 4_Frecuencia de Regiones Promotoras Asociadas a Motivos
##Determines whether motifs are enriched in specific promoter regions
####################################################################

# Count frequency of promoter region types
region_counts <- combined_data_nd2 %>%
  group_by(region) %>%
  summarise(count = n(), .groups = "drop")

# Bar plot
ggplot(region_counts, aes(x = reorder(region, -count), y = count, fill = region)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Frequency of Promoter Region Types",
    x = "Promoter Region Type",
    y = "Frequency"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

####################################################################
#### 5_Genes Más Enriquecidos en Motivos
##Identifies priority genes with higher numbers of associated motifs
####################################################################

library(ggplot2)
library(dplyr)

# Paso 1: Extraer la raíz del nombre del gen
gene_motif_counts <- combined_data_nd2 %>%
  mutate(gene_root = gsub("\\d+$", "", gene_name)) %>% # Eliminar números al final del nombre
  group_by(gene_name) %>%
  summarise(count = n(), gene_root = unique(gene_root), .groups = "drop")

# Paso 2: Seleccionar el gen con más motivos por raíz
filtered_genes <- gene_motif_counts %>%
  group_by(gene_root) %>%
  slice_max(order_by = count, n = 1) %>% # Seleccionar solo 1 gen por raíz
  ungroup()

# Paso 3: Seleccionar los 10 genes con más motivos entre los filtrados
top_10_genes <- filtered_genes %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) # Limitar a los 10 genes con más motivos

# Paso 4: Crear la gráfica
ggplot(top_10_genes, aes(x = reorder(gene_name, -count), y = count)) +
  geom_bar(stat = "identity", fill = "salmon") +
  geom_text(aes(label = count), hjust = -0.2, size = 3) + # Etiquetas encima de las barras
  coord_flip() +
  labs(
    title = "The genes most enriched in motifs from NeuroD2  ",
    x = "Gene",
    y = "Number of Associated Motifs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10), # Ajustar tamaño de las etiquetas de los ejes
    axis.title = element_text(size = 12)
  )




####MOTIF CENTRALITY 

# Calcular el centro del pico y la distancia al motivo
combined_data_nd2 <- combined_data_nd2 %>%
  mutate(
    peak_center = (start...13 + end) / 2, # Usar las columnas correctas
    motif_center = (start...4 + stop) / 2,
    distance_to_center = abs(motif_center - peak_center)
  )

ggplot(combined_data_nd2, aes(x = distance_to_center)) +
  geom_histogram(binwidth = 50, fill = "lightblue", alpha = 0.7) +
  labs(
    title = "Motif Centrality Distribution",
    x = "Distance to Peak Center",
    y = "Frequency"
  ) +
  theme_minimal()

####Binding Preference in Genomic Elements
region_counts <- combined_data_nd2 %>%
  group_by(region) %>%
  summarise(count = n(), .groups = "drop")

ggplot(region_counts, aes(x = reorder(region, -count), y = count, fill = region)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Binding Preference in Genomic Elements",
    x = "Genomic Element",
    y = "Frequency"
  ) +
  theme_minimal()

#########################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")

library(AnnotationHub)

# Query AnnotationHub for NeuroD2 domains
ah <- AnnotationHub()
domains <- query(ah, "NeuroD")
domains



# Combine motif data with peaks
combined_data_nd2 <- combined_data_nd2 %>%
  mutate(domain = ifelse(motif_id == "bHLH motif", "bHLH Domain", "Other Domain"))

library(ggplot2)

ggplot(combined_data_nd2, aes(x = domain, fill = domain)) +
  geom_bar() +
  labs(title = "Domain-Specific Binding of NeuroD2",
       x = "Domain",
       y = "Number of Peaks") +
  theme_minimal()
############################################################################


