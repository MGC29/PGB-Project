setwd("~/Documents/pgb/project_neuroD2/chipseq_neuroD2")
###########################################################
#https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
################
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.20")
install.packages("systemfonts")
install.packages("ggforce")

install.packages(c("systemfonts", "ggforce", "scatterpie"))
BiocManager::install("enrichplot")

BiocManager::install("ChIPseeker")
# Si también necesitas el paquete GenomicFeatures para cargar el archivo GFF
if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
  BiocManager::install("GenomicFeatures")
}
BiocManager::install("txdbmaker")

# Instalar txdbmaker si no está instalado
if (!requireNamespace("txdbmaker", quietly = TRUE)) {
  install.packages("txdbmaker")
}

# Cargar txdbmaker
library(txdbmaker)

library(GenomicFeatures)
library(ChIPseeker)
library(GenomicRanges)

#####################################################
#Este es el archivo del genoma completo del mouse con sus anotaciones y demas.
# Crear un objeto TxDb a partir de tu archivo GFF o GTF
txdb <- makeTxDbFromGFF("gcode.vM10.annotate.gff3", format = "gff3")

######Arreglando el archivo de peaks para hacer el analisis #####

# Leer el archivo .bed (ajusta la ruta de tu archivo)
peaks_data <- read.table("Nd2_summits_200bp.bed", header = FALSE)
write.csv(peaks_data, "Nd2_peaks2")

#############################################################################
total_peaks <- nrow(peaks_data)  # peaks_data es tu archivo de picos
cat("Total Peaks Detected:", total_peaks, "\n")

library(dplyr)

# Contar los picos por tipo de región (annotation viene de peakAnno)
peak_distribution <- as.data.frame(peakAnno@annoStat)

# Mostrar la distribución
print(peak_distribution)


# Asegurarte de que el archivo tenga las columnas correctas
# Columnas: chr, start, end, name, signal
colnames(peaks_data) <- c("chr", "start", "end", "name", "signal")

# Crear el objeto GRanges con las coordenadas y metadata
nd2 <- GRanges(seqnames = peaks_data$chr,
              ranges = IRanges(start = peaks_data$start, end = peaks_data$end),
              strand = "*",  # Aquí puedes ajustar si tienes información sobre la cadena
              name = peaks_data$name,
              signal = peaks_data$signal)

# Ver el GRanges con la metadata añadida
head(nd2)
write.csv(nd2, "Nd2_peaks2")
####################################################
#ChIP peaks coverage plot
covplot(nd2, weightCol="signal")
#specific over chrs
covplot(nd2, weightCol="signal", chrs=c("chrx", "chry"), xlim=c(4.5e7, 5e7))

####################################################
#Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#Generate the Tag Matrix
# Generate the tagMatrix using your peak file (GRanges object)
tagMatrix <- getTagMatrix(nd2, windows = promoter)

#Heatmap of ChIP binding to TSS regions#
####WARNING LO QUE VIENE OCUPA MUCHA MEMORIA OJO.######
# Plot the heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix)

####Los siguiente es para particionar el analisis anterior###
# Dividir tagMatrix en dos partes, por ejemplo
tagMatrix_part1 <- tagMatrix[1:(nrow(tagMatrix)/2), ]
tagMatrix_part2 <- tagMatrix[((nrow(tagMatrix)/2) + 1):nrow(tagMatrix), ]

# Generar un mapa de calor para cada parte
tagHeatmap(tagMatrix_part1)
tagHeatmap(tagMatrix_part2)

data("tagMatrixList")
tagMatrix <- tagMatrixList[[4]]

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

head(promoter)
##################################################
# Peak Annotation
peakAnno <- annotatePeak(nd2, tssRegion=c(-3000, 3000), TxDb=txdb)
head(peakAnno)

library(ggplot2)

plotAnnoPie(peakAnno)
library(ChIPseeker)

# Crear el gráfico de pastel con porcentajes y título
plotAnnoPie(peakAnno, title = "Distribution of genomic annotations") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Ajustar el título
  )



plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)

upsetplot(peakAnno, vennpie=TRUE)

plotDistToTSS(peakAnno,
              title="Distribution of TF NeuroD2 to TSS")

####Im trying to do others things 
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Or your organism-specific TxDb
library(org.Mm.eg.db) # For Mus musculus gene annotations
library(ReactomePA)

# Load TxDb for Mus musculus
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Link peaks to genes
gene_ids <- seq2gene(combined_data_nd2, tssRegion=c(-3000, 3000), TxDb=txdb)

library(GenomicRanges)
colnames(combined_data_nd2)
library(GenomicRanges)

# Crear el objeto GRanges a partir de combined_data_nd2
combined_data_nd2_gr <- GRanges(
  seqnames = combined_data_nd2$chr,  # Columna con los cromosomas
  ranges = IRanges(
    start = combined_data_nd2$start...13,  # Columna con el inicio
    end = combined_data_nd2$end           # Columna con el fin
  ),
  strand = combined_data_nd2$strand,       # Usa la columna 'strand'
  intensity = combined_data_nd2$intensity # Metadatos opcionales
)

# Verifica el objeto
combined_data_nd2_gr

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

gene_ids <- seq2gene(
  combined_data_nd2_gr, 
  tssRegion = c(-3000, 3000), 
  TxDb = txdb, 
  flankDistance = 5000
)

# Verifica los IDs de genes
head(gene_ids)

library(ReactomePA)

pathway1 <- enrichPathway(gene_ids, organism = "mouse")
head(pathway1)

library(clusterProfiler)

# Bar plot for top 10 pathways
barplot(pathway1, showCategory = 10, 
        title = "Enriched pathways of NeuroD2",
        font.size = 12)  # Ajusta el tamaño de letra si es necesario

# Dot plot for top 10 pathways
dotplot(pathway1, showCategory = 10, 
        title = "Top 10 Enriched Pathways (Dot Plot)",
        font.size = 12)

################################################
library(clusterProfiler)
library(org.Mm.eg.db)

# Perform GO enrichment for biological processes
go_bp <- enrichGO(gene = gene_ids, 
                  OrgDb = org.Mm.eg.db, 
                  keyType = "ENTREZID", 
                  ont = "BP",  # "BP" for Biological Process
                  pvalueCutoff = 0.05)

# View results
head(go_bp)
barplot(go_bp, showCategory = 15,title = "Enriched biological processes of NeuroD2",
        font.size = 12)  # Visualize top 10 enriched biological processes

library(clusterProfiler)
library(org.Mm.eg.db)

# Perform GO enrichment
go_bp <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",   # Biological Processes
  pvalueCutoff = 0.05
)

# View results
head(go_bp)

library(clusterProfiler)

# Recalculate pairwise similarity for GO terms
go_bp <- pairwise_termsim(go_bp)

# Check if the similarity matrix is present
go_bp@termsim  # Should now have values


library(enrichplot)

# Create an enrichment map for GO terms
emapplot(go_bp, showCategory = 15)  # Show the top 30 categories


###############################################
library(enrichplot)

emapplot(pairwise_termsim(pathway1))


cnetplot(pathway1, 
         categorySize = "geneName", 
         foldChange = NULL,
         title = "Gene-Pathway Relationship Network")

cnetplot(pathway1, 
         categorySize = "geneNum", 
         foldChange = NULL,
         showCategory = 10, # Muestra solo las 10 pathways más significativas
         title = "Gene-Pathway Relationship Network",
         node_label = "all") # Esto muestra tanto los nombres de los genes como de las categorías


library(org.Mm.eg.db)  # Base de datos de anotaciones para ratón


library(AnnotationDbi)

# Mapea los Entrez IDs a nombres de genes
gene_names <- mapIds(org.Mm.eg.db, 
                     keys = names(pathway1@gene), # Extrae los genes de pathway1
                     column = "SYMBOL",           # Obtén los nombres comunes (símbolos)
                     keytype = "ENTREZID",        # Los IDs actuales son Entrez IDs
                     multiVals = "first")         # Si hay varios nombres, toma el primero


keys <- names(pathway1@gene) # Extrae los identificadores (Entrez IDs)
head(keys)                  # Revisa los primeros valores

keys <- keys[!is.na(keys) & keys != ""]

library(AnnotationDbi)
library(org.Mm.eg.db) # Para Mus musculus

gene_names <- mapIds(org.Mm.eg.db, 
                     keys = keys, 
                     column = "SYMBOL", 
                     keytype = "ENTREZID", 
                     multiVals = "first")


head(gene_names)

keys <- names(pathway1@gene)   # Extrae los Entrez IDs
keys <- keys[!is.na(keys)]     # Elimina valores NA
keys <- keys[keys != ""]       # Elimina cadenas vacías

head(keys)  # Asegúrate de que ahora tienes identificadores válidos


library(AnnotationDbi)
library(org.Mm.eg.db)

gene_names <- mapIds(
  org.Mm.eg.db, 
  keys = keys, 
  column = "SYMBOL", 
  keytype = "ENTREZID", 
  multiVals = "first"
)

pathway1 <- enrichPathway(gene = unique(as.data.frame(peakAnno)$gene_Id), 
                          organism = "mouse")

head(combined_data_nd2)




## Ordenar los cromosomas por frecuencia
chromosome_distribution <- sort(table(seqnames(combined_data_nd2_gr)), decreasing = TRUE)

# Graficar
barplot(chromosome_distribution,
        main = "Distribution of NeuroD2 Peaks Across Chromosomes",
        xlab = "Chromosome",
        ylab = "Number of Peaks",
        las = 2,  # Rotar etiquetas del eje X
        col = "steelblue",
        cex.names = 0.8)  # Reducir tamaño de etiquetas si hay muchas


# Ajustar el límite máximo del eje Y
max_y <- max(chromosome_distribution) * 1.0  # Aumenta un 20% el límite para espacio

barplot(chromosome_distribution,
        main = "Distribution of NeuroD2 Peaks Across Chromosomes",
        xlab = "Chromosome",
        ylab = "Number of Peaks",
        las = 2,              # Rotar etiquetas del eje X
        col = "steelblue",
        cex.names = 0.8,      # Reducir tamaño de las etiquetas del eje X
        ylim = c(0, max_y),   # Ajustar límite del eje Y
        cex.main = 0.9)       # Ajustar tamaño del título


# Ajustar márgenes
par(mar = c(7, 5, 6, 5))  # Márgenes: abajo, izquierda, arriba, derecha

barplot(chromosome_distribution,
        main = "Distribution of NeuroD2 Peaks Across Chromosomes",
        xlab = "Chromosome",
        ylab = "Number of Peaks",
        las = 2,
        col = "steelblue",
        cex.names = 0.8,
        ylim = c(0, max_y),
        cex.main = 0.9)
























########################################################
# Extraer todas las anotaciones
peak_data <- as.data.frame(peakAnno)
write.csv(peak_data, "Nd2_peak_data")

# Filtrar solo los promotores
promoters <- peak_data[peak_data$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)"), ]

# extraer solo las siguientes columnas 
nd2_promoters <- promoters[, c("seqnames", "start", "end", "name", "signal", "annotation")]

# Verifica que solo tenga las columnas deseadas
head(nd2_promoters)
# Guardar en formato .bed
write.table(nd2_promoters, "nd2_promoters.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Filtrar solo las regiones Distal Intergenic como posibles enhancers
enhancers <- peak_data[peak_data$annotation == "Distal Intergenic", ]

# Asumiendo que tu data frame se llama peak_data
nd2_enhancers <- enhancers[, c("seqnames", "start", "end", "name", "signal", "annotation")]
head(nd2_enhancers)

# Guardar en formato .bed
write.table(nd2_enhancers, "nd2_enhancers.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


########################################################
# Instala rtracklayer si aún no lo tienes
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

# Cargar el paquete
library(rtracklayer)

# Leer el archivo GFF3
m10_genome_gff3 <- import("gcode.vM10.annotate.gff3")

# Filtrar genes
genes <- subset(m10_genome_gff3, type == "gene")
##############################################################################
# Filtrar transcritos
transcripts <- subset(m10_genome_gff3, type == "transcript")
######################################################################
# Seleccionar las columnas necesarias para el formato BED
genes_bed <- genes[, c("seqnames", "start", "end")]

# Guardar en archivo BED
write.table(genes_bed, "genes.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


########################################################
#Functional enrichment analysis
library(ReactomePA)
library(org.Mm.eg.db)
library(clusterProfiler)
#############################################


library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)








##############################################3
############
library(org.Mm.eg.db)
# Listar algunos identificadores válidos de ENSEMBL en el objeto org.Mm.eg.db
head(keys(org.Mm.eg.db, keytype = "ENSEMBL"))

# Usar seq2gene para obtener los genes cercanos a los picos (con el TxDb para el ratón)
gene <- seq2gene(nd2, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)

# Asegurarte de que los genes estén en un formato adecuado (por ejemplo, ENTREZID o SYMBOL)
library(org.Mm.eg.db)
head(keys(org.Mm.eg.db, keytype = "ENSEMBL"))

# Confirma que gene contiene identificadores válidos en formato de texto
gene <- as.character(gene)
head(gene)
# Extraer solo el ID de gen ENSEMBL (ENSMUSG)
gene <- sub(".*/(ENSMUSG[0-9]+).*", "\\1", gene)

# Verificar los resultados
head(gene)

# Aquí hacemos una conversión en caso de que los genes sean ENSEMBL
gene_ids <- bitr(gene, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Mm.eg.db)

# Realizar el enriquecimiento de rutas
pathway2 <- enrichPathway(gene_ids$ENTREZID)

# Ver los resultados
gene <- seq2gene(nd2, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)

head(pathway2, 2)






###########################################################
#convertir csv to bed 
# Paso 1: Leer el archivo CSV
data <- read.csv("input.csv")

# Paso 2: Seleccionar y reordenar columnas necesarias para el formato BED
# Asegúrate de que las columnas estén en el orden correcto para BED
bed_data <- data[, c("chr", "start", "end", "name", "score")]

# Paso 3: Guardar en formato BED
write.table(bed_data, "output.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

