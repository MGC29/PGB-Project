library(VennDiagram)
library(RColorBrewer)

# Here we read the three files that contain respectively:
# All the differentially expressed genes of the RNAseq analysis.
# Only the down regulated genes from the RNAseq analysis. 
# The genes from the peack of the Chip seq analysis.
rnaseq_genes <- read.csv(file="significant_genes.csv",sep="\t",header=T,dec=".", stringsAsFactors = F)
rnaseq_down <-  read.csv(file="significant_down.csv",sep="\t",header=T,dec=".", stringsAsFactors = F)
chipseq_genes <- read.csv(file="Motif_genes.csv",sep=",",header=T,dec=".", stringsAsFactors = F)

# Here we declare the colours for the Venn diagram
myCol <- c("#B3E2CD", "#FDCDAC")

# Here we plot the Venn diagram for all the RNA seq genes.
venn.diagram(
  x = list(rnaseq_genes$Genes, chipseq_genes$external_gene_name),
  category.names = c("RNA-seq genes" , "Chip-seq\ngenes" ),
  filename = 'venn_diagramm_all.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  fontface = "bold",
  cat.pos = c(-20, 20),
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  margin = 0.05
)

# Here we plot the Venn diagram only for the down regulated genes. 
venn.diagram(
  x = list(rnaseq_down$Genes, chipseq_genes$external_gene_name),
  category.names = c("RNA-seq down\nregulated genes" , "Chip-seq\ngenes" ),
  filename = 'venn_diagramm_downgenes.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  fontface = "bold",
  cat.pos = c(-20, 20),
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  margin = 0.05
)
