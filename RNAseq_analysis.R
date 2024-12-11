# The data used can be accessed with the following link.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110491
# The GEO accession number is GSE110491. 

#STEP 1: Read the count matrix and process it. 
counts<-read.csv(file="GSE110491_rawCounts.txt",sep="\t",header=T,dec=".", stringsAsFactors = F)
colnames(counts)<-gsub("^output.bam.","",colnames(counts))
colnames(counts)<-gsub(".bam","",colnames(counts))
counts.m <- as.matrix(counts[,7:ncol(counts)])
rownames(counts.m)<-counts$Geneid


#STEP 2: Perform the differential expression analysis. 
library(edgeR)
library(dplyr)
groups = c("KO", "KO", "KO", "WT", "WT", "WT")
dgList <- DGEList(counts=counts.m, group = factor(groups))
cpm.matrix <- cpm(dgList)

# We have set a threshold of 5 counts per million in at least 3 of the samples. 
countCheck <- cpm.matrix > 5
keep <- which(rowSums(countCheck) >= 3)
dgList <- dgList[keep,]
dgList$samples$lib.size <- colSums(dgList$counts)

# Here we perform the differential gene expression analysis and see which features are significant. 
dgList <- calcNormFactors(dgList, method="TMM")
dgList <- estimateDisp(dgList)
et <- exactTest(dgList, pair = c("WT", "KO"))
significant <- decideTests(et)
summary(decideTests(et,))

#STEP 3: Represent significant features with volcano plot. 
library(EnhancedVolcano)

# To properly plot the features we have to calculate the adjusted p-value and 
# add it to our "et" table. We used the Benjamini & Hochberg method as it is 
# the test used by default in the decideTests() function. 
padjust<-as.data.frame(p.adjust(et$table$PValue, method="BH"))
et$table$PValAdjust <- padjust$`p.adjust(et$table$PValue, method = "BH")`
colnames(et$table)[4] <- "PValAdjust"

# Additionally we will need a vector with the names of the significant features. 
keep <- which(significant != 0)
significant <- significant[keep,]
sig_genes <- rownames(significant@.Data)

# Here we plot the results 
EnhancedVolcano(et$table,
                lab = rownames(et$table),
                x = 'logFC',
                y = 'PValAdjust',
                ylab = expression(-~Log[10] ~ P ~ adjusted),
                axisLabSize = 14,
                selectLab = sig_genes,
                pCutoff = 0.05,
                FCcutoff = abs(0.5),
                pointSize = 3.0, 
                title = NULL,
                subtitle = NULL, 
                legendLabels = c("Not significant", expression(Log[2] ~ FC), 
                                 "p-value", expression(Log[2] ~ FC ~ and ~ p-value)),
                legendLabSize = 12,
                legendPosition = "bottom",
                )

# Because there is an outlier, we can modify the axis to eliminate it and see the other 
# features better.
EnhancedVolcano(et$table,
                lab = rownames(et$table),
                x = 'logFC',
                y = 'PValAdjust',
                ylab = expression(-~Log[10] ~ P ~ adjusted),
                selectLab = sig_genes,
                pCutoff = 0.05,
                FCcutoff = abs(0.5),
                pointSize = 3.0, 
                ylim = c(0, 21),
                xlim = c(-3,3),
                title = NULL,
                subtitle = NULL,
                legendLabels = c("Not significant", expression(Log[2] ~ FC), 
                                 "p-value", expression(Log[2] ~ FC ~ and ~ p-value)),
                legendLabSize = 12,
                legendPosition = "bottom")

# STEP 4: Represent significant features with heatmap
library(pheatmap)
library(RColorBrewer)

# First we want to obtain the values of the significant feature.
hits <- et$table[sig_genes, ]

# Then we want to order the hits from the more down expressed to the most over expressed.
order_hits <- order(hits$logFC)
final_order <- c()
for (i in 1:length(order_hits)){
  final_order[i] <- row.names(hits[order_hits[i],])
}

# Then we will read the annotations we have compiled for each significant feature.
# The file contains the type of feature, function and relation to neural function. 
# The table is available in the annexes of the written project. 
annotations <- read.csv("significant_genes_anot.csv", sep=",",
                        header=T,dec=".", stringsAsFactors = F)
row.names(annotations) <- annotations$hgnc_symbol
colnames(annotations)[4] <- "Relation to neural function"

# Here we declare the colours we want our annotations to have.  
colours <- list(
  "Relation to neural function" = c(
    "Synapse" = "#d7191c", "Neurogenesis" = "#fdae61",
    "ND" = "#ffffbf", "Immunological synapses" = "#a6d96a",
    "Action potential" = "#1a9641"
  ), 
  "Type" = c(
    "Protein coding" = "#998ec3", "ND" = "#ffffbf", "lncRNA" =  "#f1a340"
  )
)

heatcol<-colorRampPalette(c("blue", "white","red"), space = "rgb")

# Here we plot the heatmap 
heatmap <- pheatmap(as.matrix(cpm.matrix[final_order, ]),
          main = "Significantly expressed genes",
          col = heatcol(256),
          cluster_rows = F, 
          clustering_method = "ward.D2",
          clustering_distance_cols = "euclidean",
          scale = "row", colCol = heatcol, 
          cexRow = 0.5, cexCol = 0.8, show_rownames = T,
          treeheight_row = 0,
          annotation_row = annotations[, c(2,4)], 
          annotation_colors = colours
)

# With this function we can save the plot as a pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap, "Heatmap_significant_genes")

# STEP 5: GSEA
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(stringr)
library("org.Mm.eg.db", character.only = TRUE)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

# First we need the ENSEMBL IDs of our features 
genes <- getBM(attributes=c('ensembl_gene_id', "mgi_symbol"), 
               filters ='mgi_symbol', values = row.names(et$table), 
               mart = ensembl)
# However for some reason some of the IDs are duplicated and we have to remove them.
genes <- genes[!duplicated(genes$mgi_symbol),]

# Then we want to get the logFC of our features and order them from higher to lower.
logFC <- et$table[row.names(et$table) %in% genes$mgi_symbol, 1]
names(logFC) <- genes$ensembl_gene_id
logFC = sort(logFC, decreasing = TRUE)

# Here we perform the gene set expression analysis 
gse <- gseGO(geneList=logFC, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")

# Because we known NEUROD2 is involved in synapse formation, we want to obtain
# those pathways that are involved in some way in synapse. 
library(forcats)
library("viridis")
synapsis <- grep("*synap*", gse@result$Description)
synapsis <- as.data.frame(gse@result[synapsis,])

# Because we want to plot the pathways according to their enrichment score, we want
# to create a new column that will hold their "order".
order <- order(synapsis$enrichmentScore, decreasing = FALSE)
for (i in 1:length(order)){
  synapsis[order[i], 13] <- i
}
colnames(synapsis)[13] <- "Order"
synapsis$Description = str_to_title(synapsis$Description)

# Here we plot the pathways 
ggplot(data=synapsis, aes(x = reorder(Description, -Order), y = enrichmentScore, fill = p.adjust)) +
  geom_bar(stat="identity") +
  scale_y_continuous(name = "Enrichment Score") +
  scale_x_discrete(name = "") +
  scale_fill_gradient(low="red", high="blue", name = "Adjusted p-value") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"),
        axis.text = element_text(color="black"))

