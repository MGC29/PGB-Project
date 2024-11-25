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


