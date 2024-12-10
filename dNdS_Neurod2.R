# Set working directory
setwd("/home/maider/Documents/PGB/PGB-Project/dNdS_NEUROD2")

# Load necessary libraries
library(biomaRt)
library(dplyr)

# Read dN/dS table and transcription factors dataset
human_dNdS <- read.table("HumanDnDsW.txt", header = TRUE, sep = "\t")
tf_data <- read.table("HumanTFs_DBD.txt", header = TRUE, sep = "\t")
tf_data <- tf_data[, 2:4]

# Initialize Ensembl Mart
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://oct2014.archive.ensembl.org")
ddatasets= listDatasets(ensembl)
attributes <- listAttributes(ensembl)
list <-listFilters(ensembl)

# Define Ensembl IDs of the transcription factors
tf_ids <- tf_data$Ensembl.ID  

# Step 1: Identify unique new genes that are not already in the GeneId column 
# and create a new dataframe with the unique new genes
new_rows <- data.frame(
  ensembl_gene_id = setdiff(tf_ids, human_dNdS$GeneId), 
  Homo_sapiens.dN = NA, 
  Homo_sapiens.dS = NA,
  Homo_sapiens.dw = NA, 
  Homo_sapiens.w = NA
)  

# Rename column in human_dNdS to match new_rows
colnames(human_dNdS)[1] <- "ensembl_gene_id"
# Step 3: Append the new rows to the existing dataframe
human_dNdS <- rbind(human_dNdS, new_rows)
human_dNdS <- human_dNdS %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), ~ {
    # Remove NA values and collapse
    non_na_values <- na.omit(.)
    if (length(non_na_values) == 0) {
      return(NA)  # If all values are NA, return NA
    } else {
      return(non_na_values[1])  # Choose the first non-NA value
    }
  }))

human_dNdS = human_dNdS[,1]

# Fetch basic gene information
name_info <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = human_dNdS$ensembl_gene_id,
  mart = ensembl
)

human_dNdS = full_join(name_info, human_dNdS, by = "ensembl_gene_id")

# Define species and homolog attributes
species_list <- list(
  Mouse = c("mmusculus_homolog_ds", "mmusculus_homolog_dn"),
  Macaque = c("mmulatta_homolog_ds", "mmulatta_homolog_dn")
  
)

# Create an empty data frame to store the results
result_df <- data.frame(ensembl_gene_id = human_dNdS$ensembl_gene_id)

# Iterate over each species in the species_list
for (species in names(species_list)) {
  species_info <- getBM(
    attributes = c('ensembl_gene_id', species_list[[species]]),
    filters = 'ensembl_gene_id',
    values = human_dNdS$ensembl_gene_id,
    mart = ensembl
  )
  
  # Merge the species_info into the result_df by ensembl_gene_id
  result_df <- merge(result_df, species_info, by = 'ensembl_gene_id')
} 

# Step 1: Replace 0 with NA in result_df
result_df[result_df == 0] <- NA

result_df <- result_df %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), ~ {
    # Remove NA values and collapse
    non_na_values <- na.omit(.)
    if (length(non_na_values) == 0) {
      return(NA)  # If all values are NA, return NA
    } else {
      return(non_na_values[1])  # Choose the first non-NA value
    }
  }))

dNdS = full_join(human_dNdS, result_df, by = "ensembl_gene_id")

write.csv(dNdS,file="dNdS")
dNdS = read.csv(file = "dNdS.csv")
dNdS= dNdS[,2:9]

library("GO.db")
library(annotate)

GOTERM[["GO:0001813"]]
xx <- as.list(GOTERM)
names(xx)

Term(xx[[1]])
terms <- lapply(xx, Term)

immune_terms <- grep("complement", terms, ignore.case = TRUE)
immune_terms <- unlist(terms[immune_terms])
names(immune_terms)


splicing_terms <- grep("splicing", terms, ignore.case = TRUE)
splicing_terms <- unlist(terms[splicing_terms])
names(splicing_terms)

neural_terms <- grep("neuro", terms, ignore.case = TRUE)
neural_terms <- unlist(terms[neural_terms])
names(neural_terms)

tf_terms <- grep("transcription factor", terms, ignore.case = TRUE)
tf_terms <- unlist(terms[tf_terms])
names(tf_terms)



immune_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                     filters = 'go_id', values = names(immune_terms), mart = ensembl)

splicing_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                       filters = 'go_id', values = names(splicing_terms), mart = ensembl)               

neural_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                           filters = 'go_id', values = names(neural_terms), mart = ensembl)               

tf_genes <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                     filters = 'go_id', values = names(tf_terms), mart = ensembl)               


#### Subset dnds table using these two sets of genes

#Musmusculus
df_immuno <- filter(dNdS, ensembl_gene_id %in% unique(immune_gene$ensembl_gene_id))[,c(1,2,3,4,5)]
df_immuno$Geneset <- "Immunological"

df_splicing <- filter(dNdS, ensembl_gene_id %in% unique(splicing_gene$ensembl_gene_id))[,c(1,2,3,4,5)]
df_splicing$Geneset <- "Splicing"

df_tf<- filter(dNdS, ensembl_gene_id %in% unique(tf_genes$ensembl_gene_id))[,c(1,2,3,4,5)]
df_tf$Geneset <- "Transcription Factor"

df <- rbind(df_immuno, df_splicing, df_tf)

library(ggplot2)
library(ggstatsplot)
library(hrbrthemes)
library(viridis)

df %>% ggplot( aes(x=Geneset, y= mmusculus_homolog_w , fill=Geneset)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(data = df %>% filter(ensembl_gene_id == "ENSG00000171532"), aes(x = Geneset, y = mmusculus_homolog_w),
             color = "red", size = 3) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11, hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5, size = 12),  # Move x-axis label down
    axis.title.y = element_text(hjust = 0.5, size = 12)  
    
  ) + ylim(c(0,1)) +
  labs(x = "Geneset", y = "dN/dS", title = "Homos sapiens vs Mus musculus")


#Mulattata
df_immuno <- filter(dNdS, ensembl_gene_id %in% unique(immune_gene$ensembl_gene_id))[,c(1,2,6,7,8)]
df_immuno$Geneset <- "Immunological"

df_splicing <- filter(dNdS, ensembl_gene_id %in% unique(splicing_gene$ensembl_gene_id))[,c(1,2,6,7,8)]
df_splicing$Geneset <- "Splicing"


df_tf<- filter(dNdS, ensembl_gene_id %in% unique(tf_genes$ensembl_gene_id))[,c(1,2,6,7,8)]
df_tf$Geneset <- "Transcription Factor"

df <- rbind(df_immuno, df_splicing, df_tf)

library(ggplot2)
library(ggstatsplot)
library(hrbrthemes)
library(viridis)

df %>% ggplot( aes(x=Geneset, y= mmulatta_homolog_w , fill=Geneset)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(data = df %>% filter(ensembl_gene_id == "ENSG00000171532"), aes(x = Geneset, y = mmulatta_homolog_w),
             color = "red", size = 3) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11, hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5, size = 12),  # Move x-axis label down
    axis.title.y = element_text(hjust = 0.5, size = 12)  
   
  ) + ylim(c(0,1)) +
 labs(x = "Geneset", y = "dN/dS", title = "Homos sapiens vs Macaca mulatta")


library(dplyr)

# Calculate the percentile position of the target gene value within each Geneset
df <- df %>%
  group_by(Geneset) %>%
  mutate(
    # Calculate ECDF for each Geneset group
    percentile = ecdf(mmusculus_homolog_w)(mmusculus_homolog_w) * 100
  ) %>%
  ungroup()

# Filter the data to show only the target gene's percentile position
target_gene_position <- df %>%
  filter(hgnc_symbol == "NEUROD2") %>%
  dplyr::select(Geneset, mmusculus_homolog_w, percentile)

# Print the percentile position of the target gene
print(target_gene_position)

