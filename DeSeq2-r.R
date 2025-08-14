# ==============================================================================
# RNA-Seq differential expression tutorial with DESeq2
# GXA Project: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Remove unused columns (gene ID and gene name)
#genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]
head(counts)

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("genotype")
metadata

# Remove spaces in names to avoid DESeq warnings
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
metadata

# Turn genotype into a factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))
metadata$genotype

# ------------------------------------------------------------------------------
# Spot check expression for knockout gene SNAI1
# ------------------------------------------------------------------------------

gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1' ]
gene_counts = counts[gene_id, ]
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data

library(ggplot2)
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()


# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------
# Install DESeq2 if you haven't already
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

# Load the package
library(DESeq2)

# Then run your code
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ genotype)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Ignore genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5)
res

head(iris)
model = lm(Petal.Width ~ Petal.Length, iris)
plot(iris$Petal.Length, iris$Petal.Width)
abline(model)


res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)

genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
res_df[res_df$Gene.Name %in% genes_to_check, ]


# MA plot
plotMA(res)

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')


# Circos plot
BiocManager::install("biomaRt")
library(biomaRt)

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

install.packages("fastmap", lib = "C:/Users/kaurr/AppData/Local/R/win-library/4.3")
install.packages("RSQLite", lib = "C:/Users/kaurr/AppData/Local/R/win-library/4.3")

library(biomaRt)

# Find dataset name in Ensembl
ensembl <- useEnsembl(biomart="genes")
datasets = listDatasets(ensembl)
head(datasets)

dataset_nb = grep("human", datasets$description, ignore.case=TRUE)
dataset_nb


dataset = datasets$dataset[dataset_nb]
dataset

ensembl <- useDataset(dataset=dataset, mart=ensembl)

# Get coordinates of all genes
attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
values <- list(ensembl_gene_id=c())
all.genes <- getBM(attributes=attributes, values=values, mart=ensembl)
head(all.genes)

# Rename column so it matches res_df
head(res_df)
colnames(all.genes)[1] = "Gene.ID"
head(all.genes)

# Merge the DESeq output with the Ensembl gene coordinates
merged_data <- merge(all.genes, res_df, by="Gene.ID")
head(merged_data)

# 1. Make sure the folder exists
dir.create("~/Documents/rna-seq", recursive = TRUE, showWarnings = FALSE)

# 2. Then write the file
write.csv(merged_data, "~/Documents/rna-seq/deseq.csv", row.names = FALSE)

# Add chr prefix to chromosome names
merged_data$chromosome_name <- paste("chr", merged_data$chromosome_name, sep = "")
head(merged_data)
write.csv(merged_data, "~/Documents/rna-seq/deseq.csv", row.names = FALSE)

# Same for subset
merged_data_subset = merged_data[merged_data$Gene.Name %in% genes_to_check, ]
head(merged_data_subset)
write.csv(merged_data_subset, "~/Documents/rna-seq/deseq_subset.csv", row.names = FALSE)
