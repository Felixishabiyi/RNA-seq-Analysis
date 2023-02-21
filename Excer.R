#excerise

#extracting the count data
counts_exe <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")

#extracting the metadata.
coldata_exe <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")

#reading the raw count data
counts_e <- as.matrix(read.table(counts_exe, header = T, sep = '\t'))

#obtainig the dimension of the count data
dim(counts_e)

#computing the summary of the count data
summary(counts_e)

#extracting the genelength
genewidth <- as.vector(subset(counts_e, select = c(width)))

#read per kilobase
rpkb <- apply(subset(counts_e, select = c(-width)),2,
              function(x) x/(genewidth/1000))

#tpm normalization
tpm <- apply(rpkb, 2, function(x) x/sum(as.numeric(x)) * 10^6)

#sum of columns
colSums(tpm)

#to compute variable genes in the count data
varia <- apply(tpm,1,var)

#top 500 most variable genes
top_500 <- names(varia[order(varia, decreasing = T)][1:500])
top_500

#top 100 most variable genes
top_100 <- names(varia[order(varia, decreasing = T)][1:100])
top_100

#reading the raw meta data
metadata_e <- read.table(coldata_exe, header = T, sep = '\t',
                         stringsAsFactors = T)

#heatmap plot for top 500 variable genes
library(pheatmap)
pheatmap(tpm[top_500,], scale = 'row',show_rownames = F, 
         annotation_col = metadata_e)

#heatmap for top 100 variable genes
pheatmap(tpm[top_100,], scale = 'row',show_rownames = F, 
         annotation_col = metadata_e)

#drawing correlation plots
#loading the corr library
library(corrplot)

#converting the tpm to a correlation matrix
correlationMatrix <-cor(tpm)

corrplot(correlationMatrix, order = 'hclust',
         addrect = 2,addCoef.col = 'white',number.cex = 0.7)

#To derive quick and accurate clusters

#calculating the mean of the rows to extract the top 100
#hughly expressed genes
t100 <- order(rowMeans(tpm),decreasing = T)[1:100]
high_100 <- tpm[t100,]

library(pheatmap)
pheatmap(high_100, scale = 'row', show_rownames = F,
         annotation_col = metadata_e)

library(stats)
library(ggplot2)

#transposing the matrix
trans <- t(high_100)
#converting to log tables
trans <- log2(trans + 1)
#computing the PCA
pcares <- prcomp(trans)

library(ggfortify)
autoplot(pcares,data = metadata_e, colour = 'group')

#Differential expression analysis

#removing the 'width' column from the count matrix
counts_data <- as.matrix(subset(counts_e, select = c(-width)))
counts_data

#the metadata
metadata_e

#designing the formula
design <- "~ group"

library(DESeq2)
library(stats)

ddss <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metadata_e,
                              design = as.formula(design))
print(ddss)

#to obtain the row names of the dds dataset
rownames(ddss)

#to obtain the colnames of the metadata of the dds dataset
colData(ddss)

#to obtain the count matrix of the dds dataset
counts(ddss)

ddss <- ddss[rowSums(DESeq2::counts(ddss)) > 1,]

#implementing the DESeq() function on the dds dataset
ddss <- DESeq(ddss)

#computing the contrast 

DEresults = results(ddss, contrast = c("group",'CASE','CTRL'))
print(DEresults)

DESeq2::plotMA(object = ddss,ylim = c(-5, 5))

library(EnhancedVolcano)

rownames(DEresults)

#volcano plot for the up and down regulated genes
EnhancedVolcano(DEresults, x = 'log2FoldChange', y = 'pvalue', 
                lab = rownames(DEresults))

#plot for a particular label to show
selectt <- c('CYP2E1','GCKR','AADAC','FCGBP','ZG16','CA4')
EnhancedVolcano(DEresults, x = 'log2FoldChange', y = 'pvalue', 
                lab = rownames(DEresults), selectLab = selectt)

#ploting a dispersion plot
plotDispEsts(object = ddss,legend = F)

#Enrichment Analysis
library(org.Hs.eg.db)
library(clusterProfiler)

#removing rows with NA columns
DEresults <- na.omit(DEresults)

#ordering the DEG list by the stat column
DEresults <- DEresults[order(DEresults$stat,decreasing = T),]

#creating a gene list with the stat value
gl <- DEresults$stat

#assigning names to the gl list using the DEG list
names(gl) <- rownames(DEresults)

#converting the gl list to a gse format 
gse <- gseGO(gl,ont="BP",keyType = "SYMBOL",
             OrgDb = "org.Hs.eg.db",eps = 1e-300)

#obtain the dataframe format of the gse
as.data.frame(gse)

#obtaining plots for the Enrichment analysis
gseaplot(gse,geneSetID = 18)