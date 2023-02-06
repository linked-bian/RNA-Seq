#1. Input Data: 
##~/Mithrandir.training/R.data/InputData/genes.counts.matrix
gene_count <- read.table('~/Mithrandir.training/R.data/InputData/genes.counts.matrix', 
                         header = TRUE, row.names = 1)

##~/Mithrandir.training/R.data/InputData/genes.TMM.EXPR.matrix
gene_exp <- read.table('~/Mithrandir.training/R.data/InputData/genes.TMM.EXPR.matrix', 
                       header = TRUE, row.names = 1)

##~/Mithrandir.training/R.data/InputData/query_seqs.fa.emapper.annotations
library(readr)
library(tidyverse)
gene_info <- read_delim("InputData/query_seqs.fa.emapper.annotations", 
                        delim = "\t", 
                        escape_double = FALSE, 
                        col_names = FALSE, 
                        comment = "#", 
                        trim_ws = TRUE)%>%
  select(Gene_Id = X1, 
         Gene_Sumbol = X6,
         GO = X7,
         Ko = X9,
         Pathway = X10,
         COG = X21,
         Gene_Name = X22)

##~/Mithrandir.training/R.data/InputData/sample_info.txt
sample_info <- read.delim("~/Mithrandir.training/R.data/InputData/sample_info.txt", row.names=1)

#2. PCA calculate:
library(PCAtools)
p <- pca(gene_exp, metadata = sample_info, removeVar = 0.1)
p

##View Gene weights in PCA
pca_loadings <- p$loadings
view(pca_loadings)

##View Sample weight in PCA
pca_rotated <- p$rotated
view(pca_rotated)

##Plot PCA Explain Variation
screeplot(p)

##Plot PCA rotated
biplot(p, 
       x = 'PC1',
       y = 'PC2',
       colby = 'strain',
       colkey = c('KID' = 'red', 'BLO' = 'yellow'),
       shape = 'stage',
       legendPosition = 'top',
       showLoadings = TRUE,
       lab = NULL)

##Plot gene in PCA
plotloadings(p)

#3. Batch effect
boxplot(log10(gene_exp))
###library(sva)：remove batch effect if needed, exp: samples from different Cor.
