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

#2. Sample correlation analysis
##Parson correlation coefficient
sample_cor <- round(cor(gene_exp), digits = 2)
library(pheatmap)
pheatmap(sample_cor)

#3. Sample clustering analysis
##calculate distance matrix
###Transverse firstly
sample_dist <- dist(t(gene_exp))
##Clusting
sample_hc <- hclust(sample_dist)
plot(sample_hc)

