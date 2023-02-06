## -----------------------------------------------------------
sample_info <- read.table("data/rnaseq-apple/sample_info.txt", 
                       header = T, row.names = 1)
head(sample_info)


## -----------------------------------------------------------
gene_counts = read.table("data/rnaseq-apple/genes.counts.matrix", 
                         header=T, row.names=1)


## -----------------------------------------------------------
gene_exp <- read.table("data/rnaseq-apple/genes.TMM.EXPR.matrix", 
                       header = T, row.names=1)


## -----------------------------------------------------------
library(tidyverse)
library(readr)
emapper <- read_delim('data/rnaseq-apple/query_seqs.fa.emapper.annotations', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X6, 
                GO = X7, 
                KO = X9, 
                Pathway = X10, 
                OG = X21, 
                Gene_Name = X22)

gene_info <- dplyr::select(emapper,  GID, Gene_Symbol, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))


## -----------------------------------------------------------
sample_cor <- 
  cor(gene_exp, 
      method = "pearson") # 算法： pearson | kendall | spearman

sample_cor <- round(sample_cor, digits = 1)
sample_cor[1:3,1:3]


## ----eval=FALSE---------------------------------------------
## # 安装
## ## install.packages("pheatmap")
## 
## # 加载
## library(pheatmap)
## 
## # 画图
## pheatmap(sample_cor, cluster_rows = F, cluster_cols = F)


## -----------------------------------------------------------
sample_hc <- hclust(dist(t(gene_exp)))
plot(sample_hc)


## -----------------------------------------------------------
save(gene_counts, gene_exp, 
     sample_info, gene_info, 
     file = 'data/rnaseq-apple/rprepare.rdata')

