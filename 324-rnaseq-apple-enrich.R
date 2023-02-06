## -----------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
load(file = 'data/rnaseq-apple/rprepare.rdata')
load(file = 'data/rnaseq-apple/de.rdata')

# 如果后面的富集分析运行太慢，可以加载我预先分析好的结果
load(file = 'data/rnaseq-apple/enrich.rdata')


## -----------------------------------------------------------
gene <- filter(de_result, 
               abs(log2FoldChange) > 1 & padj < 0.05) %>%
  pull(id)


## -----------------------------------------------------------
geneList <- de_result$log2FoldChange
names(geneList) <- de_result$id
geneList <- sort(geneList, decreasing = T)


## Rscript emcp/create_orgdb_from_emapper.R query_seqs.fa.emapper.annotations


## -----------------------------------------------------------
# 创建一个文件夹
dir.create('R_Library', recursive = T)

# 将包安装在该文件夹下
install.packages('data/rnaseq-apple/org.My.eg.db_1.0.tar.gz', 
                 repos = NULL, #从本地安装
                 lib = 'R_Library') # 安装文件夹

# 加载 OrgDB
library(org.My.eg.db, lib = 'R_Library')


## ----eval=FALSE---------------------------------------------
## de_ego <- enrichGO(gene = gene,
##          OrgDb = org.My.eg.db,
##          keyType = 'GID',
##          ont = 'ALL',
##          qvalueCutoff = 0.05,
##          pvalueCutoff = 0.05)


## -----------------------------------------------------------
de_ego_df <- as.data.frame(de_ego)
head(de_ego_df)


## -----------------------------------------------------------
barplot(de_ego, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")


## -----------------------------------------------------------
dotplot(de_ego, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")


## -----------------------------------------------------------
cnetplot(de_ego, 
         foldChange = geneList, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)


## ----eval=FALSE---------------------------------------------
## emapplot(de_ego, showCategory = 10, pie = 'count')


## -----------------------------------------------------------
emapper <- read_delim('data/rnaseq-apple/query_seqs.fa.emapper.annotations', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X9, 
                Pathway = X10)

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))


## -----------------------------------------------------------
library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()


## ----eval=FALSE---------------------------------------------
## library(clusterProfiler)
## de_ekp <- enricher(gene,
##                 TERM2GENE = pathway2gene,
##                 TERM2NAME = pathway2name,
##                 pvalueCutoff = 0.05,
##                 qvalueCutoff = 0.05)


## -----------------------------------------------------------
de_ekp_df <- as.data.frame(de_ekp)
head(de_ekp_df)


## -----------------------------------------------------------
barplot(de_ekp, showCategory = 10)


## -----------------------------------------------------------
dotplot(de_ekp, showCategory = 10)


## -----------------------------------------------------------
cnetplot(de_ekp, 
         foldChange = geneList, 
         showCategory = 3,
         node_label = "category", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)


## -----------------------------------------------------------
emapplot(de_ekp, showCategory = 10, pie = 'count')


## -----------------------------------------------------------
save(de_ego, de_ekp, file = 'data/rnaseq-apple/enrich.rdata')

