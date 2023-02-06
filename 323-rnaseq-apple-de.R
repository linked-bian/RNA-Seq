## -----------------------------------------------------------
load(file = 'data/rnaseq-apple/rprepare.rdata')


## -----------------------------------------------------------
de_result <- read.delim("data/rnaseq-apple/genes.counts.matrix.KID_S1_vs_KID_S4.DESeq2.DE_results")


## -----------------------------------------------------------
library(tidyverse)
de_result <- 
  # 添加上下调信息
  mutate(de_result, direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')
    ))) %>%
  # 关联基因信息
  left_join(gene_info, by = c('id' = 'GID')) %>%
  # 关联表达量信息
  left_join(rownames_to_column(gene_exp, var = 'id'), by = 'id') %>%
  # 去除无用的列
  dplyr::select(-c(2:4, 6:7)) %>%
  # 按 log2FoldChange 绝对值降序排列
  arrange(desc(abs(log2FoldChange)))


## -----------------------------------------------------------
group_by(de_result, direction) %>%
  summarise(count = n())


## -----------------------------------------------------------
write.csv(de_result, 
          file = 'output/de_result.csv', 
          row.names = F, quote = F)


## -----------------------------------------------------------
save(de_result, file = 'data/rnaseq-apple/de.rdata')


## ----eval=FALSE---------------------------------------------
## BiocManager::install('EnhancedVolcano')


## -----------------------------------------------------------
library(EnhancedVolcano)
EnhancedVolcano(de_result,
                lab = de_result$id,
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Volcano Plot',
    subtitle = 'KID_S1_vs_KID_S4',
    xlim = c(-10, 10),
    pCutoff = 0.05,
    FCcutoff = 1)


## -----------------------------------------------------------
top_de_exp <- dplyr::slice(de_result, 1:20) %>%
  dplyr::select(id, contains('KID_S1'), contains('KID_S4')) %>%
  column_to_rownames(var = 'id')


## -----------------------------------------------------------
library(pheatmap)
pheatmap(top_de_exp, 
         scale = 'row',
         color = colorRampPalette(c("green", "white", "red"))(200),
         annotation_col = dplyr::select(sample_info, stage),
         annotation_colors = list(
           stage = c(S1 = '#4DBBD5FF', 
                     S4 = '#E64B35FF')),
         cutree_rows = 2)

