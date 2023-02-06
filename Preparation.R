#RNAseq相关软件包的安装
chooseCRANmirror()
install.packages("BiocManager")
BiocManager::install(c("devtools",
                       "tidyverse",
                       "gridExtra",
                       "pheatmap",
                       "ggfortify",
                       "survival",
                       "survminer",
                       "ggpubr",
                       "ggplot2",
                       "ggExtra", 
                       "DESeq2",
                       "edgeR",
                       "GEOquery",
                       "clusterProfiler",
                       "limma",
                       "org.Hs.eg.db", 
                       "TCGAbiolinks", 
                       "biomaRt", 
                       "PCAtools",
                       "EnhancedVolcano",
                       "WGCNA", 
                       "YuLab-SMU/createKEGGdb"))



