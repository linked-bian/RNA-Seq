## -----------------------------------------------------------
gene_exp <- read.table(file = 'data/pca/gene_exp.txt', 
                       sep = '\t', header = T, row.names = 1)
sample_info <- read.table(file = 'data/pca/sample_info.txt',
                          sep = '\t', header = T, row.names = 1)


## -----------------------------------------------------------
# 安装
## BiocManager::install('PCAtools')

# 加载
library(PCAtools)


## ---- message=FALSE-----------------------------------------
pca <- pca(gene_exp, metadata  = sample_info)


## -----------------------------------------------------------
pca_loadings <- pca$loadings
pca_loadings[1:4, 1:4]


## -----------------------------------------------------------
pca_rotated <- pca$rotated
pca_rotated[1:4, 1:4]


## ---- message=FALSE-----------------------------------------
screeplot(pca)


## ---- message=FALSE-----------------------------------------
biplot(pca, 
       x = 'PC1',                 # x 轴
       y = 'PC2',                 # y 轴
       colby = 'strain',          # 颜色映射
       shape = 'stage',           # 形状映射
       legendPosition = 'right',  # 图例位置
       lab = NULL                    # 样本名称显示
       )


## -----------------------------------------------------------
plotloadings(pca)


## -----------------------------------------------------------
save(pca, pca_rotated, pca_loadings, file = 'data/pca/pca.rdata')

