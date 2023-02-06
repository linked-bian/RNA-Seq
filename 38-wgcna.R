## -----------------------------------------------------------
library(WGCNA)
library(tidyverse)
gene_exp <- read.csv(file = 'data/wgcna/LiverFemale.csv',
                     row.names = 1)

sample_info <- read.csv(file = 'data/wgcna/ClinicalTraits.csv', 
                        row.names = 1)

gene_info <- read.csv(file = 'data/wgcna/GeneAnnotation.csv')


## -----------------------------------------------------------
{
  # 转置
  datExpr0 <- t(gene_exp)
  
  # 缺失数据及无波动数据过滤
  gsg <- goodSamplesGenes(
    datExpr0, 
    minFraction = 1/2 #基因的缺失数据比例阈值
    )
  datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  
  # 通过聚类，查看是否有明显异常样本, 如果有需要剔除掉
  plot(hclust(dist(datExpr)), 
       cex.lab = 1.5,
       cex.axis = 1.5, 
       cex.main = 2)
  
  
  # 如果剩余基因仍然太多(如 5 到10万条)，
  ##而电脑配置内存不够，可以进一步过滤掉
  ##波动小的基因。不要只用差异基因做。
  # library(genefilter)
  # var_gene_exp <- varFilter(
  #   as.matrix(t(datExpr)),
  #   var.func = IQR,
  #   var.cutoff = 0.5, # 实际项目中设置 0.5 以下
  #   filterByQuantile = TRUE)
  # 
  # datExpr <- t(var_gene_exp)
  datExpr[1:4,1:4]
}


## -----------------------------------------------------------
{
  datTraits <- sample_info
  datTraits[1:4, 1:4]
}


## ----message=FALSE------------------------------------------
{
  library(WGCNA)
  # 多线程
  enableWGCNAThreads(nThreads = 6)
  ## disableWGCNAThreads()
    
  # 通过对 power 的多次迭代，确定最佳 power
  sft <- pickSoftThreshold(
    datExpr, 
    powerVector = 1:20, # 尝试 1 到 20
    networkType = "unsigned" 
    )
}


## -----------------------------------------------------------
{
  # 画图
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ggthemes)
  
  fig_power1 <- ggplot(data = sft$fitIndices,
         aes(x = Power,
             y = SFT.R.sq)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    geom_hline(aes(yintercept = 0.85), color = 'red') +
    labs(title = 'Scale independence',
         x = 'Soft Threshold (power)',
         y = 'Scale Free Topology Model Fit,signed R^2') +
    theme_few() +
    theme(plot.title = element_text(hjust = 0.5))
     
  fig_power2 <- ggplot(data = sft$fitIndices,
         aes(x = Power,
             y = mean.k.)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    labs(title = 'Mean connectivity',
         x = 'Soft Threshold (power)',
         y = 'Mean Connectivity') +
    theme_few()+
    theme(plot.title = element_text(hjust = 0.5))
    
  plot_grid(fig_power1, fig_power2)
}


## ----message=FALSE------------------------------------------
{
  net <- blockwiseModules(
    # 0.输入数据
    datExpr, 
    
    # 1. 计算相关系数
    corType = "pearson", # 相关系数算法，pearson|bicor
      
    # 2. 计算邻接矩阵
    power = 6, # 前面得到的 soft power
    networkType = "unsigned", # unsigned | signed | signed hybrid
      
    # 3. 计算 TOM 矩阵
    TOMType = "unsigned", # none | unsigned | signed
    saveTOMs = TRUE,
    saveTOMFileBase = "blockwiseTOM",
  
    # 4. 聚类并划分模块
    deepSplit = 2, # 0|1|2|3|4, 值越大得到的模块就越多越小
    minModuleSize = 30,
    
    # 5. 合并相似模块
    ## 5.1 计算模块特征向量（module eigengenes， MEs），即 PC1
    ## 5.2 计算 MEs 与 datTrait 之间的相关性
    ## 5.3 对距离小于 mergeCutHeight （1-cor）的模块进行合并
    mergeCutHeight = 0.25, 
  
    # 其他参数
    numericLabels = FALSE, # 以数字命名模块
    nThreads = 0, # 0 则使用所有可用线程
    maxBlockSize = 100000 # 需要大于基因的数目
    )
  # 查看每个模块包含基因数目
  table(net$colors) 
}


## -----------------------------------------------------------
{
  library(tidyverse)
  wgcna_result <- data.frame(gene_id = names(net$colors),
                 module = net$colors) %>%
    left_join(gene_info, by = c('gene_id' = 'gene_id')) 
  head(wgcna_result)
}


## -----------------------------------------------------------
{
  # 同时绘制聚类图和模块颜色
  plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = net$colors,
    groupLabels = "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
}


## ----message=FALSE------------------------------------------
{
  # 计算相关性
  moduleTraitCor <- cor(
    net$MEs,
    datTraits,
    use = "p",
    method = 'spearman' # 注意相关系数计算方式
    )
  
  # 计算 Pvalue
  moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor, 
    nrow(datExpr))
}


## -----------------------------------------------------------
{
  # 相关性 heatmap
  sizeGrWindow(10,6)
  
  # 连接相关性和 pvalue
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  
  
  # heatmap 画图
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(net$MEs),
                 ySymbols = names(net$MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}


## -----------------------------------------------------------
{
  MET <- orderMEs(cbind(net$MEs, dplyr::select(datTraits, weight_g)))
  
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = TRUE, 
    plotHeatmaps = FALSE,
    colorLabels = TRUE,
    marHeatmap = c(3,4,2,2))
}


## -----------------------------------------------------------
{
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = FALSE, 
    plotHeatmaps = TRUE,
    colorLabels = TRUE,
    marHeatmap = c(8,8,2,2))
}


## -----------------------------------------------------------
{
  # 需要导出的模块
  my_modules <- c('blue') 
  
  # 提取该模块的表达矩阵
  m_wgcna_result <- filter(wgcna_result, module %in% my_modules)
  m_datExpr <- datExpr[, m_wgcna_result$gene_id]
  
  # 计算该模块的 TOM 矩阵  
  m_TOM <- TOMsimilarityFromExpr(
    m_datExpr,
    power = 6,
    networkType = "unsigned",
    TOMType = "unsigned")
    
  dimnames(m_TOM) <- list(colnames(m_datExpr), colnames(m_datExpr))
    
  # 导出 Cytoscape 输入文件
  cyt <- exportNetworkToCytoscape(
    m_TOM,
    edgeFile = "data/wgcna/CytoscapeInput-network.txt",
    weighted = TRUE,
    threshold = 0.2)
}

