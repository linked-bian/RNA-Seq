#gene expression mean
gene1_s1 <- 10
gene1_s2 <- 20
gene1_s3 <- 30
gene1_s4 <- 40
gene1_s5 <- 50
gene1_s6 <- 60

gene1_total<-
  gene1_s1+gene1_s2+gene1_s3+gene1_s4+gene1_s5+gene1_s6

gene1 = c(10,20,30,40,50,60)

gene1_total=0

for (i in gene1) {
  gene1_total=gene1_total+i
}
a=1
while (a<6) {
  gene1_total=gene1_total+gene1[a]
}

gene1_total=sum(gene1)
print(gene1_total)

gene_exp <- read.delim("D:/R.library/RNAseq/Rdata", row.names=1)

apply(gene_exp, 1,sum)
sample_cor<-cor(gene_exp)

write.table



















install.packages('BiocManager')
