# Load library for DESeq2
library(DESeq2)
# Load library for RColorBrewer
library(RColorBrewer)
# Load library for pheatmap
library(pheatmap)
# Load library for tidyverse
library(tidyverse)
# ggplot2 library
library(ggplot2)
# load genome
library(org.Rn.eg.db)

#load RData

load("RNA_SKM_vs_tissues.RData")
dds

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# normalize
normlzd_dds <- counts(dds, normalized=T)
head(normlzd_dds)

#save noramlized data
write.csv(normlzd_dds, "normalizedData")

# dendrogram
plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$protocol)

# Varaiance Stabilizing transformation
vsd <- vst(dds, blind = T)
# extract the vst matris from the object
vsd_mat <- assay(vsd)
# compute pairwise correlation values
vsd_cor <- cor(vsd_mat)

# -------------------

pheatmap(vsd_cor)
plotPCA(vsd, intgroup = "tissue")

# Calculating mean for each gene
mean_readCounts <- apply(read_Count[,], 1, mean)
# Calculating variance for each gene
var_readCounts <- apply(read_Count[,], 1, var)

df <- data.frame(mean_readCounts, var_readCounts)
ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")



# deseq2 pipeline

design(dds) <- ~id + tissue
dds <- DESeq(dds)

###

res <-results(dds, contrast = c("tissue", "SKM", "Heart"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_Hea_resSort.csv")
write.csv(geneinfo, "SKM_Hea_geneInfo.csv")

#

res <-results(dds, contrast = c("tissue", "SKM", "Liver"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_Liv_resSort.csv")
write.csv(geneinfo, "SKM_Liv_geneInfo.csv")

#

res <-results(dds, contrast = c("tissue", "SKM", "Hippoc"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_Hippoc_resSort.csv")
write.csv(geneinfo, "SKM_Hippoc_geneInfo.csv")

# 

res <-results(dds, contrast = c("tissue", "SKM", "Lung"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_Lun_resSort.csv")
write.csv(geneinfo, "SKM_Lun_geneInfo.csv")

#

res <-results(dds, contrast = c("tissue", "SKM", "Kidney"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_Kid_resSort.csv")
write.csv(geneinfo, "SKM_Kid_geneInfo.csv")

#

res <-results(dds, contrast = c("tissue", "SKM", "BAT"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_BAT_resSort.csv")
write.csv(geneinfo, "SKM_BAT_geneInfo.csv")

#

res <-results(dds, contrast = c("tissue", "SKM", "WAT"))
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_WAT_resSort.csv")
write.csv(geneinfo, "SKM_WAT_geneInfo.csv")