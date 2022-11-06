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

load("RNA_ERGs_8tissues.RData")
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

# design(dds) <- ~id + tissue
design (dds) <- ~0 + tissue + id
dds <- DESeq(dds)
resultsNames(dds)

###
#SKM-GN

res <-results(dds, 
              contrast = list(c("tissueSKM"),
                              c("tissueBAT", "tissueHeart", "tissueHippoc", "tissueKidney", "tissueLiver", "tissueLung", "tissueWAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "SKM_vs_All_resSort.csv")
write.csv(geneinfo, "SKM_vs_All_geneInfo.csv")

###
#Heart

res <-results(dds, 
              contrast = list(c("tissueHeart"),
                              c("tissueBAT", "tissueSKM", "tissueHippoc", "tissueKidney", "tissueLiver", "tissueLung", "tissueWAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "Heart_vs_All_resSort.csv")
write.csv(geneinfo, "Heart_vs_All_geneInfo.csv")

###
#Liver

res <-results(dds, 
              contrast = list(c("tissueLiver"),
                              c("tissueBAT", "tissueSKM", "tissueHippoc", "tissueKidney", "tissueHeart", "tissueLung", "tissueWAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "Liver_vs_All_resSort.csv")
write.csv(geneinfo, "Liver_vs_All_geneInfo.csv")

###
#BAT

res <-results(dds, 
              contrast = list(c("tissueBAT"),
                              c("tissueLiver", "tissueSKM", "tissueHippoc", "tissueKidney", "tissueHeart", "tissueLung", "tissueWAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "BAT_vs_All_resSort.csv")
write.csv(geneinfo, "BAT_vs_All_geneInfo.csv")

###
#WAT

res <-results(dds, 
              contrast = list(c("tissueWAT"),
                              c("tissueLiver", "tissueSKM", "tissueHippoc", "tissueKidney", "tissueHeart", "tissueLung", "tissueBAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "BAT_vs_All_resSort.csv")
write.csv(geneinfo, "BAT_vs_All_geneInfo.csv")

###
#Hippoc

res <-results(dds, 
              contrast = list(c("tissueHippoc"),
                              c("tissueLiver", "tissueSKM", "tissueWAT", "tissueKidney", "tissueHeart", "tissueLung", "tissueBAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "Hippoc_vs_All_resSort.csv")
write.csv(geneinfo, "Hippoc_vs_All_geneInfo.csv")

###
#Kidney

res <-results(dds, 
              contrast = list(c("tissueKidney"),
                              c("tissueLiver", "tissueSKM", "tissueWAT", "tissueHippoc", "tissueHeart", "tissueLung", "tissueBAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "Kidney_vs_All_resSort.csv")
write.csv(geneinfo, "Kidney_vs_All_geneInfo.csv")

###
#Lung

res <-results(dds, 
              contrast = list(c("tissueLung"),
                              c("tissueLiver", "tissueSKM", "tissueWAT", "tissueHippoc", "tissueHeart", "tissueKidney", "tissueBAT")),
              listvalues = c(1, -1/7))
              
summary(res)
resSort <- res[order(res$pvalue),]

geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort),
                   columns=c("ENSEMBL", "SYMBOL","GENENAME"), 
                   keytype="ENSEMBL")

write.csv(resSort, "Lung_vs_All_resSort.csv")
write.csv(geneinfo, "Lung_vs_All_geneInfo.csv")