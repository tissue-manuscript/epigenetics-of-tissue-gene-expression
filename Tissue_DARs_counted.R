# loading
load("ATAC_counted_SKM_vs_tissue.RData")

library(DiffBind)

jpeg("heatmap_counted.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotHeatmap(data_counted)
dev.off()

# normalize
cat("Normalizing...\n")
data_normalize <- dba.normalize(data_counted)
save.image("02-normalized_data.RData")

# sink(file="02-normalized_data.txt")
# data_normalize
# data_normalize$norm
# sink(file=NULL)

jpeg("pca_normalized_id.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotPCA(data_normalize, label=DBA_ID)
dev.off()

jpeg("pca_normalized_condition.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotPCA(data_normalize, label=DBA_CONDITION)
dev.off()

# contrast
cat("Contrasting data...\n")
#data_contrasted <- dba.contrast(data_normalize, reorderMeta=list(Condition=baseline_condition_name))
data_contrasted <- dba.contrast(data_normalize, reorderMeta=list(Condition="Control"))
# analyze
cat("Analyzing data...\n")
data_analyzed <- dba.analyze(data_contrasted, method=DBA_ALL_METHODS)
save.image("03-analyzed_data.RData")

sink(file="03-analyzed_data.txt")
data_analyzed

cat("\n", "EDGER REPORT:", "\n")
edger_report <- dba.report(data_analyzed, method=DBA_EDGER)
edger_report

cat("\n", "DESEQ2 REPORT:", "\n")
deseq_report <- dba.report(data_analyzed, method=DBA_DESEQ2)
deseq_report

sink(file=NULL)

dba.report(data_analyzed, method=DBA_EDGER, file="output_edger")
dba.report(data_analyzed, method=DBA_DESEQ2, file="output_deseq2")
dba.report(data_analyzed, method=DBA_ALL_METHODS, file="output_all")

jpeg("heatmap_analyzed.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotHeatmap(data_analyzed, contrast=1)
dev.off()

jpeg("volcano_analyzed.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotVolcano(data_analyzed)
dev.off()

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
jpeg("heatmap_analyzed_bindingaffinity.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotHeatmap(data_analyzed, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)
dev.off()

jpeg("venn_analyzed.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotVenn(data_analyzed,contrast=1,bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
dev.off()

jpeg("ma_analyzed.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotMA(data_analyzed)
dev.off()

jpeg("boxplot_analyzed.jpeg", units="px", width=2000, height=2000, res=300)
dba.plotBox(data_analyzed)
dev.off()

# export to BED
library(dplyr)

edger_out <- as.data.frame(edger_report)
deseq_out <- as.data.frame(deseq_report)

edger_up <- edger_out %>% filter(Fold > 0) %>% select(seqnames, start, end)
write.table(edger_up, file="edger_up.bed", sep="\t", quote=F, row.names=F, col.names=F)
system("annotatePeaks.pl edger_up.bed rn6 > edger_up_annotated.txt")

edger_down <- edger_out %>% filter(Fold < 0) %>% select(seqnames, start, end)
write.table(edger_down, file="edger_down.bed", sep="\t", quote=F, row.names=F, col.names=F)
system("annotatePeaks.pl edger_down.bed rn6 > edger_down_annotated.txt")

deseq_up <- deseq_out %>% filter(Fold > 0) %>% select(seqnames, start, end)
write.table(deseq_up, file="deseq_up.bed", sep="\t", quote=F, row.names=F, col.names=F)
system("annotatePeaks.pl deseq_up.bed rn6 > deseq_up_annotated.txt")

deseq_down <- deseq_out %>% filter(Fold < 0) %>% select(seqnames, start, end)
write.table(deseq_down, file="deseq_down.bed", sep="\t", quote=F, row.names=F, col.names=F)
system("annotatePeaks.pl deseq_down.bed rn6 > deseq_down_annotated.txt")

#Annotate output_upregulated

library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
library(clusterProfiler)


workdir = "/directory_name/"
filename = "deseq_up_annotated.txt"

setwd(workdir)
filepath = file.path(workdir, filename)

timest = toString(round(as.numeric(Sys.time())*1000))
new_folder_name = sprintf("%s_results", timest)
dir.create(new_folder_name)

# switch to results folder
setwd(new_folder_name)

peak <- readPeakFile(filepath)

# heatmap
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

jpeg('tagHeatmap.jpg', width=250, height=750)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()

# peak annotation
peakAnno <- annotatePeak(filepath, tssRegion=c(-3000,3000), TxDb=txdb, annoDb = "org.Rn.eg.db")
write.csv(peakAnno, file = "deseq_up_annotated_peakAnno.csv")

jpeg('annoPlot.jpg', width=750, height=500)
plotAnnoPie(peakAnno)
dev.off()

#Annotate output_downregulated

library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
library(clusterProfiler)


workdir = "/directory_name/"
filename = "deseq_down_annotated.txt"

setwd(workdir)
filepath = file.path(workdir, filename)

timest = toString(round(as.numeric(Sys.time())*1000))
new_folder_name = sprintf("%s_results", timest)
dir.create(new_folder_name)

# switch to results folder
setwd(new_folder_name)

peak <- readPeakFile(filepath)

# heatmap
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

jpeg('tagHeatmap.jpg', width=250, height=750)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()

# peak annotation
peakAnno <- annotatePeak(filepath, tssRegion=c(-3000,3000), TxDb=txdb, annoDb = "org.Rn.eg.db")
write.csv(peakAnno, file = "deseq_down_annotated_peakAnno.csv")

jpeg('annoPlot.jpg', width=750, height=500)
plotAnnoPie(peakAnno)
dev.off()