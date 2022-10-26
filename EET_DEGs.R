library(limma)
library(edgeR)
library(biomaRt)
library(qvalue)
library(genefilter)
library(gridExtra)
library(dplyr)
library(data.table)

setwd("~/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data")

load("November 2021 PASS1B Data Freeze/pass1b-06_v1.0_analysis_transcriptomics_transcript-rna-seq_dea_transcript_rna_seq_20211008.RData")

#############################
# gastro

gastrotimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-gastrocnemius",]
gastrosig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
gastrosigl2fcmat <- matrix(0L,nrow = length(gastrosig),ncol = 8)
rownames(gastrosigl2fcmat) <- gastrosig
colnames(gastrosigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                             "M W1","M W2","M W4","M W8")
for(i in 1:dim(gastrol2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(gastrol2fcmat)[i]
  ourgenemat <- gastrotimewise[gastrotimewise$feature_ID %in% ourgene,]
  gastrol2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  gastrol2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  gastrol2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  gastrol2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  gastrol2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  gastrol2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  gastrol2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  gastrol2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
gastrow8upsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W8"] > 0 & gastrosigl2fcmat[,"M W8"] > 0,])
gastrow8downsig <- rownames(gastrosigl2fcmat[gastrosigl2fcmat[,"F W8"] < 0 & gastrosigl2fcmat[,"M W8"] < 0,])

write.csv(gastrow8upsig,file = "gastrow8upsig.csv")
write.csv(gastrow8downsig,file = "gastrow8downsig.csv")

#############################
# heart

hearttimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-heartcnemius",]
heartsig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
heartsigl2fcmat <- matrix(0L,nrow = length(heartsig),ncol = 8)
rownames(heartsigl2fcmat) <- heartsig
colnames(heartsigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(heartl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(heartl2fcmat)[i]
  ourgenemat <- hearttimewise[hearttimewise$feature_ID %in% ourgene,]
  heartl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  heartl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  heartl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  heartl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  heartl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  heartl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  heartl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  heartl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
heartw8upsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W8"] > 0 & heartsigl2fcmat[,"M W8"] > 0,])
heartw8downsig <- rownames(heartsigl2fcmat[heartsigl2fcmat[,"F W8"] < 0 & heartsigl2fcmat[,"M W8"] < 0,])

write.csv(heartw8upsig,file = "heartw8upsig.csv")
write.csv(heartw8downsig,file = "heartw8downsig.csv")


#############################
# hippo

hippotimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-hippocnemius",]
hipposig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
hipposigl2fcmat <- matrix(0L,nrow = length(hipposig),ncol = 8)
rownames(hipposigl2fcmat) <- hipposig
colnames(hipposigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(hippol2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(hippol2fcmat)[i]
  ourgenemat <- hippotimewise[hippotimewise$feature_ID %in% ourgene,]
  hippol2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  hippol2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  hippol2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  hippol2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  hippol2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  hippol2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  hippol2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  hippol2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
hippow8upsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W8"] > 0 & hipposigl2fcmat[,"M W8"] > 0,])
hippow8downsig <- rownames(hipposigl2fcmat[hipposigl2fcmat[,"F W8"] < 0 & hipposigl2fcmat[,"M W8"] < 0,])

write.csv(hippow8upsig,file = "hippow8upsig.csv")
write.csv(hippow8downsig,file = "hippow8downsig.csv")


#############################
# kidney

kidneytimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-kidneycnemius",]
kidneysig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
kidneysigl2fcmat <- matrix(0L,nrow = length(kidneysig),ncol = 8)
rownames(kidneysigl2fcmat) <- kidneysig
colnames(kidneysigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(kidneyl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(kidneyl2fcmat)[i]
  ourgenemat <- kidneytimewise[kidneytimewise$feature_ID %in% ourgene,]
  kidneyl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  kidneyl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  kidneyl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  kidneyl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  kidneyl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  kidneyl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  kidneyl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  kidneyl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
kidneyw8upsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W8"] > 0 & kidneysigl2fcmat[,"M W8"] > 0,])
kidneyw8downsig <- rownames(kidneysigl2fcmat[kidneysigl2fcmat[,"F W8"] < 0 & kidneysigl2fcmat[,"M W8"] < 0,])

write.csv(kidneyw8upsig,file = "kidneyw8upsig.csv")
write.csv(kidneyw8downsig,file = "kidneyw8downsig.csv")


#############################
# liver

livertimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-livercnemius",]
liversig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
liversigl2fcmat <- matrix(0L,nrow = length(liversig),ncol = 8)
rownames(liversigl2fcmat) <- liversig
colnames(liversigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(liverl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(liverl2fcmat)[i]
  ourgenemat <- livertimewise[livertimewise$feature_ID %in% ourgene,]
  liverl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  liverl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  liverl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  liverl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  liverl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  liverl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  liverl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  liverl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
liverw8upsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W8"] > 0 & liversigl2fcmat[,"M W8"] > 0,])
liverw8downsig <- rownames(liversigl2fcmat[liversigl2fcmat[,"F W8"] < 0 & liversigl2fcmat[,"M W8"] < 0,])

write.csv(liverw8upsig,file = "liverw8upsig.csv")
write.csv(liverw8downsig,file = "liverw8downsig.csv")


#############################
# lung

lungtimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-lungcnemius",]
lungsig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
lungsigl2fcmat <- matrix(0L,nrow = length(lungsig),ncol = 8)
rownames(lungsigl2fcmat) <- lungsig
colnames(lungsigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(lungl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(lungl2fcmat)[i]
  ourgenemat <- lungtimewise[lungtimewise$feature_ID %in% ourgene,]
  lungl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  lungl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  lungl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  lungl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  lungl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  lungl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  lungl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  lungl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
lungw8upsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W8"] > 0 & lungsigl2fcmat[,"M W8"] > 0,])
lungw8downsig <- rownames(lungsigl2fcmat[lungsigl2fcmat[,"F W8"] < 0 & lungsigl2fcmat[,"M W8"] < 0,])

write.csv(lungw8upsig,file = "lungw8upsig.csv")
write.csv(lungw8downsig,file = "lungw8downsig.csv")


#############################
# brown

browntimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-browncnemius",]
brownsig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
brownsigl2fcmat <- matrix(0L,nrow = length(brownsig),ncol = 8)
rownames(brownsigl2fcmat) <- brownsig
colnames(brownsigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(brownl2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(brownl2fcmat)[i]
  ourgenemat <- browntimewise[browntimewise$feature_ID %in% ourgene,]
  brownl2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  brownl2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  brownl2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  brownl2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  brownl2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  brownl2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  brownl2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  brownl2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
brownw8upsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W8"] > 0 & brownsigl2fcmat[,"M W8"] > 0,])
brownw8downsig <- rownames(brownsigl2fcmat[brownsigl2fcmat[,"F W8"] < 0 & brownsigl2fcmat[,"M W8"] < 0,])

write.csv(brownw8upsig,file = "brownw8upsig.csv")
write.csv(brownw8downsig,file = "brownw8downsig.csv")


#############################
# white

whitetimewise <- transcript_rna_seq$timewise_dea[transcript_rna_seq$timewise_dea$tissue %in% "t55-whitecnemius",]
whitesig <- transcript_rna_seq$training_dea[transcript_rna_seq$training_dea$tissue_abbreviation %in% "SKM-GN" & transcript_rna_seq$training_dea$adj_p_value < 0.05,"feature_ID"]
whitesigl2fcmat <- matrix(0L,nrow = length(whitesig),ncol = 8)
rownames(whitesigl2fcmat) <- whitesig
colnames(whitesigl2fcmat) <- c("F W1","F W2","F W4","F W8",
                                "M W1","M W2","M W4","M W8")
for(i in 1:dim(whitel2fcmat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourgene <- rownames(whitel2fcmat)[i]
  ourgenemat <- whitetimewise[whitetimewise$feature_ID %in% ourgene,]
  whitel2fcmat[i,"F W1"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "1w","logFC"]
  whitel2fcmat[i,"F W2"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "2w","logFC"]
  whitel2fcmat[i,"F W4"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "4w","logFC"]
  whitel2fcmat[i,"F W8"] <- ourgenemat[ourgenemat$sex %in% "female" & ourgenemat$comparison_group %in% "8w","logFC"]
  whitel2fcmat[i,"M W1"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "1w","logFC"]
  whitel2fcmat[i,"M W2"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "2w","logFC"]
  whitel2fcmat[i,"M W4"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "4w","logFC"]
  whitel2fcmat[i,"M W8"] <- ourgenemat[ourgenemat$sex %in% "male" & ourgenemat$comparison_group %in% "8w","logFC"]
  
}
whitew8upsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W8"] > 0 & whitesigl2fcmat[,"M W8"] > 0,])
whitew8downsig <- rownames(whitesigl2fcmat[whitesigl2fcmat[,"F W8"] < 0 & whitesigl2fcmat[,"M W8"] < 0,])

write.csv(whitew8upsig,file = "whitew8upsig.csv")
write.csv(whitew8downsig,file = "whitew8downsig.csv")



