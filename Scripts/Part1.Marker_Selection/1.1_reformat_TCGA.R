
#input: XX.tsv, TCGA_450_probe and TCGA_Study_Abbreviations.txt
#output: XX_tumor.tsv and XX_normal.tsv and their index file


rm(list = ls())  
#options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)  
current_path <- args[1]  # current_path = "DeepLB path:/home/yinliang/PROJECT/DeepLB"
tumor_type <- args[2] #"STAD"  #"ESCA" #"STAD" "LIHC"  "STAD" "BRCA"
data.dir <- paste0(current_path,"/Predata/450K/")
output.dir <- paste0(current_path,"/Result/1.1_450K_reformate_data/")
purity.dir <- paste0(current_path,"/Predata/tumor_purity/")
dir.create(output.dir,showWarnings = T,recursive = TRUE)
dir.create(purity.dir,showWarnings = T,recursive = TRUE)

library(dplyr)
library(reshape2)
library(data.table)
library(dplyr)
library(reshape2)
library(InfiniumPurify)

analyze_tumor <-  c(tumor_type)

for (type.short in analyze_tumor){
  print(type.short)
  #type.short <-  "STAD"  #"ESCA" #"STAD" "LIHC"  "STAD" "BRCA"
  normal.index.filename <- paste(output.dir, "index_", tolower(type.short), "_normal", sep='')
  normal.profile.filename <- paste(output.dir, tolower(type.short), "_normal", sep='')
  tumor.index.filename <- paste(output.dir, "index_", tolower(type.short), "_tumor", sep='')
  tumor.profile.filename <- paste(output.dir, tolower(type.short), "_tumor", sep='')


  ### FUNCTION
  short_barcode <- function(sample.id){
    return(paste(unlist(strsplit(sample.id, '-'))[1:4], collapse='.'))
  }

  short_type <- function(sample.id){
    return(paste(unlist(strsplit(sample.id, '-'))[1:3], collapse='-'))
  }



  # input files: TCGA-STAD.methylation450.tsv
  # output files: index_lihc_normal, index_lihc_tumor, lihc_normal, lihc_tumor  
  # ( # ^.^ # )head index_lihc_normal
  #   1       TCGA-BC-A10Q    liver   lihc    normal
  #   2       TCGA-BC-A10R    liver   lihc    normal
  #   3       TCGA-BC-A10S    liver   lihc    normal
  #   4       TCGA-BC-A10T    liver   lihc    normal
  #   
  #   

  #step1: load the data
  probe.annotation.tb <- read.table(paste(data.dir, "TCGA_450_probe", sep=''), header=F) #TCGA_450_probe is sorted by chr and coordinates, using bash script 
  hash.tb <- read.table(paste(data.dir, "TCGA_Study_Abbreviations.txt", sep=''), header=T) #tissue.type, cancer.type
  methy.tb.tmp <- fread(paste(data.dir, "TCGA-", type.short, ".methylation450.tsv", sep=''))
  head(methy.tb.tmp)
  methy.tb.tmp[4,1]
  # prepare the methylation profiles according to the TCGA_450_probe order
  methy.tb <- methy.tb.tmp[match(probe.annotation.tb[,1], methy.tb.tmp$`Composite Element REF`), -1]
  methy.tb[1:4,1:4]
  #step2: prepare the meta table
  samples.id.tmp <- colnames(methy.tb.tmp)[2:dim(methy.tb.tmp)[2]]#"TCGA.NK.A5CR.01"
  samples.id <- unlist(lapply(samples.id.tmp, short_type))#"TCGA.NK.A5CR"
  samples.code <- unlist(lapply(samples.id.tmp, substr, 14, 15)) #"01"

  cancer.type <- rep(tolower(type.short), length(samples.id))
  tissue.short <- hash.tb[which(hash.tb$Cancer.type == tolower(type.short)), 1]
  tissue.type <- rep(tissue.short, length(samples.id))
  meta.tb <- data.frame(samples.id, tissue.type, cancer.type)
  #meta.tb
  #step3: output the normal files
  normal.index <- which(samples.code == "11")
  normal.index.tb <- meta.tb[normal.index,]
  normal.index.tb <- cbind(seq(1, length(normal.index)), normal.index.tb, rep("normal", length(normal.index)))
  ncol(methy.tb)
  normal.methy.tb <- methy.tb[, ..normal.index]
  #methy.tb[1,normal.index]
  #step4: output the tumor files
  tumor.index <- which(samples.code == "01")  
  tumor.index.tb <- meta.tb[tumor.index,]
  tumor.index.tb <- cbind(seq(1, length(tumor.index)), tumor.index.tb, rep("tumor", length(tumor.index)))
  tumor.methy.tb <- methy.tb[, ..tumor.index] 

  dim(normal.methy.tb)
  dim(tumor.methy.tb)

  write.table(normal.index.tb, file=normal.index.filename, quote=F, sep = "\t", col.names = F, row.names=F)
  write.table(normal.methy.tb, file=normal.profile.filename, quote=F, sep = "\t", col.names = F, row.names=F)
  write.table(tumor.index.tb, file=tumor.index.filename, quote=F, sep = "\t", col.names = F, row.names=F)
  write.table(tumor.methy.tb, file=tumor.profile.filename, quote=F, sep = "\t", col.names = F, row.names=F)

  beta_data <- tumor.methy.tb
  head(beta_data)
  dim(probe.annotation.tb)
  dim(beta_data)
  beta_data[,1]
  
  beta <- as.matrix(beta_data) 
  colnames(beta)
  rownames(beta) <- probe.annotation.tb$V1
  purity <- getPurity(tumor.data = beta[,1:ncol(beta)],normal.data = NULL,tumor.type= type.short)
  purity <- as.data.frame(purity)
  #View(purity)
  dim(purity)
  write.csv(file = paste0(purity.dir,type.short,"_purity.csv"),purity)
}




