# for normal blood sample (CUHK)
rm(list = ls())  
#library(readxl)
#options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)  
root.dir <- args[1]
sample_list <- args[2]
#root.dir <- "/home/yinliang/PROJECT/DeepLB"
#sample_list <- "sample_list.normal" 
input.dir <- paste0(root.dir,"/Result/1.3_extract_reads/all_reads_1sites_with_pos_and_sam")
output.dir <- paste0(root.dir,"/Result/1.3_cal_methylratio_bins/WGBS")
dir.create(output.dir,showWarnings = T,recursive = TRUE)

plasma_sample <- read.table(paste0(root.dir,"/Predata/metadata/",sample_list))
beta.matrix.cluster <- matrix(nrow=36342, ncol=nrow(plasma_sample))
sample.i <- 1
for (sample.name in plasma_sample$V1){
 #sample.name <- plasma_sample$V1[1]
  print(sample.name)
  filename <- paste(input.dir, "/", sample.name, ".methy_bins", sep='')
  methy.ratio.tb <- read.table(filename, header=F, skip=1, stringsAsFactors = F)
  #methy.ratio.tb[1:10,]
  beta.matrix.cluster[, sample.i] <- methy.ratio.tb$V2
  sample.i <- sample.i + 1
}
dim(beta.matrix.cluster)
head(beta.matrix.cluster)
write.table(cbind(plasma_sample, t(beta.matrix.cluster)), file = paste(output.dir, "/plasma_background_all.tsv", sep=''), 
quote=F, sep = "\t", col.names = F, row.names=F) 
