# select different sample by tumor purity
rm(list = ls())  
args <- commandArgs(trailingOnly = TRUE) 
root_dir <-  args[1] 
tumor <- args[2]
group <- args[3] 
select_number <- args[4] 
out_dir <- paste0(root_dir,"/Predata/tumor_purity/")
purity <- read.csv(paste0(root_dir,"/Predata/tumor_purity/",tumor,"_purity.csv"))
purity$sample <- gsub("\\.", "-", sub("^(TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}).*$", "\\1", purity$X))

library(caret)
#for TH group 

if (group == "TH"){
  select_sample <- purity[order(purity$purity, decreasing = TRUE), ][1:select_number, ]
  print(paste0("cut_off: ", round(select_sample[nrow(select_sample),]$purity,2)))
  write.table(select_sample,paste0(out_dir,tumor,"_","top",select_number,".txt"))
  png(filename =paste0(out_dir,tumor,"_","top",select_number,".png"),width = 800,height = 600, units = "px", res = 120)
  hist(purity$purity, breaks = 30, col = rgb(0,0,1,0.5), main = "Original vs Sampled Distribution", xlab = "Purity", ylab = "Frequency")
  hist(select_sample$purity, breaks = 30, col = rgb(1,0,0,0.5), add = TRUE)
  legend("topright", legend = c("Original", "Sampled"), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
  dev.off()
}else if (group == "MH") {
  original_sample <- purity
  p <- as.numeric(select_number)/nrow(purity)
  sample_indices <- createDataPartition(original_sample$purity, p = p, list = FALSE)
  sampled_sample <- original_sample[sample_indices,]
  dim(sampled_sample)
  write.table(sampled_sample,paste0(out_dir,tumor,"_","sub",select_number,".txt"))
  png(filename =paste0(out_dir,tumor,"_","sub",select_number,".png"),width = 800,height = 600, units = "px", res = 120)
  hist(original_sample$purity, breaks = 30, col = rgb(0,0,1,0.5), main = "Original vs Sampled Distribution", xlab = "Purity", ylab = "Frequency")
  hist(sampled_sample$purity, breaks = 30, col = rgb(1,0,0,0.5), add = TRUE)
  legend("topright", legend = c("Original", "Sampled"), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
  dev.off()
}
