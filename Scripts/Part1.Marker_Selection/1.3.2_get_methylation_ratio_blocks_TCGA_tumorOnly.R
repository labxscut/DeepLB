# for tumorOnly method
library("data.table")
rm(list = ls())
map_CpG_bin <- function(CpG.probe, cluster.table){
#  probe.chrs <- paste("chr", CpG.probe$V2, sep='')
  probe.chrs <- CpG.probe$V2
  cluster.chrs <- as.character(cluster.table$chr)
  chrs <- unique(probe.chrs)

  probe.start <- CpG.probe$V3
  cluster.start <- cluster.table$start
  cluster.end <- cluster.table$end

  probe.To.cluster <- matrix(nrow=length(probe.chrs), ncol=1)

  for (chr in chrs){
    print(chr)
    probe.ind <- which(probe.chrs == chr)
    tmp <- order(probe.start[probe.ind]) # sort the probe by the start position 
    probe.ind <- probe.ind[tmp] # new order of probe.ind that the correspoinding start coordinates are increasing
    
    cluster.ind <- which(cluster.chrs == chr)

    j <- 1
    cluster.point <- cluster.ind[j]

    for (probe.point in probe.ind){

      flag <- FALSE
      if (j > length(cluster.ind)) probe.To.cluster[probe.point] <- NA
      while(j <= length(cluster.ind) & flag == FALSE){
        if (probe.start[probe.point] >= cluster.start[cluster.point] & probe.start[probe.point] <= cluster.end[cluster.point]) {
          probe.To.cluster[probe.point] <- cluster.point
          flag <- TRUE
        }
        else if (probe.start[probe.point] < cluster.start[cluster.point])
          {
          probe.To.cluster[probe.point] <- NA
          flag <- TRUE
        }
        else {
          j <- j + 1
          cluster.point <- cluster.ind[j]
        }
      }
    }
  }
  #print(i)
  #print(j) 
  return(probe.To.cluster)
}

args <- commandArgs(trailingOnly = TRUE)  
root.dir <- args[1] 
target_tumor <- args[2]
group <- args[3]
select_sample <- args[4]
sampleList <- args[5]
if (group == "TH"){
  prex = "top"
}else if (group == "MH") {
   prex = "sub"
}

## bin annotation
cluster_def_file <- paste(root.dir, "/Result/1.2_define_bins/blocks_D200_S3.txt", sep='')
#df1 <- read.table("/bigData/qingjial/project/MCTD_cfWGBS/CancerDetector/data/1.define_bins/blocks_D200_S3.txt",header=T, stringsAsFactors = F)[, c(2:4, 1)] 
cluster.table <- read.table(cluster_def_file, header=T, stringsAsFactors = F)[, c(2:4, 1)] #, 5)] # chr start end are needed
head(cluster.table )
#chr  start	end	bin_ind no_CpGs
#chr1	566171	567749	

## probe annotation
TCGA.probe.annotation <-  fread(paste0(root.dir,"/Predata/450K/TCGA_450_probe"), header=F)
head(TCGA.probe.annotation)
TCGA.probe.id <- TCGA.probe.annotation$V1

## Mapping probes to each bin (index by the line number, same to cluster.table$marker_index)
probe.To.cluster <- map_CpG_bin(TCGA.probe.annotation, cluster.table)
sum(!is.na(probe.To.cluster)) # 181127 probes can map to our clusters

# map.no.vector <- NULL
# for(i in cluster.table$marker_index){
#   map.no <- length(which(probe.To.cluster == i))
#   detect.no <- cluster.table$block.probe.no[i]
#   map.no.vector <- C(map.no.vector, map.no)
#   if(map.no != detect.no) print(i)
# }

## Calculate the methlation ratio for each bin using the mapped Probes information
input.dir <- paste(root.dir, "/Result/1.1_450K_reformate_data/", sep='') 
output.dir <- paste(root.dir, "/Result/1.3_cal_methylratio_bins/450K/", sep='')

sample.annot.tb <- read.table(paste(root.dir, "/Predata/metadata/",sampleList, sep=''), stringsAsFactors = F, header=F)#SAMPLE ANNOTATION
sample.annot.tb$index <- 1:nrow(sample.annot.tb)
dim(sample.annot.tb)#2343    5
head(sample.annot.tb)

#for (tumor.type in c("coad", 'kirc', 'paad', 'sarc', "lihc", "lusc", "luad", "brca")) {
for (tumor.type in c(target_tumor)) { #,"lihc","paad"
   print(tumor.type)
   #tumor.type = "paad"
   tumor.sample.names = read.table(paste(input.dir, "index_", tumor.type, "_tumor", sep=""), stringsAsFactors = F, header=F)$V2
   sample.names <- tumor.sample.names  
   print(length(sample.names))#lihc:377
  
  # matched.normal.ind <- match(sample.names, normal.sample.names)
   matched.tumor.ind <- match(sample.names, tumor.sample.names)
   
         
  for (status in c("tumor")) {
  file_path = paste(output.dir, tumor.type, "_", status, "_only.tsv", sep='')
  if (file.exists(file_path) && file.size(file_path) > 0) {
    print(paste("file", file_path, "Exists and the size is not 0."))
  } else {
    print(paste("file", file_path, "Does not exist or has a size of 0."))
  
    #status = "tumor"
    matrix.file <- paste(input.dir, tumor.type, '_', status, sep='')
    beta.matrix <- data.matrix(fread(matrix.file, header=F))
    beta.matrix.cluster <- matrix(nrow=nrow(cluster.table), ncol=ncol(beta.matrix))
    no <- 0
    for (i in cluster.table$marker_index) {
      no <- no + 1
      #print(i)
      ind <- which(probe.To.cluster == i)
      #print(length(ind))
      if (length(ind) > 1)
        beta.matrix.cluster[no,] <- apply(beta.matrix[ind, ], 2, mean, na.rm = T)
      else if (length(ind) == 1)
        beta.matrix.cluster[no,] <- beta.matrix[ind, ]
    }
    
    ## sample index defined by Shuli's table
  status.ind <- which(sample.annot.tb$V4 == status)
  sample.annot.tb.subset <- sample.annot.tb[status.ind, ]
  tmp.ind <- match(sample.names, sample.annot.tb.subset$V1)
  sample.index.name <- sample.annot.tb.subset$V1[tmp.ind]
  print(sample.index.name)
  write.table(cbind(sample.index.name, t(beta.matrix.cluster)), file = paste(output.dir, tumor.type, "_", status, "_only.tsv", sep=''), quote=F, sep = "\t", col.names = F, row.names=F)
  }
  tumor_only <- fread(paste(output.dir, tumor.type, "_", status, "_only.tsv", sep=''))
  tumor_only[1:4,1:4]
  purity <- read.table(paste0(root.dir ,"/Predata/tumor_purity/",toupper(tumor.type),"_",prex,select_sample,".txt"))
  #purity <- read.table("/home/yinliang/PROJECT/DeepLB/Predata/tumor_purity/BRCA_top30.txt")
  purity$sample <- sub("-[^-]+$", "", purity$sample)
  tumor_purity <- tumor_only[tumor_only$V1 %in% purity$sample,]
  tumor_purity[1:4,1:4]
  dim(tumor_purity)
  write.table(tumor_purity, file =  paste0(output.dir, tumor.type, "_", prex,select_sample, ".tsv", sep=''),quote=F, sep = "\t", col.names = F, row.names=F)
  }
}

