rm(list=ls())
options(stringsAsFactors=F)

## This function define the blocks
find_block <- function(chr.probe.coordinates, window.size=100, min.probe.number=3, chr.probe.names, chr.probe.genes){
  
  start.coordinate <- NULL
  end.coordinate <- NULL
  probe.no <- NULL
  probe.list <- NULL
  gene.list <- NULL
  gene.no <- NULL
  
  block.ind <- 1
  i <- 2
  while(i <= length(chr.probe.coordinates)) {
    #print(i)
    if (chr.probe.coordinates[i] <= chr.probe.coordinates[i-1] + 2*window.size) { #adjacent probes are within the window, merge them
      block.ind <- c(block.ind, i)
      i <- i + 1
    }
    else if (length(block.ind) >= min.probe.number) {
        # check the CpG Clusters among the probes
        print(block.ind)
        start.coordinate <- c(start.coordinate, chr.probe.coordinates[block.ind[1]] - window.size)
        end.coordinate <- c(end.coordinate, chr.probe.coordinates[block.ind[length(block.ind)]] + window.size)
        probe.no <- c(probe.no, length(block.ind))
        probe.list <- c(probe.list, paste(chr.probe.names[block.ind], collapse=','))
        
        gene.str.tmp <- paste(chr.probe.genes[block.ind], collapse=',') # one CpG might have multiple genes
        gene.str.tmp.2 <- setdiff(unique(unlist(strsplit(gene.str.tmp, ","))), ".") #unique them, then remove .
        gene.str <- paste(gene.str.tmp.2, collapse = ',') # then collapse as a string seperated by ','
        gene.list <- c(gene.list, gene.str)
        gene.no <- c(gene.no, length(gene.str.tmp.2))
        
        block.ind <- i # restart a block
        i <- i + 1
    } else {
      block.ind <- i # restart a block
      i <- i + 1
    }
  }
  block.list = data.frame(start=start.coordinate, end=end.coordinate, probe.no = probe.no, probe.list, gene.list, gene.no)
  return(block.list)
}

require(gplots)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)  
current_path <- args[1] 
window.size <- args[2] #100
min.probe.number <- args[3] #3

output.dir <- paste0(current_path,"/Result/1.2_define_bins/")
dir.create(output.dir)
CpG.probe.tmp <- read.table(paste0(current_path,"/Predata/450K/illuminaMethyl450_hg38_GDC"), header=F, stringsAsFactors = F) 
CpG.probe <- CpG.probe.tmp %>% filter(V3 != '*')
probe.reorder.ind <- order(CpG.probe[,3], CpG.probe[,4]) # sort 450k probes by chr+coordinate
CpG.probe.sorted <- CpG.probe[probe.reorder.ind, ]

chr.names <- setdiff(unique(CpG.probe[,3]), '*')

## for each chr, using "sliding window" method to merge the adjacent probes to block (at least 3 probes)

block.chr <- NULL
block.start <- NULL
block.end <- NULL
block.probe.no <- NULL
block.probe.list <- NULL
block.gene.list <- NULL
block.gene.no <- NULL

for (chr in chr.names){
  print(chr)
  ind <- which(CpG.probe.sorted[, 3] == chr)
  #print(length(ind))
  chr.probe.coordinates <- CpG.probe.sorted[ind, 4]
  chr.probe.names <- CpG.probe.sorted[ind, 1]
  chr.probe.genes <- CpG.probe.sorted[ind, 2]
  
  chr.block.list <- find_block(chr.probe.coordinates, window.size, min.probe.number, chr.probe.names, chr.probe.genes)
  #block.size <- block.list$end - block.list$start 
  print(dim(chr.block.list)[1])
  block.chr <- c(block.chr, rep(chr, dim(chr.block.list)[1]))
  block.start <- c(block.start, chr.block.list$start)
  block.end <- c(block.end, chr.block.list$end)
  block.probe.no <- c(block.probe.no, chr.block.list$probe.no)
  
  block.probe.list <- c(block.probe.list, chr.block.list$probe.list)
  block.gene.list <- c(block.gene.list, chr.block.list$gene.list)
  block.gene.no <- c(block.gene.no, chr.block.list$gene.no)
}

block.df <- data.frame(marker_index=seq(1, length(block.start)), chr=block.chr, start=block.start, end=block.end, block.probe.no=block.probe.no, 
                       block.probes=block.probe.list, block.genes=block.gene.list, block.gene.no=block.gene.no)

block.size <- block.df$end - block.df$start
summary(block.size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 204     352     438     479     559    3724 

round(table(block.df$block.gene.no) / dim(block.df)[1] * 100, digits = 3)
# 0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22 
# 9.471 61.562 24.272  3.858  0.561  0.102  0.036  0.022  0.014  0.008  0.011  0.003  0.014  0.008  0.008  0.008  0.008  0.006  0.008  0.006  0.006  0.003  0.006 

summary(block.df$block.probe.no)
sum(block.df$block.probe.no)

write.table(block.df, file=paste(output.dir, "blocks_D", window.size*2, "_S", min.probe.number, "_full_information.txt", sep=""), quote=F, sep = "\t", col.names = T, row.names=F)      
      
#output the simple format for add_complementary_bins.py
write.table(block.df[,1:4], file=paste(output.dir, "blocks_D", window.size*2, "_S", min.probe.number, ".txt", sep=""), quote=F, sep = "\t", col.names = T, row.names=F) 


