library(MASS)
library(parallel)
set.seed(2016)
args <- commandArgs(trailingOnly = T)
train <- read.delim(args[1],head=F)
outputDir <- args[2]
nCores = as.numeric(args[3])
dir.create(outputDir,recursive=T,showWarnings=F)

paraEstMoM <- function(meths) {
    mu = mean(meths,na.rm = T)
    var = var(meths,na.rm = T)
    mu[mu == 0] = 1e-5
    mu[mu == 1] = 1-1e-5
    var[var == 0] = 1e-9
    momAlpha = -mu*(var+mu*mu-mu)/var
    momBeta = (mu-1)*(var+mu*mu-mu)/var
    return(c('shape1'=momAlpha,'shape2'=momBeta))
}
#train2 <- train
#train <- train1
lables <- train[,1]
lables <- gsub('_tumor','',lables)
train <- train[,-1]
#train[train == 0] = 1e-9
#train[train == 1] = 1-1e-9
typeTrain <- split(train,lables)
dim(train)
train[1,1:10]
typePara <- list()
for (type in names(typeTrain)) {
    typePara[[type]] <- do.call('cbind',mclapply(as.list(typeTrain[[type]]),paraEstMoM, mc.cores = nCores))
}

alphas = do.call(cbind,lapply(typePara, function(x) x[1,]))
write.table(t(alphas), paste0(outputDir,'/alphas'),sep='\t',row.names=T,col.names=F, quote = F)
betas = do.call(cbind,lapply(typePara, function(x) x[2,]))
write.table(t(betas), paste0(outputDir,'/betas'),sep='\t',row.names=T,col.names=F, quote = F)

#feature selection
# type specific features
normType = 'plasma_background'
tumorTypes = names(typeTrain)
tumorTypes = tumorTypes[!(tumorTypes %in% c(normType))]

# Difference between expect values
indexes <- 1:nrow(alphas)
dist2norm <- sapply(tumorTypes,function(m) {
                    sapply(indexes, function(n) {
                      paras <- c(typePara[[m]][,n], typePara[[normType]][,n])
                      if (any(is.na(paras))) {
                        return(NaN)
                      } else {
                        return(paras[1]/sum(paras[1:2])-paras[3]/sum(paras[3:4]))
                      }
                    })
              })

write.table(t(dist2norm), paste0(outputDir,'/dist2norm'),sep='\t',row.names=T,col.names = F,quote = F)

