library("zinbwave")
library("DESeq2")
library("edgeR")
library("variancePartition")

setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/sectional-temporal-celltype/Simulations/simulation.res2")
setwd("C:/Users/YUXUANXUAN/Dropbox/2021_Researches/sectional-temporal-celltype/Simulations/simulation.res2")

sim.SCRIP8<-readRDS("sim.SCRIP8.rds")
sim<-sim.SCRIP8


counts <- counts(sim)
colnames(counts) <- paste0("Cell",1:ncol(counts))
rownames(counts) <- paste0("Gene",1:nrow(counts))

sum=apply(counts,1,sum)
counts=counts[which(sum>0),]
# zinb-WAVE weights
group <- factor(sim$group.comp$group)
time=factor(sim$group.comp$time)
coldata<-data.frame(condition=group,
                    individual=factor(sim$group.comp$sample),
                    time=time)
rownames(coldata) <- colnames(counts)
example.data<-list(counts=counts,coldata=coldata)
setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/sectional-temporal-celltype/Rcode/modDream/data")
save(example.data,file="example.data.rda")



setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/sectional-temporal-celltype/Rcode/modDream/data")
example.data<-readRDS("example.data.RDS")
counts  <- example.data$counts
coldata <- example.data$coldata
fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition + time)
zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, observationalWeights = TRUE)
weights <- assay(fluidigm_zinb, "weights")
counts[1:5,1:5]
weights[1:5,1:5]
saveRDS(weights,"ZINB.WaVE.weight.RDS")
### Dream ###
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)


form <- ~condition + time + (1 + time |individual)
vobjDream = voomWithDreamWeights( d0, form, coldata )
vobjDream$weights[1:5,1:5]
# vobjDream.weight<-vobjDream
# vobjDream.weight$weights <- weights*vobjDream.weight$weights
setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/sectional-temporal-celltype/Rcode/modDream/data")
setwd("C:/Users/YUXUANXUAN/Dropbox/2021_Researches/sectional-temporal-celltype/Rcode/modDream/data")

saveRDS(vobjDream.weight,"vobjDream.weight.RDS")
vobjDream.weight<-readRDS("vobjDream.weight.RDS")
vobjDream.weight[["weights"]][1:5,1:5]
L = getContrast( vobjDream, form, coldata, c("condition1", "condition2"))
fitmm.weight = modDream( vobjDream.weight, form, coldata, L)

dream.weight <- topTable(fitmm.weight, coef="L1", number=nrow(counts) )
dream.weight <- dream.weight[rownames(counts),]



library("DESeq2")

library("zinbwave")
library("edgeR")

library("lme4")      # lmerControl
library("iterators") # nextElem function
library("pbkrtest") 
library("foreach")
library("variancePartition")
library("stats")
library("plyr")
library("dplyr")



# library("BiocParallel")
# library("Biobase")
# library("splatter")
# library("scater")
# 
# library("checkmate")



#' @importFrom lmerControl  lme4
#' @importFrom nextElem     iterators
#' @importFrom pbkrtest     get_SigmaG
#' @import     foreach
#' @importFrom variancePartition     colinearityScore
# ========================================================== #
library(devtools)
install_github("xuanxuanyu-bios/SCUTA",force=TRUE)
library("SCUTA")
data(example.data)
data(fit.res)
library("zinbwave")
library("DESeq2")
library("edgeR")
library("variancePartition")

counts  <- example.data$counts
coldata <- example.data$coldata
coldata$time<-as.numeric(coldata$time)
fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition + time)


zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000, observationalWeights = TRUE)
weights <- assay(fluidigm_zinb, "weights")


d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
form <- ~ condition + time + (1 + time |individual)
vobjDream = voomWithDreamWeights( d0, form, coldata )
vobjDream.weight<-vobjDream
vobjDream.weight$weights <- weights*vobjDream.weight$weights


fit.scMLLM     <- scMLLM( vobjDream.weight, form, coldata)
fit.scMLLM.res <- topTable(fit.scMLLM, coef="condition2", number=nrow(counts) )
head(fit.scMLLM.res)
