suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(CpGassoc))
suppressPackageStartupMessages(library(SmartSVA))
library(isva)
suppressPackageStartupMessages(library(IHW))
suppressPackageStartupMessages(library(ggplot2))
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ggpubr)
library(RColorBrewer)
#library(missMethyl)
library(openxlsx)
library(minfiData)
library(CpGFilter)
library(grid)
library(structSSI)
library(GEOquery)
library(hexbin)
#library(qvalue)
#library(adaptMT)
#library(FDRreg)

Args <- commandArgs()
print(Args)
experiment <- Args[6]

#### load ref data (xReactiveProbes, ann450k, iccref) =====================
load("/home/bailing/projects/ewas/analysis/reffile/ref.RData")

#### set data directory and working directory =============================
outputDir <- "/home/wlw/projects/ewas/analysis"
baseDir <- "/home/wlw/projects/ewas/data"
setwd(paste0(outputDir, "/", experiment))
datadir <- paste0(baseDir, "/", experiment)

beta2m <- function(beta.values){
  log2(beta.values/(1 - beta.values))
}
m2beta <- function(m.values){
  2^m.values/(2^m.values + 1)
}

#### read in data and QC for samples ======================================
if(file.exists(paste0(datadir, "/SampleSheet.csv"))){
  targets <- read.metharray.sheet(datadir, pattern="SampleSheet.csv")
  rgSet <- read.metharray.exp(targets = targets)
  ## QC for samples
  detP <- detectionP(rgSet)
  keep <- colMeans(detP) < 0.01
  if(sum(keep) != nrow(targets)){
    cat("==========", sum(!keep)," sample was filted==========\n")
    rgSet <- rgSet[,keep]
    targets <- targets[keep,]
    detP <- detP[,keep]
  }
  # Normalization
  mSetSq <- preprocessQuantile(rgSet)
}else{
  gset <- getGEO(filename = paste0(datadir, "/", experiment, "_family.soft")) 
  ## beta values and m values
  bVals <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                  ncol = length(gset@header$sample_id))
  rownames(bVals) <- Table(gset@gsms[[1]])$ID_REF
  colnames(bVals) <- gset@header$sample_id
  for (i in 1:length(gset@gsms)) {
    bVals[, i] <- Table(gset@gsms[[i]])[[2]]
  }
  if (all(bVals>=0 & bVals<=1)){
    mVals <- beta2m(bVals)
  } else {
    mVals <- bVals
    bVals <- m2beta(mVals)
  }
  
  ## P values (if exists)
  if (dim(Table(gset@gsms[[1]]))[2]>=3){
    detP <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                   ncol = length(gset@header$sample_id))
    rownames(detP) <- rownames(bVals)
    colnames(detP) <- colnames(bVals)
    for (i in 1:length(gset@gsms)) {
      detP[, i] <- Table(gset@gsms[[i]])[[3]]
    }
  }
 
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets.xlsx"), sheet = experiment)
  ## match phenotype info with values(beta, M, P) and probe info in ann450k
  bVals <- bVals[is.element(rownames(bVals), ann450k$Name), targets$Sample_Name]
  mVals <- mVals[is.element(rownames(mVals), ann450k$Name), targets$Sample_Name]
  if (exists("detP")){
    detP <- detP[is.element(rownames(detP), ann450k$Name), targets$Sample_Name]
    ## QC for samples
    keep <- colMeans(detP,na.rm=TRUE) < 0.01
    if(sum(keep) != nrow(targets)){
      cat("==========", sum(!keep)," sample was filted==========\n")
      targets <- targets[keep,]
      detP <- detP[,keep]
      bVals <- bVals[, keep]
      mVals <- mVals[, keep]
    }
  }
  ## convert to GenomicRatioSet class
  mSetSq <- makeGenomicRatioSetFromMatrix(mat = bVals, pData = targets, what = "Beta")
}

#### QC for probes ========================================================
if (exists("detP")){
  keep <- rowSums(detP < 0.01,na.rm=TRUE) >= ncol(detP) * 0.5
  cat(paste(rep("=",10),collapse = ""),sum(keep==FALSE)," probes have bad quality ",
      paste(rep("=",10),collapse = ""),"\n")
  mSetSqFlt <- mSetSq[keep,]
} else {
  mSetSqFlt <- mSetSq
}

#filter out the probes from the X and Y chromosomes
keep <- !(featureNames(mSetSqFlt) %in%
            ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
#filter out probes that have shown to be cross-reactive
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
#filter out the probes where common SNPs may affect the CpG or single base extension (SBE) site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))

#### get M-values and beta values =========================================
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

#Remove the lines with NA values
is.na.vector <- function(x){
  y <- sum(as.numeric(is.na(x)))>0
  return(y)
}
keep <- !apply(mVals,1,is.na.vector)
mVals <- mVals[keep,]
keep <- !apply(bVals,1,is.na.vector)
bVals <- bVals[keep,]

#### SVA ==================================================================
mod <- model.matrix( ~ targets$Sample_Group)

#for m value
Y.r <- t(resid(lm(t(mVals) ~ targets$Sample_Group)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
sva.res <- smartsva.cpp(mVals, mod, mod0=NULL, n.sv=n.sv)

#### CpGassoc analysis and BH====================================================
cpgRes <- cpg.assoc(bVals, targets$Sample_Group)
cpgResDf <- cpgRes$results
cpgResDf <- na.omit(cpgResDf)

cpgRes.sv <- cpg.assoc(bVals,targets$Sample_Group,as.data.frame(sva.res$sv))
cpgResDf.sv <- cpgRes.sv$results
cpgResDf.sv <- na.omit(cpgResDf.sv)

#Ploting the distribution of p-values
pval_beforesva <- ggplot(NULL,aes(x=cpgResDf$P.value))+ 
  geom_histogram(binwidth=0.01,fill = "#00AFBB",colour="black",size=0.2)+
  theme(panel.background = element_rect(fill = "transparent",colour = NA), #Making the background white
        axis.line=element_line(colour="black"), #Making the axis black
        axis.text= element_text(size=12))+
  scale_y_continuous(expand=c(0,0))+
  labs(title=paste0(experiment,"   P-value_before SVA"))+
  ylab("Count")+
  xlab("P-value")

pval_aftersva <- ggplot(NULL,aes(x=cpgResDf.sv$P.value))+ 
  geom_histogram(binwidth=0.01,fill = "#00AFBB",colour="black",size=0.2)+
  theme(panel.background = element_rect(fill = "transparent",colour = NA), #Making the background white
        axis.line=element_line(colour="black"), #Making the axis black
        axis.text= element_text(size=12))+
  scale_y_continuous(expand=c(0,0))+
  labs(title=paste0(experiment,"   P-value_after SVA"))+
  ylab("Count")+
  xlab("P-value")

pdf.filename <- paste0(experiment,"_pval_hist.pdf")
pdf(pdf.filename,width=8,height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(7,6)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
print(pval_beforesva,vp = vplayout(2:3,2:5))
print(pval_aftersva,vp = vplayout(5:6,2:5))
dev.off()

#Calculation of GIF, pi0, and significant count 
chisq <- qchisq(1-cpgResDf$P.value,1)
lambda <- median(chisq)/qchisq(0.5,1)
lambda <- round(lambda,3)

chisq.sv <- qchisq(1-cpgResDf.sv$P.value,1)
lambda.sv <- median(chisq.sv)/qchisq(0.5,1)
lambda.sv <- round(lambda.sv,3)

pi0 <- pi0est(cpgResDf.sv$P.value,lambda = seq(0.05, 0.95, 0.05),pi0.method = "bootstrap")$pi0
pi0 <- round(pi0,3)
sig.count <- sum(cpgResDf.sv$FDR<0.05)

data.qc <- data.frame(ID = experiment,
                      lambda_before_sva = lambda,
                      lambda_after_sva = lambda.sv,
                      pi0_est = pi0,
                      significant_cpg = sig.count)
write.csv(data.qc,file="QC.csv")

