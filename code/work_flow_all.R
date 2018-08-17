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

Args <- commandArgs()
print(Args)
experiment <- Args[6]

#### load ref data (xReactiveProbes, ann450k, iccref) =====================
load("/home/bailing/projects/ewas/analysis/reffile/ref.RData")

#### set data directory and working directory =============================
outputDir <- "/home/bailing/projects/ewas/analysis"
baseDir <- "/home/bailing/projects/ewas/data"


setwd(paste0(outputDir, "/", experiment))
datadir <- paste0(baseDir, "/", experiment)

#### read in data and QC for samples ======================================

if(file.exists(paste0(datadir, "/SampleSheet.csv"))){
  ### starting from raw data ----------------------------------------------
  targets <- read.metharray.sheet(datadir, pattern="SampleSheet.csv")
  rgSet <- read.metharray.exp(targets = targets)
  
  ## QC for samples
  
  # identifies failed positions
  detP <- detectionP(rgSet)
  
  # fig 1.1 qcReport
  qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Sample_Group, 
           pdf="1.1_qcReport_before_filter.pdf")
  
  keep <- colMeans(detP) < 0.01
  if(sum(keep) != nrow(targets)){
    cat("==========", sum(!keep)," sample was filted==========\n")
    rgSet <- rgSet[,keep]
    targets <- targets[keep,]
    detP <- detP[,keep]
  }
  pal <- brewer.pal(8,"Dark2")
  pdf("2.1_barplot_meanPvalue_by_sample.pdf")
  barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
          cex.names=0.8, ylab="Mean detection p-values")
  abline(h=0.05,col="red")
  dev.off()
  
  # preprocess data
  mSetSq <- preprocessQuantile(rgSet)
}else{
  ### starting form GEO format file ---------------------------------------
  gset <- getGEO(filename = paste0(baseDir, "/", experiment, "/", experiment, "_family.soft")) 
  ## b values data and p value
  bVals <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                  ncol = length(gset@header$sample_id))
  rownames(bVals) <- Table(gset@gsms[[1]])$ID_REF
  colnames(bVals) <- gset@header$sample_id
  detP <- bVals
  for (i in 1:length(gset@gsms)) {
    bVals[, i] <- Table(gset@gsms[[i]])[[2]]
    detP[, i] <- Table(gset@gsms[[i]])[[3]]
  }
  ## m values
  beta2m <- function(beta.values){
    log2(beta.values/(1 - beta.values))
  }
  mVals <- beta2m(bVals)
  
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets.xlsx"), sheet = experiment)
  
  ## match phenotype info with values(beta, M, P) and probe info in ann450k
  bVals <- bVals[is.element(rownames(bVals), ann450k$Name), targets$Sample_Name]
  mVals <- mVals[is.element(rownames(mVals), ann450k$Name), targets$Sample_Name]
  detP <- detP[is.element(rownames(detP), ann450k$Name), targets$Sample_Name]
  
  ## QC for samples
  keep <- colMeans(detP) < 0.01
  if(sum(keep) != nrow(targets)){
    cat("==========", sum(!keep)," sample was filted==========\n")
    targets <- targets[keep,]
    detP <- detP[,keep]
    bVals <- bVals[, keep]
    mVals <- mVals[, keep]
  }
  pal <- brewer.pal(8,"Dark2")
  pdf("2.1_barplot_meanPvalue_by_sample.pdf")
  barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
          cex.names=0.8, ylab="Mean detection p-values")
  abline(h=0.05,col="red")
  dev.off()
  
  ## convert to GenomicRatioSet class
  mSetSq <- makeGenomicRatioSetFromMatrix(mat = bVals, pData = targets, what = "Beta")
}

#### QC for probes ========================================================
keep <- rowSums(detP < 0.01) >= ncol(detP) * 0.5
cat(paste(rep("=",10),collapse = ""),summary(keep)[2]," probes have bad quality ",
    paste(rep("=",10),collapse = ""),"\n")
mSetSqFlt <- mSetSq[keep,]

#filter out the probes from the X and Y chromosomes
keep <- !(featureNames(mSetSqFlt) %in%
            ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]

#filter out probes that have shown to be cross-reactive
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] 

#filter out the probes where common SNPs may affect the CpG or single base extension (SBE) site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))

# fig 2.2 qcReport (The workflow for GEO format data doesn't have RGChannelSet object)
# qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Sample_Group,
#          pdf="2.2_qcReport_after_filter.pdf")

#### get M-values and beta values =========================================
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

m2beta <- function(m.values){
  2^m.values/(2^m.values + 1)
}
beta2m <- function(beta.values){
  log2(beta.values/(1 - beta.values))
}

# fig 3
pdf("3_densityPlot_after_nor.pdf")
densityPlot(bVals, sampGroups=targets$Sample_Group,
            main="Density Plot", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# fig 4
pdf("4.1_MDSplot_by_SampleGroup.pdf")
plotMDS(mVals, top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

#### SVA ==================================================================
mod <- model.matrix( ~ targets$Sample_Group)

#for m value
Y.r <- t(resid(lm(t(mVals) ~ targets$Sample_Group)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
sva.res <- smartsva.cpp(mVals, mod, mod0=NULL, n.sv=n.sv)
modsv <- cbind(mod, sva.res$sv)
fitsv <- lm.fit(modsv, t(mVals))
mVals.sva <- t(fitsv$fitted.values)

#for b value
bVals.sva <- m2beta(mVals.sva)

#### ICC analysis by cpgfilter ============================================

# only calculate icc when replicates exist
if(length(unique(targets$Rep_Design)) != length(targets$Rep_Design)){
  covariate.name.list <- c("icc.m", "icc.b", "var.m","var.b","icc.ref","Relation_to_Island","pos","chr")
  icc.m <- CpGFilterICC(mVals.sva, targets$Rep_Design, logit.transform = FALSE)
  icc.b <- CpGFilterICC(bVals.sva, targets$Rep_Design, logit.transform = FALSE)
}else{
  covariate.name.list <- c("var.m","var.b","icc.ref","Relation_to_Island","pos","chr")
}

#### calculate variance ===================================================

var.m <- rowSds(mVals.sva)
var.b <- rowSds(bVals.sva)

#### CpGassoc analysis ====================================================

cpgRes <- cpg.assoc(bVals.sva, targets$Sample_Group)

# fig 5
pdf("5_CpGassoc_analysis.pdf")
plot(cpgRes)
dev.off()

#### IHW ==================================================================

anno.item <- c("chr", "pos", "strand", "Name", "AddressA", "AddressB", "Type", "NextBase",
               "Color", "Probe_rs", "Probe_maf", "CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf",
               "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name",
               "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "Phantom", "DMR",
               "Enhancer", "HMM_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS")
categorical.cov <- c("chr","strand","Type","Color","Relation_to_Island","NextBase")

# ihw Diagnostic plots for covariate
diagno_plot <- function(data = ihwResDf, covariate.name, phenotype = "Sample_Group", remark = "-", method.name){
  ihw.p1 <- gghistogram(data, x = "pvalue",fill = "#00AFBB", binwidth = 0.01,
                        main ="Original p-value")
  
  ihw.p2 <- gghistogram(data, x = "pvalue",fill = "grey25", binwidth = 0.01,
                        color = "grey25",main ="Original p-value by ihw group") +
    facet_wrap( ~ ihwGroup, nrow = 2)
  
  ihw.p3 <- ggplot(data,aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100)
  
  ihw.p4 <- ggplot(data, aes(x = pvalue, col = ihwGroup)) + stat_ecdf(geom = "step") 
  
  ihw.title <- ggtexttable(data.frame(PhenoType = phenotype,
                                      CovariateName = covariate.name,
                                      Method = method.name,
                                      Remark = remark),
                           theme = ttheme("blank"),rows = NULL)
  pdf.filename <- paste0(method.name, "_Diagnostic_plots_for_", covariate.name,".pdf")
  pdf(pdf.filename, width = 8, height = 12)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(8,2)))
  vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
  print(ihw.title, vp = vplayout(1,1))
  print(ihw.p1, vp = vplayout(2:3,1:2))
  print(ihw.p2, vp = vplayout(4:6,1:2))
  print(ihw.p3, vp = vplayout(7:8,1))
  print(ihw.p4, vp = vplayout(7:8,2))
  dev.off()
  cat(covariate.name,method.name,"plot done\n")
}

ihwResPlot <- function(covariate.name, outputDir){
  # get covariate values
  print(covariate.name)
  cpgResDf <- cpgRes$results
  if( is.element(covariate.name, anno.item) ){
    ann450k.iterm <- ann450k[, which(colnames(ann450k)==covariate.name)]
    if(!is.element(class(ann450k.iterm), c("integer","numeric")) ){
      ann450k.iterm <- as.factor(ann450k.iterm)
    }
    cpgResDf$covariate <- ann450k.iterm[match(cpgResDf$CPG.Labels, rownames(ann450k))]
  }else if( is.element(covariate.name, c("icc.m","icc.b","var.m","var.b")) ){
    cpgResDf$covariate <- get(covariate.name)
  }else if(covariate.name == "icc.ref"){
    cpgResDf$covariate <- iccref$rangeEPIC[match(cpgResDf$CPG.Labels, iccref$CpGname)]
  }else{
    stop("wrong covariate.name")
  }
  cpgResDf <- na.omit(cpgResDf)
  
  # calculate the adjust p-values base on the type of variable
  if(is.element(covariate.name, categorical.cov)){
    #covariate is categorical variable, use GBH
    method.name <- "GBH"
    unadj_p <- cpgResDf$P.value
    names(unadj_p) <- cpgResDf$CPG.Labels
    gbhRes <- Adaptive.GBH(unadj_p, cpgResDf$covariate, alpha = 0.05, method = "storey")
    gbhResDf <- gbhRes@p.vals
    ihwResDf <- gbhResDf #just unify the R variable name.
    colnames(ihwResDf)[1:2] <- c("pvalue", "adj_pvalue")
    ihwResDf$covariate <- ihwResDf$group
    
    gbhRes.plot <- plot(gbhRes)
    ggsave(paste0("6_gbh_result_",covariate.name,".pdf"),
           gbhRes.plot, width = 8, height = 8)
  }else{
    #covariate is continuous variable, use IHW
    method.name <- "IHW"
    ihwRes <- ihw( P.value ~ covariate, alpha = 0.05, data = cpgResDf)
    ihwResDf <- as.data.frame(ihwRes)
    rownames(ihwResDf) <- cpgResDf$CPG.Labels
    ihwResDf <- ihwResDf[order(ihwResDf$pvalue),]
    
    ihwRes.ggplot2 <- plot(ihwRes)
    ggsave(paste0("6_ihw_result_",covariate.name,".pdf"),
           ihwRes.ggplot2,width = 8,height = 8)
  }
  
  ihwResDf$ihwGroup <- groups_by_filter(ihwResDf$group, 8)
  write.csv(ihwResDf, paste0(method.name,"_table_for_",covariate.name,".csv"),row.names = T)
  cat(covariate.name, method.name, "calc done\n")
  
  # ihw Diagnostic plots for covariate
  diagno_plot(data = ihwResDf, covariate.name = covariate.name, method.name = method.name)
  
}

# RUN IHW
for (covariate.name in covariate.name.list) {
  ihwResPlot(covariate.name = covariate.name,
             outputDir = outputDir)
}

#### compare results of different methods =================================

### compare discovery numbers 

## counts for probes under different p-value range
getNum <- function(p.value, cut.seq = seq(from = 0, to = 0.05, length.out = 50)){
  #to figure out how many p.value < each cut.seq
  #p.value required a number vector
  #output is a data.fram which colnames are "seq" and "num"
  res <- data.frame(seq = cut.seq, num = NA)
  for (i in 1:length(cut.seq)) {
    s <- cut.seq[i]
    res[i,2] <- sum(p.value < s)
  }
  return(res)
}

## extract adjust p-values for each covariate
get_pvalue <- function(covariate.name){
  if(covariate.name %in% categorical.cov){
    adj_pvalue <- read.csv(paste0("GBH_table_for_", covariate.name, ".csv"))$adj_pvalue
  }else{
    adj_pvalue <- read.csv(paste0("IHW_table_for_", covariate.name, ".csv"))$adj_pvalue
  }
  return(adj_pvalue)
}
p.list <- list()
for(i in covariate.name.list){
  p.list[[i]] <- get_pvalue(covariate.name = i)
}

cpgResDf <- cpgRes$results
p.list[["raw"]] <- cpgResDf$P.value
p.list[["fdr"]] <- cpgResDf$FDR


## Num-Pvalue plot
PvalueNumPlot <- function(p.list, cut.seq = seq(from = 0, to = 0.05, length.out = 50),
                          output.data.file = NULL, output.plot.file = NULL,
                          col.palette = brewer.pal(12, "Paired")){
  
  #input: a list included every p.value of class you want plot
  #cut.seq is same as cut.seq in getNum()
  #output.data.file/output.plot.file: path to file
  class.name <- names(p.list)
  class.n <- length(p.list)
  plot_data <- data.frame()
  for(i in 1:length(p.list)){
    p.value <- p.list[[i]]
    plot_df <- getNum(p.value, cut.seq = cut.seq)
    plot_df$class <- class.name[i]
    plot_data <- rbind(plot_data, plot_df)
  }
  
  if(!is.null(output.data.file)){
    #save data as a csv file.
    num.data <- cbind(data.frame(seq = cut.seq), unstack(plot_data, num ~ class))
    write.csv(num.data, output.data.file, row.names = F, quote = F)
  }
  
  if(!is.null(output.plot.file)){
    #plot and save
    #g <- ggplot(plot.data, aes(x = seq, col = class, y = num))+geom_line()
    g <- ggline(data = plot_data, x = "seq", y = "num",
                color = "class", numeric.x.axis = T, plot_type = "l",
                size= 1,xlab = "P-values", ylab = "Discoveries",
                legend.title = "Methods", palette = col.palette[as.factor(class.name)])
    ggsave(output.plot.file, g, width = 8, height = 6)
  }
}
PvalueNumPlot(p.list, output.plot.file = "num-p.pdf")

### concordance rate between FDR and other methoads

## extract probes with p-value < 0.05 for each covariate
get_probes <- function(covariate.name){
  if(covariate.name %in% categorical.cov){
    df <- read.csv(paste0("GBH_table_for_", covariate.name, ".csv"), row.names = 1)
  }else{
    df <- read.csv(paste0("IHW_table_for_", covariate.name, ".csv"), row.names = 1)
  }
  index <- df$adj_pvalue < 0.05
  probes <- rownames(df)[index]
  return(probes)
}

probes.list <- list()
for(i in covariate.name.list){
  probes.list[[i]] <- get_probes(covariate.name = i)
}

cpgResDf <- cpgRes$results
probes.list[["raw"]] <- cpgResDf$CPG.Labels[cpgResDf$P.value < 0.05]
probes.list[["fdr"]] <- cpgResDf$CPG.Labels[cpgResDf$FDR < 0.05]

ratio.df <- data.frame(methods = rep(NA, (length(probes.list)-2)), concordance.rate = NA)
for(i in 1:(length(probes.list) - 2)){
  index <- probes.list[[i]] %in% probes.list[["fdr"]]
  proportion <- sum(index)/length(probes.list[[i]])
  ratio.df[i, "methods"] <- names(probes.list)[i]
  ratio.df[i, "concordance.rate"] <- proportion
}
ratio.df <- ratio.df[order(ratio.df$concordance.rate),]
p <- ggbarplot(ratio.df, "methods", "concordance.rate", label = TRUE, lab.pos = "out", lab.nb.digits = 3, 
          color = "methods", fill = "methods")
ggsave("concordance_rate_to_fdr.pdf", plot = p, width = 8, height = 6)

#### save objects as RData ================================================
save.image("analysis.RData")
