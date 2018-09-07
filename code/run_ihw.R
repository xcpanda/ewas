suppressPackageStartupMessages(library(IHW))
suppressPackageStartupMessages(library(ggplot2))
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(openxlsx)
library(grid)
library(structSSI)
library(diptest)
library(qvalue)
library(cowplot)
# library(adaptMT)
# library(FDRreg)

Args <- commandArgs()
print(Args)
experiment <- Args[6]

outputDir <- "/home/bailing/projects/ewas/analysis"
setwd(paste0(outputDir, "/", experiment))

cov.df <- read.csv("Covariates.csv", header = TRUE, row.names = 1)
cov.df$DHS.cv <- str_replace_na(cov.df$DHS.cv, replacement = "FALSE")
covariate.name.list <- colnames(cov.df)
contin.cv <- c("sd.b.cv", "sd.m.cv", "var.b.cv","var.m.cv","mean.b.cv",
               "MAD.cv","DIP.cv","precision.cv","pos.cv", "icc.b", "icc.m")
cate.cv <- c("refgene_pos.cv","CpGlocation.cv","chr.cv", "DHS.cv", "direction.cv")

for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    cov.df[[i]] <- as.numeric(cov.df[[i]])
  }else{
    cov.df[[i]] <- as.factor(cov.df[[i]])
  }
}
cpgResDf.sv <- read.csv("cpgResDf.sv.csv", header = TRUE, row.names = 1)

### method BH ===============================
starttime <- Sys.time()
BH_pvalue <- p.adjust(cpgResDf.sv$P.value, method = "fdr")
endtime <- Sys.time()
BH_time <- as.numeric(endtime - starttime, units = "secs")

### mthod ST ===============================
starttime <- Sys.time()
ST_qvalue <- qvalue(cpgResDf.sv$P.value)$qvalues
endtime <- Sys.time()
ST_time <- as.numeric(endtime - starttime, units = "secs")

### method IHW ===============================
IHW_result <- list()
IHW_time <- list()
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    covariate_type <- "ordinal"  ## for continous covariates
    starttime <- Sys.time()
    result <- ihw(cpgResDf.sv$P.value, cov.df[[i]], alpha = 0.05, covariate_type = covariate_type)
    endtime <- Sys.time()
  }else{
    covariate_type <- "nominal"  ## for category covariates
    starttime <- Sys.time()
    result <- ihw(cpgResDf.sv$P.value, as.factor(cov.df[[i]]), alpha = 0.05, covariate_type = covariate_type)
    endtime <- Sys.time()
  }
  IHW_result[[i]] <- adj_pvalues(result)
  IHW_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
}

write.csv(IHW_result, file = "IHW_result.csv", row.names = FALSE, quote = FALSE)
write.csv(IHW_time, file = "IHW_time.csv", row.names = FALSE, quote = FALSE)

IHW_result <- as.data.frame(IHW_result)
IHW_time <- as.data.frame(IHW_time)

col.palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#949494", "#E6AB02", "#A6761D", 
                 "#A6CEE3", "#1F78B4", "#B2DF8A", "#FFFF00", "#FB9A99", "#E31A1C", "#FDBF6F", 
                 "#FF7F00", "#CAB2D6", "#000000", "#7FFFD4")

## compare probes' number ==================
IHW_probes <- list()
for(i in colnames(IHW_result)){
  IHW_probes[[paste0("ihw_", i)]] <- sum(IHW_result[[i]] < 0.05)
}
probes <- as.data.frame(IHW_probes)
probes$BH <- sum(BH_pvalue < 0.05)
probes$ST <- sum(ST_qvalue < 0.05)

probes <- as.data.frame(t(probes))
colnames(probes) <- "discoveries"
probes$method <- rownames(probes)
p1 <- ggplot(probes, aes(x = method, y = discoveries)) + geom_bar(aes(fill = method), stat = "identity") +
       theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
       scale_fill_manual(values = col.palette) + ggtitle(" Number of Discovered CpGs")
# ggsave("IHW_BH_ST.pdf")

## discoveries VS p-value ====================
colnames(IHW_result) <- paste0("ihw_", colnames(IHW_result))
IHW_result$BH <- BH_pvalue
IHW_result$ST <- ST_qvalue
p.list <- as.list(IHW_result)

getNum <- function(p.value, cut.seq = seq(from = 0, to = 0.05, length.out = 100)){
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

class.name <- names(p.list)
class.n <- length(p.list)
plot_data <- data.frame()
for(i in 1:length(p.list)){
  p.value <- p.list[[i]]
  plot_df <- getNum(p.value)
  plot_df$class <- class.name[i]
  plot_data <- rbind(plot_data, plot_df)
}

p2 <- ggline(data = plot_data, x = "seq", y = "num",
            color = "class", numeric.x.axis = T, plot_type = "l",
            size= 1,xlab = "P-values", ylab = "Discoveries",
            legend.title = "Methods", palette = col.palette)

#ggsave("discoveries_vs_pvalue.pdf", g, width = 8, height = 6)

### compare time ============================
times_df <- IHW_time
colnames(times_df) <- paste0("ihw_", colnames(times_df))
times_df$BH <- BH_time
times_df$ST <- ST_time

times_df <- as.data.frame(t(times_df))
colnames(times_df) <- "run_time"
times_df$method <- rownames(times_df)
p3 <- ggplot(times_df, aes(x = method, y = run_time)) + geom_bar(aes(fill = method), stat = "identity") +
        theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())  + 
        scale_fill_manual(values = col.palette) + ggtitle("Run Time") + ylab("Run Time (secs)")
# ggsave("runtime.pdf")


### save plots
title <- ggtexttable(data.frame(dataset = experiment),
                         theme = ttheme("blank", base_size = 30),rows = NULL)
file_name <- paste0(experiment, "_IHW_ST_BH.pdf")
pdf(file_name, width = 12, height = 16)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6,2)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
print(title, vp = vplayout(1,1))
print(p1, vp = vplayout(2:3,1))
print(p3, vp = vplayout(2:3,2))
print(p2, vp = vplayout(4:6,1:2))
dev.off()

