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
library(splines)

Args <- commandArgs()
print(Args)
experiment <- Args[6]

outputDir <- "/home/bailing/projects/ewas/analysis"
setwd(paste0(outputDir, "/", experiment))

load("/home/bailing/projects/ewas/analysis/reffile/ref.RData")

cov.df <- read.csv("Covariates.csv", header = TRUE, row.names = 1)
cov.df$DHS.cv <- str_replace_na(cov.df$DHS.cv, replacement = "FALSE")
cov.df$probe_type <- ann450k[rownames(cov.df), "Type"]
colnames(cov.df) <- str_replace_all(colnames(cov.df), pattern = ".cv", replacement = "")
covariate.name.list <- colnames(cov.df)

## classfication of covariates
contin.cv <- c("sd.b", "sd.m", "var.b","var.m","mean.b",
               "MAD","DIP","precision","pos", "icc.b", "icc.m")
cate.cv <- c("refgene_pos","CpGlocation","chr", "DHS", "direction", "probe_type")
statistic.cv <- c("sd.b", "sd.m", "var.b","var.m","mean.b", "icc.b", "icc.m",
                  "MAD","DIP","precision","direction")
CpGs.cv <- c("pos", "refgene_pos","CpGlocation","chr", "DHS", "probe_type")

for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    cov.df[[i]] <- as.numeric(cov.df[[i]])
  }else{
    cov.df[[i]] <- as.factor(cov.df[[i]])
  }
}
cpgResDf.sv <- read.csv("cpgResDf.sv.csv", header = TRUE, row.names = 1)
rownames(cpgResDf.sv) <- cpgResDf.sv$CPG.Labels
cpgResDf.sv <- cpgResDf.sv[rownames(cov.df), ]     

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

### method CAMT ===============================
CAMT_result <- list()
CAMT_time <- list()
source("/home/bailing/projects/ewas/code/CAMT/camt.cor.func.R")
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    X <- ns(cov.df[[i]], df = 6)
    starttime <- Sys.time()
    camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = X, f1.var = X,
                     control.method = 'knockoff')
    endtime <- Sys.time()
    CAMT_result[[i]] <- camt.obj$fdr
    CAMT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }else{
    starttime <- Sys.time()
    camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = cov.df[[i]], f1.var = cov.df[[i]],
                     control.method = 'knockoff')
    endtime <- Sys.time()
    CAMT_result[[i]] <- camt.obj$fdr
    CAMT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }
}

write.csv(CAMT_result, file = "CAMT_result.csv", row.names = FALSE, quote = FALSE)
write.csv(CAMT_time, file = "CAMT_time.csv", row.names = FALSE, quote = FALSE)

CAMT_result <- as.data.frame(CAMT_result)
CAMT_time <- as.data.frame(CAMT_time)

col.palette <- c("BH"="#1B9E77", "ST"="#D95F02", "sd.b"="#7570B3", "sd.m"="#E7298A", 
                 "var.b"="#E6AB02",  "var.m"="#A6761D", "mean.b"="#A6CEE3", 
                 "MAD"="#1F78B4", "DIP"="#B2DF8A", "precision"="#FFFF00", "icc.b"="#FB9A99", 
                 "icc.m"="#E31A1C", "pos"="#949494", "refgene_pos"="#FDBF6F", "chr"="#FF7F00", 
                 "CpGlocation"="#CAB2D6", 
                 "DHS"="#2F4F4F", "direction"="#7FFFD4", "probe_type" ="#FFE4E1")

## compare probes' number ==================
CAMT_probes <- list()
for(i in colnames(CAMT_result)){
  CAMT_probes[[i]] <- sum(CAMT_result[[i]] < 0.05)
}
probes <- as.data.frame(CAMT_probes)
probes$BH <- sum(BH_pvalue < 0.05)
probes$ST <- sum(ST_qvalue < 0.05)

probes <- as.data.frame(t(probes))
colnames(probes) <- "discoveries"
probes$method <- rownames(probes)

group1 <- c("BH", "ST")
group2 <- rownames(probes)[rownames(probes) %in% statistic.cv]
group3 <- rownames(probes)[rownames(probes) %in% CpGs.cv]
method_rank <- c(group1, group2, group3)
probes <- probes[method_rank, ]
# probes <- probes[order(probes$discoveries), ]

write.csv(probes, file = "CAMT_Discoveries.csv", quote = FALSE)
p1 <- ggplot(probes, aes(x = method, y = discoveries)) + geom_bar(aes(fill = method), stat = "identity") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  scale_fill_manual(values = col.palette) + ggtitle(" Number of Discovered CpGs") +
  scale_x_discrete(limits = probes$method)
# p1 <- ggplot(probes, aes(x = method, y = discoveries)) + geom_bar(aes(fill = method), stat = "identity", position = "dodge2") +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5)) + 
#   scale_fill_manual(values = col.palette) + ggtitle(" Number of Discovered CpGs") + scale_x_discrete(limits = probes$method)


## discoveries VS p-value ====================
CAMT_result$BH <- BH_pvalue
CAMT_result$ST <- ST_qvalue
p.list <- as.list(CAMT_result)

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


### compare time ============================
times_df <- CAMT_time
times_df$BH <- BH_time
times_df$ST <- ST_time

times_df <- as.data.frame(t(times_df))
colnames(times_df) <- "run_time"
times_df$method <- rownames(times_df)
write.csv(times_df, file = "run_time.csv", quote = FALSE)
p3 <- ggplot(times_df, aes(x = method, y = run_time)) + geom_bar(aes(fill = method), stat = "identity") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))  + 
  scale_fill_manual(values = col.palette) + ggtitle("Run Time") + ylab("Run Time (secs)") + scale_x_discrete(limits = probes$method)
# ggsave("runtime.pdf")


### save plots
title <- ggtexttable(data.frame(dataset = experiment),
                     theme = ttheme("blank", base_size = 30),rows = NULL)
file_name <- paste0(experiment, "_CAMT_ST_BH.pdf")
pdf(file_name, width = 12, height = 16)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6,2)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
print(title, vp = vplayout(1,1))
print(p1, vp = vplayout(2:3,1))
print(p3, vp = vplayout(2:3,2))
print(p2, vp = vplayout(4:6,1:2))
dev.off()

