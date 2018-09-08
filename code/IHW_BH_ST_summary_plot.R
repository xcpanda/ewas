setwd("/home/bailing/projects/ewas/analysis")
dataset <- c("GSE102468", "GSE104293", "GSE109914", "GSE114753", "GSE50660", "GSE53045",
             "GSE40279", "GSE40576", "GSE42861", "GSE87095", "GSE49149", "GSE50798", 
             "GSE59250", "GSE60787", "GSE63704", "GSE64380", "GSE66210", "GSE73515", 
             "GSE76938", "GSE95049")
result <- data.frame()
for(s in dataset){
  file_path <- paste0(s, "/Discoveries.csv")
  if(file.exists(file_path)){
    df <- read.csv(file_path, header = TRUE)
    df <- df[, c(2,3)]
    df$dataset <- s
    df$rank <- rank(df$discoveries, ties.method = "average")
    result <- rbind(result, df)
  }
}

col.palette <- c("BH"="#1B9E77", "ST"="#D95F02", "sd.b"="#7570B3", "sd.m"="#E7298A", 
                 "var.b"="#E6AB02",  "var.m"="#A6761D", "mean.b"="#A6CEE3", 
                 "MAD"="#1F78B4", "DIP"="#B2DF8A", "precision"="#FFFF00", "icc.b"="#FB9A99", 
                 "icc.m"="#E31A1C", "pos"="#949494", "refgene_pos"="#FDBF6F", "chr"="#FF7F00", 
                 "CpGlocation"="#CAB2D6", 
                 "DHS"="#2F4F4F", "direction"="#7FFFD4", "probe_type" ="#FFE4E1")
method_rank <- c("BH", "ST", "sd.b", "sd.m", "var.b", "var.m", "mean.b", "MAD", "DIP", 
                 "precision", "icc.m", "icc.b", "direction", "pos", "refgene_pos",
                 "CpGlocation", "chr", "DHS", "probe_type")
library(ggpubr)
p1 <- ggboxplot(result, x = "method", y = "rank",
          fill = "method",
          add = "dotplot",
          order = method_rank,
          add.params = list(color = "black", fill = "black", dotsize = 0.3),
          palette = col.palette,
          title = "Rank by Covariate"
          ) + theme_pubr(x.text.angle = 60)

### run_time boxplot =============================================
result <- data.frame()
for(s in dataset){
  file_path <- paste0(s, "/run_time.csv")
  if(file.exists(file_path)){
    df <- read.csv(file_path, header = TRUE)
    df <- df[, c(2,3)]
    df$dataset <- s
    result <- rbind(result, df)
  }
}
result$lg_runtime <- log10(result$run_time)
p2 <- ggboxplot(result, x = "method", y = "lg_runtime",
          fill = "method",
          add = "dotplot",
          order = method_rank,
          add.params = list(color = "black", dotsize = 0.3, fill = "black"),
          palette = col.palette,
          ylab = "log10(run time)/secs",
          title = "Run Times Across 20 Datasets"
          ) + theme_pubr(x.text.angle = 60)
p2 <- ggpar(p2, ylim = c(-1.5, 3))

## save plots
library(grid)
pdf("Summary_20datasets.pdf", width = 12, height = 18)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1)))
vplayout = function(x,y)viewport(layout.pos.row = x,layout.pos.col = y)
print(p1, vp = vplayout(1:2,1))
print(p2, vp = vplayout(3,1))
dev.off()


