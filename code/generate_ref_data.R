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

# xReactiveProbes object
xReactiveProbesDir <- "/home/cbw/projects/ewas/data/ref/48639-non-specific-probes-Illumina450k.xlsx"
xReactiveProbes <- read.xlsx(xReactiveProbesDir,sheet = "nonspecific cg probes")

# iccref object
iccrefDir <- "/home/cbw/projects/ewas/data/ref/Supplementary_Table_5.xlsx"
iccref <- read.xlsx(iccrefDir)

# ann450k object
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## save objects =====================
save(xReactiveProbes, ann450k, iccref, file = "/home/bailing/projects/ewas/analysis/reffile/ref.RData")
