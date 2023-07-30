setwd("/mnt/wgs")
library(readr)
library(RColorBrewer)
library(dplyr)
library(eulerr)
myCol <- brewer.pal(3, "Pastel2")
#grch38_snp <- read.table("/mnt/wgs/wgs_analysis_GRCh38/Var_Recal/only_snp_liftover_chm13.tsv")
#chm13_snp <- read.table("/mnt/wgs/wgs_analysis_chm13/Var_Recal/only_snp.tsv")

#grch38_chm13_snp <- intersect(grch38_snp,chm13_snp)
jpeg(file="snp_venn.jpeg")

VennDiag <- euler(c("GRCH38" = 1014585, "CHM13" = 919801, "GRCH38&CHM13" = 3566943 ), factor_names = TRUE)
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5, quantities = TRUE,
     fill=c("lightgreen", "lightblue"))

dev.off()
##############################################Venn diagram for indels#######################################
#grch38_indel <- read.table("/mnt/wgs/wgs_analysis_GRCh38/Var_Recal/only_indel_liftover_chm13.tsv")
#chm13_indel <- read.table("/mnt/wgs/wgs_analysis_chm13/Var_Recal/only_indel.tsv")

#grch38_chm13_indel <- intersect(grch38_indel,chm13_indel)
jpeg(file="indel_venn.jpeg")
VennDiag <- euler(c("GRCH38" = 282077, "CHM13" = 545680, "GRCH38&CHM13" = 437994 ), factor_names = TRUE)
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5, quantities = TRUE,
     fill=c("lightgreen", "lightblue"))

dev.off()
