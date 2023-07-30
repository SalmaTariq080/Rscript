library(VennDiagram)
setwd("/mnt/Exome_Sequencing")
library(readr)
library(RColorBrewer)
library(dplyr)
library(eulerr)

myCol <- brewer.pal(3, "Pastel2")
#grch38_snp <- read.table("/mnt/Exome_Sequencing/Map_WXS_GRCH38/var_recal/only_snp_liftover_chm13.tsv")
#decoy_snp <- read.table("/mnt/Exome_Sequencing/Map_WXS_Decoy/var_recal/only_snp_liftover_chm13.tsv")
#chm13_snp <- read.table("/mnt/Exome_Sequencing/Map_WXS_CHM13/var_recal/only_snp.tsv")
#grch38_decoy_snp <- intersect(grch38_snp,decoy_snp)
#decoy_chm13_snp<-intersect(decoy_snp,chm13_snp)
#grch38_chm13_snp<-intersect(grch38_snp,chm13_snp)
#grch38_decoy_c_hm13snp<- intersect(intersect(grch38_snp,decoy_snp),chm13_snp)
jpeg(file="snp_venn.jpeg")

VennDiag <- euler(c("GRCH38" = 52798, "Decoy" = 18025, "CHM13" = 490352, "GRCH38&Decoy" = 479853, "Decoy&CHM13" = 1929, 
                    "GRCH38&CHM13" = 47191, "GRCH38&Decoy&CHM13" = 923508), factor_names = TRUE)
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5, quantities = TRUE,
     fill=c("lightblue", "lightpink", "lightgreen"))


dev.off()
##############################################Venn diagram for indels#######################################
#grch38_indel <- read.table("/mnt/Exome_Sequencing/Map_WXS_GRCH38/var_recal/only_indel_liftover_chm13.tsv")
#decoy_indel <- read.table("/mnt/Exome_Sequencing/Map_WXS_Decoy/var_recal/only_indel_liftover_chm13.tsv")
#chm13_indel <- read.table("/mnt/Exome_Sequencing/Map_WXS_CHM13/var_recal/only_indel.tsv")

#grch38_decoy_indel <- intersect(grch38_indel,decoy_indel)
#decoy_chm13_indel<-intersect(decoy_indel,chm13_indel)
##grch38_chm13_indel<-intersect(grch38_indel,chm13_indel)
#grch38_decoy_chm13_indel<-intersect(intersect(grch38_indel,decoy_indel),chm13_indel)
jpeg(file="indel_venn.jpeg")
VennDiag <- euler(c("GRCH38" = 69666, "Decoy" = 57706, "CHM13" = 109481, "GRCH38&Decoy" = 45100, "Decoy&CHM13" = 485, 
                    "GRCH38&CHM13" = 13988, "GRCH38&Decoy&CHM13" = 160299), factor_names = TRUE)
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5, quantities = TRUE,
     fill=c("lightblue", "lightpink", "lightgreen"))
dev.off()

#############################################another way to make venn diagrams###################################################################
draw.triple.venn(911332, 853737, 1462981, 
                 844394, 819356, 859721, 817929, 
                 category=c("GRCh38","Decoy","CHM13"),
                 col=myCol,
                 cex = rep(1, 7), fontface = rep("bold", 7), fontfamily = rep("serif", 7), sep.dist = 0.1, rotation.degree = 80, euler.d = TRUE, scaled = TRUE)
dev.off()
