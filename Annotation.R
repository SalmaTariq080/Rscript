#BiocManager::install("maftools",force = TRUE)
#library(maftools)
#BiocManager::install("VariantAnnotation")
#install.packages("ggplot2")
library(VariantAnnotation)
library(ggplot2)
setwd("/mnt/Exome_Sequencing/annovar")
#vcf <- readVcf("/home/salmatariq/exome_sequencing/GenotypeGVCFs/merged.vcf.gz", "hg38")
#vcf <- readVcf("/mnt/Exome_Sequencing/annovar/myanno.hg38_multianno.vcf", "hg38")
vcf <- readVcf("/home/salmatariq/exome_sequencing/Filtered_indel.vcf.gz", "hg38")
rowRanges(vcf)
info(vcf)
header(vcf)
sampleID <- samples(header(vcf))
sampleID
meta(header(vcf))
meta(header(vcf))$fileformat
meta(header(vcf))$source
meta(header(vcf))$contig
rd <- rowRanges(vcf)
rd[1:2]
as.vector(seqnames(rd)[1:2])
start(rd)[1:2]
end(rd)[1:2]
refBase <- ref(vcf)
refBase[1:2]
refBase <- as.character(refBase)
refBase[1:2]
altBase <- alt(vcf)
alt(vcf)[1:2]
altBase <- unlist(lapply(altBase, as.character))
###########################################################
info(header(vcf))[1:2, ]
info(vcf)[1:2, ]
geno(header(vcf))[1:2, ]
paste0("GT: ", geno(header(vcf))[1, 3])
matGT <- geno(vcf)$GT
matGT[1:2, ]
tbl <- table(geno(vcf)$GT)
tbl_dat <- as.data.frame(tbl)
#pdf("exome_sequencing.pdf")
jpeg(file="variant_frequency_plot_indel_af.jpeg")
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme_classic()
dev.off()
paste0("DP: ", geno(header(vcf))[3, 3])
matDP <- geno(vcf)$DP
matDP[1:2, ]
summary(as.vector(matDP))
jpeg(file="depth_plot_indel_af.jpeg")
ggplot(as.data.frame(matDP),aes(x=Sample1))+geom_histogram()+
  labs(x="",y="Counts")+
  scale_x_log10()+
  theme_classic()
dev.off()
paste0("GQ: ", geno(header(vcf))[4, 3])
matGQ <- geno(vcf)$GQ
matGQ[1:2, ]
summary(as.vector(matGQ))
jpeg(file="Genotype_Quality_plot_indel_af.jpeg")
ggplot(as.data.frame(matGQ),aes(x=Sample1))+geom_histogram()+
  labs(x="",y="Counts")+
  scale_x_log10()+
  theme_classic()
dev.off()
var_1 <- rownames(geno(vcf)$GT)[
  geno(vcf)$GT=="0/1" | 
    geno(vcf)$GT=="1/1"]
var_1
varTab1 <- data.frame(variant=names(rd)[names(rd) %in% var_1],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_1]),
                      start=start(rd)[names(rd) %in% var_1],
                      end=end(rd)[names(rd) %in% var_1],
                      stringsAsFactors = FALSE)
ref_base <- ref(vcf)[rownames(vcf) %in% var_1]
ref_base[1:2]
varTab1$refBase <- as.character(ref_base)
##############################################################################
alt_base <- lapply(alt(vcf)[rownames(vcf) %in% var_1],`[[`,1)
alt_base[1]
alt_base <- lapply(alt_base,as.character)
alt_base[1]
varTab1$altBase <- unlist(alt_base)
adCount <- geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_1]
adCount[1]
varTab1$refCount <- unlist(lapply(adCount,`[[`,1))
varTab1$altCount <- unlist(lapply(adCount,`[[`,2))
varTab1$genoType <- geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_1]
varTab1$gtQuality <- geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_1]
var_2 <- rownames(geno(vcf)$GT)[geno(vcf)$GT=="1/2"]
varTab2 <- data.frame(variant=names(rd)[names(rd) %in% var_2],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_2]),
                      start=start(rd)[names(rd) %in% var_2],
                      end=end(rd)[names(rd) %in% var_2],
                      refBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,1),as.character)),
                      altBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,2),as.character)),
                      refCount=unlist(lapply(
                        geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,2)),
                      altCount=unlist(
                        lapply(geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,3)),
                      genoType=geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_2],
                      gtQuality=geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_2],
                      stringsAsFactors = FALSE)
varTab <- rbind(varTab1, varTab2)
varTab[1:2, ]
for (k in 1:length(varTab$variant)) {
  if (width(varTab$refBase[k]) < width(varTab$altBase[k])) {
    varTab$mutType[k] <- "INS"
  } else if (width(varTab$refBase[k]) > width(varTab$altBase[k])) {
    varTab$mutType[k] <- "DEL"
  } else if (width(varTab$refBase[k]) == 1 & width(varTab$altBase[k]) == 1) {
    varTab$mutType[k] <- "SNP"
  } else {
    varTab$mutType[k] <- "Others"
  }
}
########################################################################################
for (k in 1:length(varTab1$variant)) {
  if (width(varTab1$refBase[k]) < width(varTab1$altBase[k])) {
    varTab1$mutType[k] <- "INS"
  } else if (width(varTab1$refBase[k]) > width(varTab1$altBase[k])) {
    varTab1$mutType[k] <- "DEL"
  } else if (width(varTab1$refBase[k]) == 1 & width(varTab1$altBase[k]) == 1) {
    varTab1$mutType[k] <- "SNP"
  } else {
    varTab1$mutType[k] <- "Others"
  }
}
tbl <- table(varTab1$mutType)
tbl_dat <- as.data.frame(tbl)
tbl
jpeg(file="Variations_summary_indel_af.jpeg")
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+
  theme_classic()
dev.off()

ti <- c("A>G","G>A","C>T","T>C")
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
varTab1$nuSub <- paste0(varTab1$refBase,">",varTab1$altBase)
varTab1$TiTv[varTab1$nuSub %in% ti] <- "Ti"
varTab1$TiTv[varTab1$nuSub %in% tv] <- "Tv"
varTab1[1:2,]
varX <- varTab1[varTab1$mutType == "SNP", ]
tbl <- table(varX$nuSub)
tbl_dat <- as.data.frame(tbl)
tbl
jpeg(file="transitions_summary_indel_af.jpeg")
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+
  theme(legend.position = "none")
dev.off()
tbl <- table(varX$TiTv)
tbl_dat <- as.data.frame(tbl)
tbl
jpeg(file="Ti_Tv_indel_af.jpeg")
ggplot(as.data.frame(table(varX$TiTv)),aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
  theme(legend.position = "none")
dev.off()
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(GenomicFeatures)
#library(stringr)
########################################################################################

#dev.off()
