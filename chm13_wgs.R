library(vcfR)
setwd("/mnt/wgs/wgs_analysis_chm13")
vcf <- read.vcfR( "/mnt/wgs/wgs_analysis_chm13/Var_Recal/only_snp.vcf", verbose = FALSE )
dna <- ape::read.dna("/mnt/wgs/Reference/ref2/chm13v2.0.fa", format = "fasta")
gff <- read.table("/mnt/wgs/Reference/chm13.draft_v2.0.gene_annotation.gff3", sep="\t", quote="")
chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)
dev.off()
plot(chrom)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 100, max_DP = 200, min_MQ = 58,  max_MQ = 59)
plot(chrom)
chrom <- proc.chromR(chrom, verbose=TRUE)
jpeg("vcfR_plot1.jpg", width = 550, height = 550)
plot(chrom)
dev.off()
jpeg("vcfR_plot2.jpg", width = 550, height = 550)
chromoqc(chrom, dp.alpha=20)
dev.off()
chromoqc(chrom, xlim=c(5e+05, 6e+05))
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)

heatmap.bp(dp[1001:1500,])
is.na(dp[na.omit(dp == 0)]) <- TRUE
jpeg("heatmap.jpg", width = 550, height = 550)
heatmap.bp(dp[1001:1500,])
par(mar=c(8,4,4,2))
dev.off()
barplot(apply(dp, MARGIN=2, mean, na.rm=TRUE), las=3)
par(mar=c(5,4,4,2))
gt <- extract.gt(vcf)
gt
gt <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)
gt
gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE)
gt
gt <- extract.gt(vcf, element = 'HQ')
gt
myHQ1 <- masplit(gt[,1:2], sort = 0)
myHQ1
vcf_field_names(vcf, tag = "FORMAT")
Z <- vcfR2tidy(vcf, format_fields = c("GT", "DP"))
names(Z)
Z$meta
Z$fix
Z$gt
queryMETA(vcf, element = 'DP')
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", log='y', las=2)
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
abline(h=seq(0,1e4, by=100), col="#C0C0C088")
par(mar=c(5,4,4,2))
if( require(reshape2) & require(ggplot2) ){
  dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
  dpf <- dpf[ dpf$Depth > 0,]
  p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                         scale = "count", trim=TRUE)
  p <- p + theme_bw()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.text.x = element_text(angle = 60, hjust = 1, size=12))
  #  p <- p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
  p <- p + scale_y_continuous(trans=scales::log2_trans(), 
                              breaks=c(1, 10, 100, 800),
                              minor_breaks=c(1:10, 2:10*10, 2:8*100))
  p <- p + theme(axis.title.y = element_text(size=12))
  p <- p + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
  p <- p + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
  p <- p + stat_summary(fun.y=median, geom="point", shape=23, size=2)
  p
} else {
  message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed.  Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
}
sums <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[1,])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[2,])
dp[dp2 > 0] <- NA
dp[dp < 4] <- NA
par(mar=c(8,4,1,1))
jpeg("vcfR_plot3.jpg", width = 550, height = 550)
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth")
abline(h=seq(0,200, by=20), col="#C0C0C088")
par(mar=c(5,4,4,2))
dev.off()
is.na( vcf@gt[,-1][ is.na(dp) ] ) <- TRUE
vcf
