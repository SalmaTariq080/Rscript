remotes::install_github("UofABioinformaticsHub/shinyNgsreports")
library(ngsReports)
install.packages('remotes')
setwd("/mnt/RNAseq-samples/assemble_transcript/trim_galore_results/fastqc_after_TrGl")
files <- list.files(getwd(), pattern = "fastqc.zip", full.names = TRUE)
fdl <- FastqcDataList(files)
getModule(fdl[[1]], "Summary")
reads <- readTotals(fdl)
library(dplyr)
library(pander)
reads %>%
  dplyr::filter(grepl("R1", Filename)) %>% 
  pander(
    big.mark = ",",
    caption = "Read totals from R1 libraries", 
    justify = "lr"
  )
#plot summary for all samples 
jpeg("Summary.jpg", width = 550, height = 550)
plotSummary(fdl)
dev.off()
#Visulizing total reads
jpeg("reads_total.jpg", width = 550, height = 550)
plotReadTotals(fdl)
dev.off()
plotReadTotals(fdl) +
  theme(
    legend.position = c(1, 1), 
    legend.justification = c(1, 1),
    legend.background = element_rect(colour = "black")
  )
#Plotting per base sequence quality
#plotBaseQuals(fdl[[1]])
jpeg("Base_quality.jpg", width = 550, height = 550)
plotBaseQuals(fdl)
dev.off()
#plotBaseQuals(fdl[1:4], plotType = "boxplot")
#Plot mean sequence qulaities
jpeg("Seq_quality.jpg", width = 550, height = 550)
plotSeqQuals(fdl)
dev.off()
#Plot per base sequence content
jpeg("Seq_content.jpg", width = 550, height = 550)
plotSeqContent(fdl)
dev.off()
#Plotting adapter content
jpeg("Adapter_content.jpg", width = 550, height = 550)
plotAdapterContent(fdl)
dev.off()
#Plotting sequence duplication levels
jpeg("Duplication_levels.jpg", width = 550, height = 550)
plotDupLevels(fdl)
dev.off()
gcAvail(gcTheoretical, "Genome")
#Plotting GC Content
jpeg("GC_content.jpg", width = 550, height = 550)
plotGcContent(fdl)
dev.off()
#Plotting overrepresented Sequence
jpeg("over_rep.jpg", width = 550, height = 550)
plotOverrep(fdl)
dev.off()

