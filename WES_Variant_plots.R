setwd("/mnt/Exome_Sequencing")
#install.packages("hrbrthemes")
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.csv("Difference_EXV.csv")
data <- data.frame(data)
E_GRCh38 <- data$GRCh38
E_CHM13 <- data$CHM13
E_Decoy <- data$Decoy

jpeg(file="Error.jpeg")
boxplot(E_CHM13, E_GRCh38, E_Decoy,
        main = "Error rate",
        at = c(1,2,3),
        names = c("CHM13", "GRCh38", "GRCH38+Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
data <- read.csv("TITVratiowes.csv")
data <- data.frame(data)
T_GRCh38 <- data$GRCh38
T_CHM13 <- data$CHM13
T_Decoy <- data$Decoy

jpeg(file="transition.jpeg")
boxplot(T_CHM13, T_GRCh38, T_Decoy,
        main = "Ti/Tv",
        at = c(1,2,3),
        names = c("CHM13", "GRCh38", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()

data <- data.frame(data)
SNP_CHM13 <- data$CHM13_DSNP
SNP_Decoy <- data$Decoy_DSNP
jpeg(file="SNPs.jpeg")
boxplot(SNP_CHM13, SNP_Decoy,
        main = "Reduction in SNPs",
        at = c(1,2),
        names = c("CHM13",  "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="Insertions.jpeg")
Ins_CHM13 <- data$CHM13_DINS
Ins_Decoy <- data$Decoy_DINS
boxplot(Ins_CHM13, Ins_Decoy,
        main = "Reduction in Insertions",
        at = c(1,2),
        names = c("CHM13",  "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="Deletions.jpeg")
del_CHM13 <- data$CHM13_DDEL
del_Decoy <- data$Decoy_DDEL
boxplot(del_CHM13, del_Decoy,
        main = "Reduction in Deletions",
        at = c(1,2),
        names = c("CHM13", "Decoy"),                                                                                                                                                                         
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="Indels.jpeg",width=600, height=600)
indel_CHM13 <- data$CHM13_DIND
indel_Decoy <- data$Decoy_DIND
boxplot(indel_CHM13, indel_Decoy,
         names = c("CHM13",  "GRCh38+Decoy"),
        main = "Reduction in indels",
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
