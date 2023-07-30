setwd("/mnt/Exome_Sequencing")
#install.packages("hrbrthemes")
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.csv("Boxplots_value.csv")
data <- data.frame(data)
E_CHM13 <- data$CHM13_IMP
E_Decoy <- data$DECoy_IMP

jpeg(file="Mapped_reads.jpeg")
boxplot(E_CHM13, E_Decoy,
        main = "Difference in mapped reads",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()


T_CHM13 <- data$CHM13IM.P
T_Decoy <- data$Decoy_IM.P

jpeg(file="MPR.jpeg")
boxplot(T_CHM13, T_Decoy,
        main = "Mapped and paired",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()

SNP_CHM13 <- data$CHM13_RUMR
SNP_Decoy <- data$GRCH38_RUMR
jpeg(file="UMR.jpeg")
boxplot(SNP_CHM13, SNP_Decoy,
        main = "Reduction in un mapped reads",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="Duplications.jpeg")
Ins_CHM13 <- data$MDR_CHM13
Ins_Decoy <- data$LDR_Decoy
boxplot(Ins_CHM13, Ins_Decoy,
        main = "Duplications",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="NPR.jpeg")
del_CHM13 <- data$I_NPA_CHM13

del_Decoy <- data$R_NPA_DECOY
boxplot(del_CHM13, del_Decoy,
        main = "Non_Primary Reads",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="PPR.jpeg")
indel_CHM13 <- data$CHM13_IPPR
indel_Decoy <- data$Decoy_RPPR
boxplot(indel_CHM13, indel_Decoy,
        main = "Properly paired reads",
        at = c(1,2),
        names = c("CHM13", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
jpeg(file="SAR.jpeg")
SAR_CHM13 <- data$CHM13_SAR
SAR_GRCH38 <- data$GRCh38_SAR
SAR_Decoy <- data$Decoy_SAR
boxplot(SAR_CHM13, SAR_GRCH38, SAR_Decoy,
        main = "SAR",
        at = c(1,2,3),
        names = c("CHM13", "GRCh38", "Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()