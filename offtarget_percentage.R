setwd("/mnt/Exome_Sequencing")
#install.packages("hrbrthemes")
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.csv("offtarget.csv")
data <- data.frame(data)
E_GRCh38 <- data$X..Off_target_GRCh38
E_CHM13 <- data$X.Off.target_CHM13
E_Decoy <- data$X..Off_target_GRCh38.Decoy

jpeg(file="offtarget.jpeg")
boxplot(E_CHM13, E_GRCh38, E_Decoy,
        main = "% Off-Target Bases",
        at = c(1,2,3),
        names = c("CHM13", "GRCh38", "GRCH38+Decoy"),
        las = 2,
        col = c("skyblue","seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE
)
dev.off()
