setwd("/mnt/Exome_Sequencing/Map_WXS_CHM13/var_recal")
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.table("Quality.txt")
data2 <- data[,c('V1','V2')]
jpeg(file="Qualityplot.jpeg")
boxplot(V1 ~ V2, data = data2,
        main = "Quality",
        at = c(1,2),
        names = c("Quality"),
        las = 2,
        border = "black",
        horizontal = TRUE,
        notch = FALSE
)
dev.off()

