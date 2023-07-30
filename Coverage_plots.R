setwd("/mnt/Exome_Sequencing")
#install.packages("hrbrthemes")
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.csv("WES_Coverage.csv")
data <- data.frame(data)
a = data$X31
jpeg(file="Coverage.jpeg")
boxplot(a, 
        main = "Coverage",
        at = c(1),
        names = c("DEL"),
        las = 2,
        col = c("seagreen"),
        border = "black",
        horizontal = FALSE,
        notch = FALSE,
        plot = TRUE
)
dev.off()
