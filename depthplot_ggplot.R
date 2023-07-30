setwd("/mnt/wgs")
library(ggplot2)
df<-data.frame(mean=c(37.729,48.4094),sd=c(44.4629,50.4249),Category=as.factor(c("ERR1395554","ERR1395557")))
jpeg(file="coverage_plot1.jpeg")
ggplot(df,aes(x=Category))+geom_boxplot(aes(lower=mean-sd,upper=mean+sd,middle=mean,ymin=mean-3*sd,ymax=mean+3*sd),stat="identity") + coord_cartesian( ylim = c(-10, 100))
dev.off()
?ggplot
