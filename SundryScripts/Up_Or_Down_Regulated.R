rm(list=ls())
cat("\014")
loadPackages <- function(){
  library(dplyr)
  library(cowplot)
  library(ggsignif)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(ineq)
  library(DESeq2)
  library(reshape2)
  library(GenomicFeatures)
  library(tidyverse)
  library(data.table)
  library(ggbiplot)
  library(RColorBrewer)
  library(VennDiagram)
  library(rtracklayer)
  library(ggridges)
  library(moments)
  library(e1071)
  library(gdata)
  library(readxl)
  library(GGally)
  library(corrplot)
  library(dendroextras)
  library(ggpubr)
  library(gplots)
  library(CliquePercolation) #version 0.3.0
  library(qgraph)            #version 1.6.5
  library(Matrix)            #version 1.2-18
  library(topGO)
  #library(rmarkdown)
}

loadPackages()
setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
data <- read_excel(paste0(getwd(),'/HeatMap.xlsx'))
colnames(data)<-make.names(names(data))
temp<-make.names((data$Biological.Processes))
data$Biological.Processes<-NULL
rownames(data)<-temp
data<-data.frame(data)
headers<-c('Ammonia', 'Lactate','Osmolality','Everything')
data1<-cbind(data$Ammonia_Up,data$Lactate_Up, data$Osmo_Up, data$Everything_Up)
colnames(data1)<-headers
rownames(data1)<-temp
data2<-cbind(data$Ammonia_Down,data$Lactate_Down, data$Osmo_Down, data$Everything_Down)
colnames(data2)<-headers
rownames(data2)<-temp
datatemp<-as.matrix(data2)
tiff("DownregulatedProcesses.tiff",bg = "transparent",width = 540, height = 480)
heatmap <- heatmap.2(datatemp, scale = 'none',
                     #add.expr=makeRects(selection),
                     labCol=colnames(datatemp),
                     labRow=gsub(".", " ", row.names(datatemp),fixed=TRUE),
                     #labCol= labelCol,
                     #labRow=labelRow,
                     Colv = FALSE, 
                     Rowv = FALSE,
                     dendrogram = 'none',
                     #ColSideColors = geneGrouping$coloring(colorsForGroups),
                     breaks=seq(0,0.5,0.005),
                     col = colorRampPalette(brewer.pal(8, "Blues"))(100),
                     density.info='none', denscol="white", 
                     symkey=FALSE, trace='none',
                     margins = c(6,10), cexCol=1, cexRow=1,
                     main="Downregulated Processes")
dev.off()

plot5<-ggplot(data,
              aes(x=Threshold, y=Time.to.Run))+
  geom_point(alpha=0.75,size=35)+
  geom_line(size=10)+
  xlim(0.6,1.0)+
  ylim(0,110)+
  theme_classic()+
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 75),
    axis.title.x = element_text(color = "black", size = 150, face = "bold"),
    axis.title.y = element_text(color = "black", size = 150, face = "bold"),
    axis.text.x = element_text(color = "black", size = 100, face = "bold"),
    axis.text.y = element_text(color = "black", size = 100, face = "bold")) +
  labs(x = expression(Threshold~Value), y = expression(Time~to~Run~Algorithm))
save_plot("Upregulated_Processes.jpeg", heatmap, ncol = 1, base_height=48, base_width=48,limitsize=FALSE)

