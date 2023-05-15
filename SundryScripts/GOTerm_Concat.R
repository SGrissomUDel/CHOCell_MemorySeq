rm(list=ls())
cat("\014")

loadPackages <- function(){
  library(dplyr)
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
  #library(rmarkdown)
}

loadPackages()

setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
#setwd('/Users/SpencerGrissom/Downloads')
getwd()

UnitProt_GO <- read_excel(paste0(getwd(),"/GO_Terms_UnitProt_CriGriPICR.xlsx"),guess_max = 20000)
DAVID_GO <- read_excel(paste0(getwd(),"/GO_Terms_DAVID_CriGriPICR.xlsx"),guess_max = 20000)
Molly_GO<- read_excel(paste0(getwd(),"/genes2GO_CHO_topGO.xlsx"),guess_max = 20000)

UnitProt_GO$geneID <- as.double(UnitProt_GO$geneID)
DAVID_GO$geneID <- as.double(DAVID_GO$geneID)
Molly_GO$geneID <- as.double(Molly_GO$geneID)

GO_Table <-left_join(DAVID_GO, UnitProt_GO,by='geneID')
GO_Table <-left_join(GO_Table, Molly_GO,by='geneID')



