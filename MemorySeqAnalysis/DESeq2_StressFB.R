rm(list=ls())
cat("\014")

######################################################################
######################################################################
##############        Functions and Such         #####################
######################################################################
######################################################################
#Functions include melting of data from count tables structured as 
#rows=gene count number, columns=replicate number, filtering functions,
#functions for summarizing data, and plotting gene expression

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
  library(hexbin)
  library(apeglm)
  #library(rmarkdown)
}
readDataIn <- function(ExcelFile,ColumnNames){
  loadPackages()
  setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
  #setwd('/Users/SpencerGrissom/Downloads')
  getwd()
  GTF <- readGFF("GCF_003668045.1_CriGri-PICR_genomic.gff")
  GTF_dat<-data.frame(GTF)
  GTF_dat$Dbxref <- as.vector(GTF$Dbxref)
  gene_id_to_biotype <- GTF_dat[-which(duplicated(GTF_dat$gene)),c('gene','gene_biotype')]
  
  data <- read_excel(paste0(getwd(),"/",ExcelFile))
  data <- data %>% filter(gene != "__alignment_not_unique" & gene != "__ambiguous" & gene != "__no_feature" & gene != "__not_aligned" & gene != "__too_low_aQual")
  colnames(data) <- ColumnNames
  
  data <-merge(data, gene_id_to_biotype, by='gene')
  data <- filter(data, gene_biotype=='protein_coding')
  data$gene_biotype<-NULL
  data2 <- data[,-1]
  rownames(data2) <- t(data[,1])
  return(data2)
}
loadPackages()


######################################################################
######################################################################
##############              Read Data In         #####################
######################################################################
######################################################################
#

######################################################################
ColumnNames <- c("gene","Control1","Control2","HighAmmonia1","HighAmmonia2","HighAmmonia3","HighLactate1","HighLactate2","HighLactate3","HighOsmolality1","HighOsmolality2","HighOsmolality3","HighEverything1","HighEverything2")
data <- readDataIn("Day5_CountTable.xlsx",ColumnNames)
######################################################################

dataInfo <- as.data.frame(matrix(nrow=(ncol(data)),ncol=2))
colnames(dataInfo) <- c("condition","type")
rownames(dataInfo) <- ColumnNames[-1]

######################################################################
dataInfo$condition <- c(rep("Control",2),rep("HighAmmonia",3),rep("HighLactate",3),rep("HighOsmolality",3),rep("HighEverything",2))
dataInfo$type <- rep("paired-end",13)
######################################################################

dataInfo$condition <-factor(dataInfo$condition)
dataInfo$type <-factor(dataInfo$type)


dds <- DESeqDataSetFromMatrix(countData = data, colData = dataInfo, design = ~ condition)
#keep <- rowSums(counts(dds)) > 1
#keep <- rowSums(counts(dds) >= 5) >= 3
#dds <- dds[keep,]


######################################################################
######################################################################
##############  Check the Transformation Methods #####################
######################################################################
######################################################################
#Transformation and visualization of the transformed data

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

dds1 <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds1, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


######################################################################
######################################################################
##############          PCA Analysis             #####################
######################################################################
######################################################################
#Generate PCA plot for the conditions tested

plotPCA(vsd, intgroup = c("condition"))

######################################################################
######################################################################
##############        DESeq Analysis             #####################
######################################################################
######################################################################
#Execute Differnential gene expression analysis and plot the spread of the data

dds<-DESeq(dds)
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()



######################################################################
######################################################################
#########  Isolate Results for Control Comparison ####################
######################################################################
######################################################################
#Store data based off the condition

resultsControl_HighAmmonia <- results(dds, contrast = c("condition", "HighAmmonia","Control"), alpha=0.05 )
resultsControl_HighLactate <- results(dds, contrast = c("condition" ,"HighLactate", "Control"), alpha=0.05  )
resultsControl_HighOsmolality <- results(dds, contrast = c("condition","HighOsmolality", "Control"), alpha=0.05  )
resultsControl_HighEverything <- results(dds, contrast = c("condition","HighEverything", "Control"), alpha=0.05  )


######################################################################
######################################################################
######  Find Significantly Differentially Expressed Genes ############
######################################################################
######################################################################
#Isolate genes that were significantly (p-value<0.05) differentially expressed
#(log base 2 fold change>0.58 or 1.5 times change in expression)
#Also transforms data to be just (Gene, L2FC, p-value,p-value adjusted)
#and isolates genes that were up or down-regulated



######################################################################
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
######################################################################

All_Amm <-resultsControl_HighLactate$log2FoldChange
All_Amm<-data.frame(cbind(All_Amm,resultsControl_HighLactate$pvalue,resultsControl_HighLactate$padj))
rownames(All_Amm) <- rownames(resultsControl_HighLactate)
colnames(All_Amm)<-c("log2FoldChange","pvalue",'padj')

All_Lac <-resultsControl_HighAmmonia$log2FoldChange
All_Lac<-data.frame(cbind(All_Lac,resultsControl_HighAmmonia$pvalue,resultsControl_HighAmmonia$padj))
rownames(All_Lac) <- rownames(resultsControl_HighAmmonia)
colnames(All_Lac)<-c("log2FoldChange","pvalue",'padj')

All_Osmo <-resultsControl_HighOsmolality$log2FoldChange
All_Osmo<-data.frame(cbind(All_Osmo,resultsControl_HighOsmolality$pvalue,resultsControl_HighOsmolality$padj))
rownames(All_Osmo) <- rownames(resultsControl_HighOsmolality)
colnames(All_Osmo)<-c("log2FoldChange","pvalue",'padj')

All_Every <-resultsControl_HighEverything$log2FoldChange
All_Every<-data.frame(cbind(All_Every,resultsControl_HighEverything$pvalue,resultsControl_HighEverything$padj))
rownames(All_Every) <- rownames(resultsControl_HighEverything)
colnames(All_Every)<-c("log2FoldChange","pvalue",'padj')

resultsControl_HighAmmonia_tb <- resultsControl_HighAmmonia %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

resultsControl_HighLactate_tb <- resultsControl_HighLactate %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

resultsControl_HighOsmolality_tb <- resultsControl_HighOsmolality %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

resultsControl_HighEverything_tb <- resultsControl_HighEverything %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_HighAmmonia <- resultsControl_HighAmmonia_tb %>%
  filter(pvalue < padj.cutoff & abs(log2FoldChange) > lfc.cutoff & !is.na(padj))

sig_HighLactate <- resultsControl_HighLactate_tb %>%
  filter(pvalue < padj.cutoff & abs(log2FoldChange) > lfc.cutoff & !is.na(padj))

sig_HighOsmolality <- resultsControl_HighOsmolality_tb %>%
  filter(pvalue < padj.cutoff & abs(log2FoldChange) > lfc.cutoff & !is.na(padj))

sig_HighEverything <- resultsControl_HighEverything_tb %>%
  filter(pvalue < padj.cutoff & abs(log2FoldChange) > lfc.cutoff & !is.na(padj))

HighAmmoniaSet <- sig_HighAmmonia$gene
HighLactateSet <- sig_HighLactate$gene
HighOsmolalitySet <- sig_HighOsmolality$gene
HighEverythingSet <- sig_HighEverything$gene

HighAmmoniaSetNeg <- filter(sig_HighAmmonia,log2FoldChange<0)$gene
HighLactateSetNeg <- filter(sig_HighLactate,log2FoldChange<0)$gene
HighOsmolalitySetNeg <- filter(sig_HighOsmolality,log2FoldChange<0)$gene
HighEverythingSetNeg <- filter(sig_HighEverything,log2FoldChange<0)$gene

HighAmmoniaSetPos <- filter(sig_HighAmmonia,log2FoldChange>0)$gene
HighLactateSetPos <- filter(sig_HighLactate,log2FoldChange>0)$gene
HighOsmolalitySetPos <- filter(sig_HighOsmolality,log2FoldChange>0)$gene
HighEverythingSetPos <- filter(sig_HighEverything,log2FoldChange>0)$gene


######################################################################
######################################################################
##############  Venn Diagram for Overlap         #####################
######################################################################
######################################################################
#Plots the overlap of differentially expressed genes and the ones that were up
#or down-regulated

v2 <- venn.diagram(list(HighAmmonia=HighAmmoniaSet, HighLactate=HighLactateSet, HighOsmolality=HighOsmolalitySet, HighEverything=HighEverythingSet),
                   fill = c("red", "green", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("DifferentialExpression_FB.jpeg")
grid.newpage()
grid.draw(v2)
dev.off()

v2 <- venn.diagram(list(HighAmmonia=HighAmmoniaSetNeg, HighLactate=HighLactateSetNeg, HighOsmolality=HighOsmolalitySetNeg, HighEverything=HighEverythingSetNeg),
                   fill = c("red", "green", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("DifferentialExpression_Neg_FB.jpeg")
grid.newpage()
grid.draw(v2)
dev.off()

v2 <- venn.diagram(list(HighAmmonia=HighAmmoniaSetPos, HighLactate=HighLactateSetPos, HighOsmolality=HighOsmolalitySetPos, HighEverything=HighEverythingSetPos),
                   fill = c("red", "green", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("DifferentialExpression_Pos_FB.jpeg")
grid.newpage()
grid.draw(v2)
dev.off()
######################################################################
######################################################################
#########  Heritable Genes and Differentially Expressed ##############
######################################################################
######################################################################
#Using a list of the heritable genes, also look at the overlap of differntially
#expressed genes and those with heritable patterns

######################################################################
HeritableList_MemorySeq<-read.table("HeritableGenesCHO.txt")
######################################################################

colnames(HeritableList_MemorySeq)<-"gene"
commonHighAmmonia <- intersect(sig_HighAmmonia$gene, HeritableList_MemorySeq$gene)
commonHighLactate <- intersect(sig_HighLactate$gene, HeritableList_MemorySeq$gene)
commonHighOsmolality <- intersect(sig_HighOsmolality$gene, HeritableList_MemorySeq$gene)
commonHighEverything <- intersect(sig_HighEverything$gene, HeritableList_MemorySeq$gene)

temp1<-data.frame(commonHighAmmonia)
colnames(temp1)<-'gene'
temp2<-data.frame(commonHighLactate)
colnames(temp2)<-'gene'
temp3<-data.frame(commonHighOsmolality)
colnames(temp3)<-'gene'
temp4<-data.frame(commonHighEverything)
colnames(temp4)<-'gene'
temp<-rbind(temp1,temp2,temp3,temp4)
UniqueDEG<-unique(temp)

v2 <- venn.diagram(list(HighAmmonia=commonHighAmmonia, HighLactate=commonHighLactate, HighOsmolality=commonHighOsmolality, HighEverything=commonHighEverything),
                   fill = c("red", "green", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("DifferentialExpression_FB_MemorySeq.jpeg")
grid.newpage()
grid.draw(v2)
dev.off()




######################################################################
######################################################################
#######  Generate Heat Maps Based on Z-Score and LFC #################
######################################################################
######################################################################
#Generate whole transcriptome heat maps

log2Normalization <- function(x){
  return(log(x,2)+1)
}
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

HeritableList<-read.table("HeritableGenesCHO.txt")
colnames(HeritableList)<-"gene"
RawCount<-data.frame(counts(dds,normalized=FALSE))
RawCount_Heritable<-filter(RawCount, row.names(RawCount) %in% HeritableList$gene)
NormalizedCount<-data.frame(counts(dds,normalized=TRUE))
NormalizedCount_Heritable<-filter(NormalizedCount, row.names(NormalizedCount) %in% HeritableList$gene)


log2_NormalizedCount_Heritable<-data.frame(lapply(NormalizedCount_Heritable,log2Normalization))
ZScore_Heritable <- t(apply(RawCount_Heritable, 1, cal_z_score))

heatmap <- heatmap.2(ZScore_Heritable, scale = 'none',
                     #add.expr=makeRects(selection),
                     labCol=colnames(ZScore_Heritable),
                     #labRow=row.names(ZScore_Heritable),
                     #labCol= labelCol,
                     #labRow=labelRow,
                     #Colv = as.dendrogram(geneClustering_expSample), 
                     #Rowv = as.dendrogram(geneClustering_expSample),
                     dendrogram = 'none',
                     #ColSideColors = geneGrouping$coloring(colorsForGroups),
                     col = redgreen(75),
                     density.info='density', denscol="white", 
                     symkey=FALSE, trace='none',
                     margins = c(6,10), cexCol=1, cexRow=1,
                     Colv=FALSE,
                     main="Row Z-Score For Heritable Genes")


temp_resultsControl_HighAmmonia_tb<-resultsControl_HighAmmonia_tb[complete.cases(resultsControl_HighAmmonia_tb), ]
temp_resultsControl_HighLactate_tb<-resultsControl_HighLactate_tb[complete.cases(resultsControl_HighLactate_tb), ]
temp_resultsControl_HighOsmolality_tb<-resultsControl_HighOsmolality_tb[complete.cases(resultsControl_HighOsmolality_tb), ]
temp_resultsControl_HighEverything_tb<-resultsControl_HighEverything_tb[complete.cases(resultsControl_HighEverything_tb), ]

AmmoniaLFC<-cbind(temp_resultsControl_HighAmmonia_tb$gene,temp_resultsControl_HighAmmonia_tb$log2FoldChange)
AmmoniaLFC<-data.frame(AmmoniaLFC)
colnames(AmmoniaLFC)<-c('gene','LogFold2Change')
rownames(AmmoniaLFC)<-AmmoniaLFC$gene
AmmoniaLFC$gene<-NULL
AmmoniaLFC<-AmmoniaLFC[HeritableList$gene, ]

LactateLFC<-cbind(temp_resultsControl_HighLactate_tb$gene,temp_resultsControl_HighLactate_tb$log2FoldChange)
LactateLFC<-data.frame(LactateLFC)
colnames(LactateLFC)<-c('gene','LogFold2Change')
rownames(LactateLFC)<-LactateLFC$gene
LactateLFC$gene<-NULL
LactateLFC<-LactateLFC[HeritableList$gene, ]

OsmolalityLFC<-cbind(temp_resultsControl_HighOsmolality_tb$gene,temp_resultsControl_HighOsmolality_tb$log2FoldChange)
OsmolalityLFC<-data.frame(OsmolalityLFC)
colnames(OsmolalityLFC)<-c('gene','LogFold2Change')
rownames(OsmolalityLFC)<-OsmolalityLFC$gene
OsmolalityLFC$gene<-NULL
OsmolalityLFC<-OsmolalityLFC[HeritableList$gene, ]

EverythingLFC<-cbind(temp_resultsControl_HighEverything_tb$gene,temp_resultsControl_HighEverything_tb$log2FoldChange)
EverythingLFC<-data.frame(EverythingLFC)
colnames(EverythingLFC)<-c('gene','LogFold2Change')
rownames(EverythingLFC)<-EverythingLFC$gene
EverythingLFC$gene<-NULL
EverythingLFC<-EverythingLFC[HeritableList$gene, ]

LFC_matrix<-cbind(AmmoniaLFC,LactateLFC,OsmolalityLFC,EverythingLFC)
LFC_matrix <-data.frame(LFC_matrix)
LFC_matrix <- sapply(LFC_matrix, as.numeric )
rownames(LFC_matrix)<-HeritableList$gene
colnames(LFC_matrix)<-c('HighAmmonia','HighLactate','HighOsmolality','HighEverything')


heatmap <- heatmap.2(LFC_matrix, scale = 'none',
                     #add.expr=makeRects(selection),
                     labCol=colnames(LFC_matrix),
                     #labRow=row.names(ZScore_Heritable),
                     dendrogram = 'row',
                     hclustfun=hclust,
                     Colv = FALSE, 
                     #Rowv = as.dendrogram(geneClustering_expSample_c),
                     col = redgreen(75),
                     density.info='density', denscol="white", 
                     symkey=FALSE, trace='none',
                     margins = c(6,10), cexCol=1, cexRow=1,
                     main="Log2Fold Change For Heritable Genes")



