rm(list=ls())
cat("\014")

######################################################################
######################################################################
##############        Functions and Such         #####################
######################################################################
######################################################################
print("Loading all Packages Now")  
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
  library(chromoMap)
  library(Rtsne)
}
meltData <- function(ExcelFileNoise, ExcelFileMemory){
  loadPackages()
  setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
  #setwd('/Users/SpencerGrissom/Downloads')
  getwd()
  #GTF <- readGFF(paste0(getwd(),"/Cricetulus_griseus_picr.CriGri-PICR.104.gtf"))
  GTF <- readGFF(paste0(getwd(),"/GCF_003668045.1_CriGri-PICR_genomic.gff"))
  #GTF <- readGFF(paste0(getwd(),"/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff"))
  GTF_dat<-data.frame(GTF)
  GTF_dat$Dbxref <- as.vector(GTF$Dbxref)
  gene_id_to_biotype <- GTF_dat[-which(duplicated(GTF_dat$gene)),c('gene','gene_biotype','Dbxref')]
  dataNoise <- read_excel(paste0(getwd(),"/",ExcelFileNoise))
  dataNoise <- dataNoise %>% filter(gene != "__alignment_not_unique" & gene != "__ambiguous" & gene != "__no_feature" & gene != "__not_aligned" & gene != "__too_low_aQual")
  
  SampleSize <- dataNoise %>% dplyr::select(starts_with("Bulk")) %>% ncol(.)
  GeneCount <- dataNoise %>% dplyr::select(starts_with("Bulk")) %>% nrow(.)
  collapseNoise <- dataNoise %>% dplyr::select(starts_with("Bulk")) %>% data.matrix(.) %>% t(.) %>% as.vector(.)
  collapseNoise <- cbind("gene"=rep(dataNoise$gene,each=SampleSize),"controls"=rep("Controls",each=SampleSize*GeneCount),"number"=rep(1:SampleSize,times=GeneCount),"total_counts"=collapseNoise)
  
  dataMemory <- read_excel(paste0(getwd(),"/",ExcelFileMemory))
  dataMemory <- dataMemory %>% filter(gene != "__alignment_not_unique" & gene != "__ambiguous" & gene != "__no_feature" & gene != "__not_aligned" & gene != "__too_low_aQual")
  
  SampleSize <- dataMemory %>% dplyr::select(starts_with("Single")) %>% ncol(.)
  GeneCount <- dataMemory %>% dplyr::select(starts_with("Single")) %>% nrow(.)
  collapseMemory <- dataMemory %>% dplyr::select(starts_with("Single")) %>% data.matrix(.) %>% t(.) %>% as.vector(.)
  collapseMemory <- cbind("gene"=rep(dataMemory$gene,each=SampleSize),"controls"=rep("ExpSample",each=SampleSize*GeneCount),"number"=rep(1:SampleSize,times=GeneCount),"total_counts"=collapseMemory)
  
  collapseData <- rbind(collapseNoise,collapseMemory)
  collapseData <- data.frame(collapseData)
  collapseData <- arrange(collapseData,gene,controls,number)
  
  collapseData <-merge(collapseData, gene_id_to_biotype, by='gene')
  collapseData <- filter(collapseData, gene_biotype=='protein_coding')
  collapseData <- collapseData %>% separate(Dbxref,into=c("TBR","geneID"), remove=FALSE, sep=":")
  #txdb <- makeTxDbFromGFF(file = paste0(getwd(),"/Cricetulus_griseus_picr.CriGri-PICR.104.gtf"), format="gtf")
  txdb <- makeTxDbFromGFF(file = paste0(getwd(),"/GCF_003668045.1_CriGri-PICR_genomic.gff"))
  lengthsPergeneid <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene"))))
  lengthtbl<-data.frame(gene = names(lengthsPergeneid), length = lengthsPergeneid)
  collapseData <- left_join(collapseData, lengthtbl, by='gene')
  collapseData$total_counts <- as.integer(collapseData$total_counts)
  
  collapseData<-collapseData %>%
    dplyr::group_by(controls, number) %>%
    dplyr::mutate(totalMappedReads = sum(total_counts), 
                  rpm = 1000000*total_counts/totalMappedReads, 
                  rpk = 1000*total_counts/length,
                  rpkScalePerMillion = sum(rpk)/1000000,
                  tpm = rpk/rpkScalePerMillion,
                  rpkm = 1000*rpm/length) 
  collapseData$gene_version <- NULL
  collapseData$gene_source <- NULL
  collapseData$gene_biotype <- NULL
  collapseData$TBR <- NULL
  collapseData$Dbxref <- NULL
  GO_ID <- read_excel(paste0(getwd(),"/Final_GO_ID_CriGriPICR.xlsx"))
  collapseData$geneID <- as.double(collapseData$geneID)
  collapseData <- left_join(collapseData, GO_ID,by='geneID')
  
  #dataNoise_rpm <- dataNoise %>%
  #dplyr::transmute_at(.vars=vars(starts_with("Bulk")), .funs=funs(1000000*./sum(.)))
  #dataNoise_rpm <- cbind("gene_id"=dataNoise$gene_id, dataNoise_rpm, "length"=dataNoise$length)
  #dataNoise %>% dplyr::select(starts_with("Bulk")) %>% apply(.,2,sum)
  
  return(collapseData)
}
remove_genes_less1_5tpm <- function(data){
  data <- dplyr::filter(data, totalMappedReads>0.5*10^6)
  data2 <- data.table(data)
  geneExpr <- data2[, list(medianTPM = median.default(tpm), length=length), 
                    by = c('gene')]
  geneExpr <- transform(geneExpr, 
                        keepGene = (medianTPM > 2.5 & length>1000) | 
                          (medianTPM > 3 & length>1000) | 
                          (medianTPM > 10 & length>500) |
                          (medianTPM > 50 & length>100))
  #data2 <- merge(data2, subset(geneExpr, keepGene, select=c('gene_id')),by='gene_id')
  
  keepGenes <- filter(geneExpr, keepGene=='TRUE')
  data2 <- left_join(unique(keepGenes), data2, by=c('gene','length'))
  data2$keepGene <- NULL
  
  return(data2)
}
remove_genes_less1_rpm <- function(data){
  data <- collapseData
  data <- dplyr::filter(data, totalMappedReads>0.5*10^6)
  data2 <- data.table(data)
  geneExpr <- data2[, list(meanRPM = mean(rpm)), 
                    by = c('gene')]
  geneExpr <- transform(geneExpr, 
                        keepGene = meanRPM > 1)
  #data2 <- merge(data2, subset(geneExpr, keepGene, select=c('gene')),by='gene')
  
  keepGenes <- filter(geneExpr, keepGene=='TRUE')
  data2 <- left_join(unique(keepGenes), data2, by=c('gene'))
  data2$keepGene <- NULL
  return(data2)
}
plotRPMsForGene <- function(geneID,data){
  dataSub <- filter(data, gene == geneID)
  temp <- ks.test(x=filter(dataSub,controls=='ExpSample')$rpm, 
                  y=filter(dataSub,controls=='Controls')$rpm)
  g <- ggplot(dataSub, aes(x=controls, y=rpm, color=controls))+
    geom_boxplot(outlier.alpha = 0)+
    geom_jitter()+
    theme_classic()+
    geom_signif(comparisons=list(c("Controls","ExpSample")),map_signif_level = TRUE)+
    ggtitle(geneID)
  return(g)
}
calc_variability_metrics <- function(data){
  # Filter out samples that have less than 1 million reads ---------------------------
  data <- filter(data, totalMappedReads>0.5*10^6)
  
  # Calculate CVs, fanos, etc. ---------------------------
  data_with_metrics <- data %>% dplyr::group_by(gene,controls,geneID) %>% 
    dplyr::summarise(sd(rpm) / mean(rpm), # CV from rpm
                     sd(tpm) / mean(tpm), # CV from tpm
                     mean(rpm), # mean rpm
                     mean(tpm), # mean tpm) 
                     moments::skewness(tpm),
                     moments::kurtosis(tpm),
                     e1071::kurtosis(tpm,,3),
                     sd(log2(tpm)) / mean(log2(tpm)),
                     mean(log2(tpm)),
                     sd(log2(tpm+1)) / mean(log2(tpm+1)),
                     mean(log2(tpm+1))
    )
  colnames(data_with_metrics) <- c('gene','controls','geneID',
                                   'CV_rpm','CV_tpm','Mean_rpm','Mean_tpm', 'skewness','kurtosis',
                                   'kurtosis_e1071', 'CV_log2tpm', 'Mean_log2tpm','CV_log2_tpm_plus1','Mean_log2_tpm_plus1')
  return(data_with_metrics)
}

######################################################################
######################################################################
##############        Collect Data and Summarize #####################
######################################################################
######################################################################
print("Loading data in Now")  

#collapseDataFull <- meltData("Noise_Control_FinalTable.xlsx","Memory_Seq_Samples_FinalTable.xlsx")
collapseData <- meltData("Noise_Control_FinalTable.xlsx","Memory_Seq_Samples_FinalTable_Reduced.xlsx")
#collapseData2 <- meltData("Noise_Control_FinalTable.xlsx","Memory_Seq_Samples_FinalTable_Reduced2.xlsx")

#collapseData<-collapseData2
TrimmedData <- remove_genes_less1_5tpm(collapseData)

rm(collapseData)
#plotRPMsForGene("Aamdc",TrimmedData)
summary <- calc_variability_metrics(TrimmedData)


######################################################################
######################################################################
#########Poisson Regression Fitting for Heritability##################
######################################################################
######################################################################
exp_sample_data <- filter(summary, controls=='ExpSample')
control_data <- filter(summary, controls=='Controls')
percentile_cutoff <- 0.98
mean_tpm_cutoff <- 2.5

model_exp_sample <- glm(CV_tpm ~ log2(Mean_tpm), family = poisson(link="log"), 
                        data = exp_sample_data)
model_controls <- glm(CV_tpm ~ log2(Mean_tpm), family = poisson(link="log"), 
                      data = control_data)

fit_data_exp <- data.frame(gene = exp_sample_data$gene, 
                           #GeneSymbol = exp_sample_data$GeneSymbol,
                           controls = exp_sample_data$controls,
                           Mean_tpm = exp_sample_data$Mean_tpm, 
                           CV_tpm = exp_sample_data$CV_tpm,
                           fitted_y= fitted(model_exp_sample),
                           resid  = resid(model_exp_sample))

fit_data_controls <- data.frame(gene = control_data$gene, 
                                #GeneSymbol = control_data$GeneSymbol,
                                controls = control_data$controls,
                                Mean_tpm = control_data$Mean_tpm,
                                CV_tpm = control_data$CV_tpm,
                                fitted_y= fitted(model_controls),
                                resid = resid(model_controls))

fit_data <- rbind(fit_data_exp, fit_data_controls)
rm(fit_data_controls)
rm(fit_data_exp)

cutoff099 <- quantile(filter(fit_data, controls=='ExpSample' & Mean_tpm>mean_tpm_cutoff)$resid, percentile_cutoff)

fit_data$cutoff <- factor((fit_data$Mean_tpm > mean_tpm_cutoff) & (fit_data$resid > cutoff099))

HeritableList_MemorySeq<-filter(fit_data,cutoff==TRUE & controls=='ExpSample')
HeritableList_Noise<-filter(fit_data,cutoff==TRUE & controls=='Controls')
HeritableList_MemorySeq<- HeritableList_MemorySeq[!(HeritableList_MemorySeq$gene %in% HeritableList_Noise$gene),]
GO_List<-unique(cbind(TrimmedData$gene,TrimmedData$GO_Terms))
colnames(GO_List)<-(c('gene','GO_Terms'))
HeritableList_MemorySeq<-merge(HeritableList_MemorySeq,GO_List,by='gene')



Temp<-dcast(data = TrimmedData,formula = number~gene,fun.aggregate = sum,value.var = "rpm")
Labels<-Temp$number
Temp$number<-as.factor(Temp$number)
## for plotting
colors = rainbow(length(unique(Temp$number)))
names(colors) = unique(Temp$number)

## Executing the algorithm on curated data
tsne <- Rtsne(Temp[,-1], dims = 2, perplexity=floor((nrow(Temp) - 1) / 3), verbose=TRUE, max_iter = 10000)
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=Temp$number, col=colors[Temp$number])

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scran")

#Cluster analysis was conducted for both CHO-K1 and HEK datasets using BuildSNNGraph function from scran, setting nearest neighbors k = 20. Then, the Walktrap method (Pons & Latapy, 2005) from igraph computed and identified each cluster using random walks.  To detect cluster marker genes we utilized the findMarkers function from scran which carries out a Welch t-test to perform pairwise comparisons between pairs of clusters (Duò et al., 2020; Soneson & Robinson, 2018), setting pval.type= “all” and direction= “up”. The top 10 candidate marker genes with the lowest p-value for each cluster were detected and reported in Tables 1, 2, and 3 for CHO, and in Tables 4 and 5 for HEK (Supplement 5).
