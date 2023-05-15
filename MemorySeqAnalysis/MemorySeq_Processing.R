rm(list=ls())
cat("\014")

######################################################################
######################################################################
##############  Functions for Data Processing    #####################
######################################################################
######################################################################
#Functions include melting of data from count tables structured as 
#rows=gene count number, columns=replicate number, filtering functions,
#functions for summarizing data, and plotting gene expression

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
}
meltData <- function(ExcelFileNoise, ExcelFileMemory,GCF_FileName,GOTermsID_FileName,WDAddress,ProteinCoding){
  loadPackages()
  setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
  getwd()
  GTF <- readGFF(paste0(getwd(),"/",GCF_FileName))
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
  if(ProteinCoding=="Yes"){
  collapseData <- filter(collapseData, gene_biotype=='protein_coding')
  }
  collapseData <- collapseData %>% separate(Dbxref,into=c("TBR","geneID"), remove=FALSE, sep=":")
  txdb <- makeTxDbFromGFF(file = paste0(getwd(),"/",GCF_FileName))
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
  collapseData$TBR <- NULL
  collapseData$Dbxref <- NULL
  GO_ID <- read_excel(paste0(getwd(),"/",GOTermsID_FileName))
  collapseData$geneID <- as.double(collapseData$geneID)
  collapseData <- left_join(collapseData, GO_ID,by='geneID')
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
  data_biotype<-unique(data[,c('gene','gene_biotype')])
  data_with_metrics<-left_join(data_with_metrics, data_biotype, by='gene')
  colnames(data_with_metrics) <- c('gene','controls','geneID',
                                   'CV_rpm','CV_tpm','Mean_rpm','Mean_tpm', 'skewness','kurtosis',
                                   'kurtosis_e1071', 'CV_log2tpm', 'Mean_log2tpm','CV_log2_tpm_plus1','Mean_log2_tpm_plus1','gene_biotype')
  return(data_with_metrics)
}


######################################################################
######################################################################
##############   Collect Data and Summarize      #####################
######################################################################
######################################################################
#Change out the NoiseControl and MemorySeq files names to the xlsx count table
#as well as working directory, GCF file, and GO ID terms
#TrimmedData filters data to be TPM > 1.5 and summary 
#includes all the Mean, CoV, skew, etc.


print("Loading data in Now")  

######################################################################
WorkingDirectory<-'/Users/SpencerGrissom/Downloads/RNA_Seq_Results'
Noise_FileName<-"Noise_Control_FinalTable.xlsx"
MemorySeq_FileName<-"Memory_Seq_Samples_FinalTable.xlsx"
GFF_FileName<-"GCF_003668045.1_CriGri-PICR_genomic.gff"
GOTermsID_FileName<-"Final_GO_ID_CriGriPICR.xlsx"
ProteinCodingOnly<-"No"
######################################################################

collapseData <- meltData(Noise_FileName,MemorySeq_FileName,GFF_FileName,GOTermsID_FileName,WorkingDirectory,ProteinCodingOnly)
TrimmedData <- remove_genes_less1_5tpm(collapseData)
summary <- calc_variability_metrics(TrimmedData)


######################################################################
######################################################################
#########Poisson Regression Fitting for Heritability##################
######################################################################
######################################################################
#Identification for heritability, two most important parameters include 
#percentile for deviation from the fit and the mean TPM cut-off
#This section also extracts the heritable genes

print("Performing Poisson regression fitting and plotting data")  

######################################################################
percentile_cutoff <- 0.98
mean_tpm_cutoff <- 2.5
######################################################################

exp_sample_data <- filter(summary, controls=='ExpSample')
control_data <- filter(summary, controls=='Controls')
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
                           resid  = resid(model_exp_sample),
                           gene_biotype = exp_sample_data$gene_biotype)
fit_data_controls <- data.frame(gene = control_data$gene, 
                                #GeneSymbol = control_data$GeneSymbol,
                                controls = control_data$controls,
                                Mean_tpm = control_data$Mean_tpm,
                                CV_tpm = control_data$CV_tpm,
                                fitted_y= fitted(model_controls),
                                resid = resid(model_controls),
                                gene_biotype = control_data$gene_biotype)
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

######################################################################
######################################################################
#########     Plot Heritable Genes and Residuals    ##################
######################################################################
######################################################################


print("Plotting Heritability Results Now")  

######################################################################
CallOutGenes<-c('Hmox1','Serpine1','Ier3','Nnmt','Tp63','Hmgcs2')
######################################################################

######################################################################
#Plots the scatter of Noise/MemorySeq samples relative to the 
#fitted poisson regression
plot1<-ggplot(fit_data_controls,aes(x=log2(Mean_tpm),y=CV_tpm)) +
  stat_bin2d(binwidth = c(0.2, 0.07)) +
  xlim(-1,15) +
  ylim(-0.2,4) +
  geom_line(aes(log2(Mean_tpm),fitted_y),size=5,color = "black") +
  scale_fill_continuous(low='dimgray',high="navyblue") + 
  theme_classic() +
  theme(
    axis.title.x = element_text(color = "black", size = 56, face = "bold"),
    axis.title.y = element_text(color = "black", size = 56, face = "bold"),
    axis.text.x = element_text(color = "black", size = 48, face = "bold"),
    axis.text.y = element_text(color = "black", size = 48, face = "bold")
  ) +
  labs(x = expression(log[2](Mean[TPM])), y = expression(CV[TPM]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=32), #change legend title font size
        legend.text = element_text(size=32)) #change legend text font size
plot2<-ggplot(fit_data_exp,aes(x=log2(Mean_tpm),y=CV_tpm)) +
  stat_bin2d(binwidth = c(0.2, 0.07)) +
  xlim(-1,15) +
  ylim(-0.2,4) +
  geom_line(aes(log2(Mean_tpm),fitted_y),size=5,color = "black") +
  scale_fill_continuous(low='dimgray',high="navyblue") + 
  theme_classic() +
  theme(
    axis.title.x = element_text(color = "black", size = 56, face = "bold"),
    axis.title.y = element_text(color = "black", size = 56, face = "bold"),
    axis.text.x = element_text(color = "black", size = 48, face = "bold"),
    axis.text.y = element_text(color = "black", size = 48, face = "bold")
  ) +
  labs(x = expression(log[2](Mean[TPM])), y = expression(CV[TPM]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=32), #change legend title font size
        legend.text = element_text(size=32)) #change legend text font size
p<-plot_grid(plot1, plot2, labels = "AUTO")
save_plot("Control_MemorySeq_PoissonFit.png", p, ncol = 2, base_height=18, base_width=24,limitsize=FALSE)


######################################################################
#Plots the scatter of Noise/MemorySeq samples with Poisson fit and 
#includes genes that are emphasized
plot3<-ggplot(filter(fit_data,cutoff == FALSE), aes(x=log2(Mean_tpm),y=CV_tpm)) +
  stat_bin2d(binwidth = c(0.1, 0.02)) +
  scale_fill_continuous(low='dimgray',high="navyblue") + 
  xlim(-1,15) +
  ylim(-0.2,3) +
  geom_point(data=filter(fit_data,cutoff == TRUE),
             aes(x=log2(Mean_tpm), y=CV_tpm),color='green',size=5) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(color = "black", size = 72, face = "bold"),
    axis.title.y = element_text(color = "black", size = 72, face = "bold"),
    axis.text.x = element_text(color = "black", size = 48, face = "bold"),
    axis.text.y = element_text(color = "black", size = 48, face = "bold")
  ) +
  labs(x = expression(log[2](Mean[TPM])), y = expression(CV[TPM]))+
  theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
        legend.key.height = unit(4.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=72), #change legend title font size
        legend.text = element_text(size=72),
        strip.text = element_text(size = 72))+ #change legend text font size
  geom_text(data=filter(fit_data,gene %in% CallOutGenes),
            aes(x=log2(Mean_tpm), y=CV_tpm+0.125,label=gene), size=16, color='red',fontface="bold",position=position_jitter(width=0,height=0.15)) + 
  geom_point(data=filter(fit_data,gene %in% CallOutGenes),
             aes(x=log2(Mean_tpm), y=CV_tpm),color='red',size=5) + 
  facet_wrap(~controls)
save_plot("SignificantGenes.png", plot3, ncol = 1, base_height=24, base_width=36,limitsize=FALSE)
  
######################################################################
#Distribution of residuals for noise and MemorySeq samples
plot4<-ggplot(fit_data, aes(x=resid,fill=controls,color=controls))+
  geom_histogram(bins=100,position="identity", alpha=0.5)+
  geom_vline(xintercept = cutoff099,size=6)+
  geom_rug()+
  theme_classic()+
  theme(
    axis.title.x = element_text(color = "black", size = 72, face = "bold"),
    axis.title.y = element_text(color = "black", size = 72, face = "bold"),
    axis.text.x = element_text(color = "black", size = 72, face = "bold"),
    axis.text.y = element_text(color = "black", size = 72, face = "bold")) +
  labs(x = expression(Residuals), y = expression(Count))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=38), #change legend title font size
        legend.text = element_text(size=38))

save_plot("Histogram.png", plot4, ncol = 1, base_height=18, base_width=20,limitsize=FALSE)
 
######################################################################
#Distribution of skewness for each gene
plotSkew<-ggplot(summary, aes(x=skewness,fill=controls,color=controls))+
  geom_histogram(bins=100,position="identity", alpha=0.5)+
  geom_rug()+
  theme_classic()+
  theme(
    axis.title.x = element_text(color = "black", size = 64, face = "bold"),
    axis.title.y = element_text(color = "black", size = 64, face = "bold"),
    axis.text.x = element_text(color = "black", size = 52, face = "bold"),
    axis.text.y = element_text(color = "black", size = 52, face = "bold")) +
  labs(x = expression(Heritability~Index), y = expression(Count))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=38), #change legend title font size
        legend.text = element_text(size=38))
save_plot("HistogramSkewness.png", plotSkew, ncol = 1, base_height=18, base_width=20,limitsize=FALSE)


######################################################################
#Representative plots for Heritable/Non-heritable gene expression
Heritable_Genes<-HeritableList_MemorySeq
Top6_Resid <- Heritable_Genes[(nrow(Heritable_Genes)-5):nrow(Heritable_Genes),'gene']
NonHeritable_Genes <- filter(fit_data[order(fit_data$resid),],controls=='ExpSample' & cutoff=='FALSE')
mid_resid <- median.default(NonHeritable_Genes$resid)
NonHeritable_Genes<-NonHeritable_Genes %>%
  dplyr::mutate(dist_mid = abs(resid-mid_resid)) %>%
  dplyr::arrange(dist_mid)

Mid6_Resid <- NonHeritable_Genes[1:6,'gene']

genesToPlot <- filter(TrimmedData,
                      gene%in%CallOutGenes | gene%in%Mid6_Resid)
genes<-c('Gapdh','Hmox1','Nmi','Serpine1','Rras2','Ier3','Fbxo31','Nnmt','Vps26a','Tp63','Znf512b','Hmgcs2')
genesToPlot$gene <- factor(genesToPlot$gene, genes)
plot5<-ggplot(genesToPlot,
       aes(x=log2(tpm+1), color=controls, fill=controls))+
  geom_density(alpha=0.5)+
  xlim(0,14)+
  facet_wrap(~gene,nrow=6, scales = 'free_y')+
  theme_classic()+
  geom_rug()+
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 80))+
  theme(
    axis.title.x = element_text(color = "black", size = 72, face = "bold"),
    axis.title.y = element_text(color = "black", size = 72, face = "bold"),
    axis.text.x = element_text(color = "black", size = 32, face = "bold"),
    axis.text.y = element_text(color = "black", size = 48, face = "bold")) +
  labs(x = expression(log[2](TPM+1)), y = expression(Density))
  
save_plot("SampleHistogram.png", plot5, ncol = 1, base_height=24, base_width=36,limitsize=FALSE)

######################################################################
#Ratio of the coefficient of variation for the noise vs memorySeq samples
plotCV_exp <- filter(fit_data,controls=='ExpSample')
plotCV_con <- filter(fit_data,controls=='Controls')
plotCV <- merge(plotCV_exp,plotCV_con,by='gene')

plotRatio<-ggplot(plotCV,aes(x=CV_tpm.y,y=CV_tpm.x))+
  xlim(0,0.7) +
  ylim(0,3) +
  stat_bin2d(binwidth = c(0.005, 0.02)) +
  geom_abline(intercept = 0, slope = 1)+
  #stat_bin2d(binwidth = c(0.05*2, 0.01*2)) +
  geom_point(data=filter(plotCV,cutoff.x == TRUE),
             aes(x=(CV_tpm.y), y=CV_tpm.x),color='blue',size=10) + 
  geom_point(data=filter(plotCV,cutoff.y == TRUE),
              aes(x=(CV_tpm.y), y=CV_tpm.x),color='red',size=10)+   
  theme_classic() +
  theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
        legend.key.height = unit(4.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=72), #change legend title font size
        legend.text = element_text(size=72),
        strip.text = element_text(size = 72))+ #change legend text font size
  theme(
    axis.title.x = element_text(color = "black", size = 72, face = "bold"),
    axis.title.y = element_text(color = "black", size = 72, face = "bold"),
    axis.text.x = element_text(color = "black", size = 72, face = "bold"),
    axis.text.y = element_text(color = "black", size = 72, face = "bold")
  ) +
  labs(x = expression(CV[Noise~Control]), y = expression(CV[MemorySeq~Samples]))
save_plot("RatioCV.png", plotRatio, ncol = 1, base_height=24, base_width=36,limitsize=FALSE)

######################################################################
######################################################################
#########Clustering Functions and Heatmap/Correlations################
######################################################################
######################################################################
#Series of functions required for generating Pearson pairwise correlations
#and plotting the resulting heatmap

pearsonDist <- function(matdata){
  return(as.dist(1-cor(t(matdata))))
}
pearsonCompleteClustering <- function(x){
  hclust(pearsonDist(x))
}
Grouping <- function(hclustObj, groupsNeeded) {
  groupAssignment <- dendroextras::slice(
    hclustObj, k = groupsNeeded)
  # restore to same order as data that was clustered.
  groupAssignment <- groupAssignment[order(hclustObj$order)]
  
  splitRowsByGroup <- function(dat){
    if (is.matrix(dat)) {
      splitDfs <- split(as.data.frame(dat), groupAssignment)
      splitdata <- lapply(splitDfs, as.matrix)
    } else {
      splitdata <- split(dat, groupAssignment)
    }
    return(splitdata)
  }
  coloring <- function(colorChoices = rainbow(7)) {
    numColors <- length(colorChoices)
    cycledIfNecessary <- ((groupAssignment-1) %% numColors) + 1
    geneColors <- colorChoices[cycledIfNecessary]
    return(geneColors)
  }
  
  output <- list(splitRowsByGroup = splitRowsByGroup,
                 coloring = coloring)
  return(output)
}
geneClusterer <- pearsonCompleteClustering
colorsForGroups <- brewer.pal(n=3,name='Paired') 
color4heatmap <- rev(brewer.pal(11,'RdYlBu'))
plotScatterGenePair <- function(gene1, gene2, dataset){
  temp <- filter(dataset, gene==gene1 | gene==gene2)
  temp_cast <- dcast(temp, 'number + controls ~ gene', value.var = 'tpm',
                     fun.aggregate = max)
  xname <- colnames(temp_cast)[4]
  yname <- colnames(temp_cast)[3]
  colnames(temp_cast) <- c('number','controls','gene1','gene2')
  q <- ggplot(temp_cast, aes(x = gene1, 
                             y = gene2, 
                             color=controls)) + 
    geom_point() + theme_classic()+
    ylab(xname)+
    xlab(yname)
  return(q)
}
calculate_correlation_coeff <- function(gene1, gene2, dataset){
  gene1 <- filter(dataset, gene==gene1 & controls == 'ExpSample')
  gene2 <- filter(dataset, gene==gene2 & controls == 'ExpSample')
  R <- cor(gene1$tpm, gene2$tpm)
  return(R)
}
makeRects <- function(cells){
  nr <- nrow(cells)
  nc <- ncol(cells)
  coords = expand.grid(nr:1, 1:nc)[cells,]
  xl=coords[,2]-0.25
  yb=coords[,1]-0.25
  xr=coords[,2]+0.25
  yt=coords[,1]+0.25
  rect(xl,yb,xr,yt,border="black",lwd=1)
}
corrMatrix <- function(gene_list, dataSet,cor_function){ 
  
  gene_list2 <- gene_list[,c('gene','CV_tpm')]
  
  temp <- dplyr::left_join(gene_list2, 
                           dataSet, 
                           by = 'gene')
  temp_exp <- filter(temp, controls == 'ExpSample')
  temp_controls <- filter(temp, controls == 'Controls')
  
  data_for_corrmat <- data.frame(temp_exp)
  data_for_corrmat_cast <- acast(data = data_for_corrmat, 
                                 formula = 'number~gene', 
                                 value.var = 'tpm',
                                 fun.aggregate = mean)
  test <- cor(data_for_corrmat_cast, method = cor_function)
}

######################################################################
######################################################################
#########  Create Correlation Matrix and Plot         ################
######################################################################
######################################################################
#Running the correlation functions and also highlighting the most positive and
#negative correlation pair for demonstrating the co-fluctuations

print("Preparing Correlation Matrix Now")

test <- corrMatrix(dataSet = TrimmedData, 
                          gene_list = filter(Heritable_Genes, cutoff==TRUE & controls=='ExpSample'), 
                          cor_function = 'pearson')
num_groups = 10
  RowCount <- nrow(test)
  ColCount <- ncol(test)
  geneClustering_expSample <- geneClusterer(test)
  ForPlotMin=which(test == min(test), arr.ind = TRUE)[1,]
  ForPlotMax=which(abs(test-0.95) == min(abs(test-0.95)), arr.ind = TRUE)[1,]
  selection <- matrix(rep(F), nrow=RowCount, ncol=ColCount)
  selection[ForPlotMin[1]:RowCount,ForPlotMin[2]] <- T
  selection[ForPlotMin[1],ForPlotMin[2]:ColCount] <- T
  selection[ForPlotMax[1]:RowCount,ForPlotMax[2]] <- T
  selection[ForPlotMax[1],ForPlotMax[2]:ColCount] <- T
  labelCol <- rep(NA, ColCount)
  labelCol[c(ForPlotMin[2],ForPlotMax[2])] <- colnames(test)[c(ForPlotMin[2],ForPlotMax[2])]
  labelRow <- rep(NA, RowCount)
  labelRow[c(ForPlotMin[1],ForPlotMax[1])] <- row.names(test)[c(ForPlotMin[1],ForPlotMax[1])]
  geneGrouping <- Grouping(geneClustering_expSample, 
                           groupsNeeded = num_groups) # changed for # groups
  colorsForGroups <- brewer.pal(n=6,name='Paired')
  
  print("Plotting Correlation Heatmap Now")  
  heatmap <- heatmap.2(t(as.matrix(test)), scale = 'none',
                                 #add.expr=makeRects(selection),
                                 #labCol=colnames(test),
                                 #labRow=row.names(test),
                                 labCol= labelCol,
                                 labRow=labelRow,
                                 Colv = as.dendrogram(geneClustering_expSample), 
                                 Rowv = as.dendrogram(geneClustering_expSample),
                                 dendrogram = 'both',
                                 ColSideColors = geneGrouping$coloring(colorsForGroups),
                                 col = color4heatmap,
                                 density.info='density', denscol="black", 
                                 symkey=FALSE, trace='none',
                                 margins = c(6,10), cexCol=1, cexRow=1,
                                 main="Correlations between all pairs of heritable genes")

MinGene1 <- row.names(test)[ForPlotMin[2]]
MinGene2 <- colnames(test)[ForPlotMin[1]]

MaxGene1 <- row.names(test)[ForPlotMax[2]]
MaxGene2 <- colnames(test)[ForPlotMax[1]]

plotScatterGenePair(MaxGene1,MaxGene2,TrimmedData)

######################################################################
######################################################################
#########             Sensitivty Analysis             ################
######################################################################
######################################################################
#Functions and execution of the Cooks Distance analysis for measuring the data's
#senstivity to outliers significantly influencing the measured correlations

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
plotCooks <- function(data,corMatrix,correlationThresh,gene1,gene2){
  temp_exp <- data.frame(filter(data, controls=='ExpSample'))
  seq_flat <- flattenCorrMatrix(corMatrix)
  seq_flat <- dplyr::filter(seq_flat,cor>correlationThresh)
  #seq_flat <- dplyr::filter(seq_flat,row==gene1 & column==gene2)

  for (i in 1:dim(seq_flat)[1]){
    genes <- dplyr::filter(temp_exp, gene==seq_flat$row[i] | gene==seq_flat$column[i])
    genes_cast <- data.frame(acast(genes, formula = 'number~gene', 
                                  value.var = 'tpm',
                                  fun.aggregate = mean))
    colnames(genes_cast) <- c('gene1', 'gene2')
   fit <- lm(gene1 ~ gene2 + 0,
              data = genes_cast)
    par(mfrow = c(2,2))
    pdf(paste0(getwd(),
               seq_flat$row[i], '_', seq_flat$column[i], '.pdf'))
    plot(fit)
    dev.off()
  }
}

plotCooks(TrimmedData,test,0.97,'Gnl1','Mepce')


######################################################################
######################################################################
#########            k-Clique Percolation             ################
######################################################################
######################################################################
#Uses the correlation matrix and the CliquePercolation package to generate the
#community network clusters. Can change the k and threshold value, which can
#be optimized using cpThreshold. Then colllects the genes in each cluster in a
#generated text file

library(CliquePercolation) #version 0.3.0
library(qgraph)            #version 1.6.5
library(Matrix)            #version 1.2-18


Correlation<-test
start_time <- Sys.time()
W <- qgraph::qgraph(Correlation, theme = "colorblind", layout = "spring", cut = 0.6,threshold=0.2)
#thresholds <- cpThreshold(W, method = "weighted", k.range = c(4),
#                          I.range = c(0.7),
#                          threshold = c("largest.components.ratio","chi",'entropy'))
end_time <- Sys.time()
end_time-start_time
#thresholds
dev.off()
result<-cpAlgorithm(W, k = 4, method = "weighted", I = 0.7)
png(file="K_Clique_Cluster",
    width=1800, height=1600)
g3 <- cpColoredGraph(W, list.of.communities = result$list.of.communities.numbers,layout = "spring", theme='colorblind',
                    labels=row.names(Correlation), vsize=3, cut=0.95, border.width=0.5, 
                     border.color='black',legend.cex=.35,
                     edge.width = 0.01,alpha=0.5, title ="Network Community Identification")
dev.off()

Conversion<-W$graphAttributes$Nodes$labels
Conversion<-data.frame(Conversion)
Conversion <-data.frame(cbind(row.names(Conversion),Conversion[,1]))
colnames(Conversion) <- c('gene','label')
for(i in 1:length(result$list.of.communities.labels)){
  temp<-result$list.of.communities.labels[[i]]
  temp<-filter(Conversion,label %in% temp)
  temp<-left_join(temp,TrimmedData,by='gene')
  temp<-filter(temp,controls=='ExpSample',number==1)
  write.table(data.frame(temp$geneID),file='test.txt',append=TRUE, row.names=FALSE)
  cat('\n',file='test.txt',append = TRUE)  
}
