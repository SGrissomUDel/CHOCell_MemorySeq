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

print("Loading all Packages Now")  
loadPackages <- function(){
  library(topGO)
  library(GO.db)
  library(biomaRt)
  library(Rgraphviz)
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
  library(CliquePercolation) #version 0.3.0
  library(qgraph)            #version 1.6.5
  library(Matrix)            #version 1.2-18
  library(topGO)
  #library(rmarkdown)
}
meltData <- function(ExcelFileNoise, ExcelFileMemory){
  loadPackages()
  setwd('/Users/SpencerGrissom/Downloads/RNA_Seq_Results')
  #setwd('/Users/SpencerGrissom/Downloads')
  getwd()
  #GTF <- readGFF(paste0(getwd(),"/Cricetulus_griseus_picr.CriGri-PICR.104.gtf"))
  GTF <- readGFF(paste0(getwd(),"/GCF_003668045.1_CriGri-PICR_genomic.gff"))
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
  txdb <- makeTxDbFromGFF(file = paste0(getwd(),"/GCF_003668045.1_CriGri-PICR_genomic.gff"), format="gff")
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
  # Filter out samples that have less than 1 million reads
  data <- filter(data, totalMappedReads>0.5*10^6)
  
  # Calculate CVs, fanos, etc.
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
#Input the MemorySeq and Noise Control count tables to generate the relevant
#genes and generate a GO Table for the gene IDs

print("Loading data in Now")  

collapseData <- meltData("Noise_Control_FinalTable.xlsx","Memory_Seq_Samples_FinalTable_Reduced.xlsx")
collapseData <-filter(collapseData,controls=="ExpSample",number==1)
RelevantInfo<-data.frame(cbind(collapseData$gene,collapseData$geneID,collapseData$GO_Terms))
colnames(RelevantInfo)<-c('gene','geneID','GO_Terms')
bg_genes=as.vector(RelevantInfo$gene)



######################################################################
######################################################################
#########     Use this Block for DESeq Results        ################
######################################################################
######################################################################
#Change the Name and the matrix used and the minimum node count
#In the case here, this is for describing the GO enrichment of the stress
#associated genes for the with L2FC>0.58 and p-value<0.05

######################################################################
sampleMatrix<-resultsControl_HighAmmonia_tb
RelevantName<-'HighAmmonia_LFC058Full'
nodeCount<-10
######################################################################

candidate_list =cbind(sampleMatrix$gene,sampleMatrix$log2FoldChange,sampleMatrix$pvalue,sampleMatrix$padj)
candidate_list<-data.frame(candidate_list)
colnames(candidate_list)<-c('gene','log2FoldChange','pvalue','padj')
candidate_list$log2FoldChange <- as.numeric(candidate_list$log2FoldChange)
candidate_list$pvalue <- as.numeric(candidate_list$pvalue)
candidate_list$padj <- as.numeric(candidate_list$padj)
candidate_list<-candidate_list[complete.cases(candidate_list), ]
candidate_list=filter(candidate_list,abs(log2FoldChange)>0.58)
candidate_list=filter(candidate_list,padj<0.05)
#candidate_list<-filter(candidate_list, gene %in% HeritableList_MemorySeq$gene)
candidate_list<-as.vector(candidate_list$gene)
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


######################################################################
######################################################################
####### Use this Block for a List of Genes from GeneID ###############
######################################################################
######################################################################
#This is used to check a series of GeneID, such as those from the
#network communities

######################################################################
# candidate_list<-read_excel("genes_2_assess.xlsx")
# nodeCount<-10
# RelevantName<-'FourthCluster_68_Threshold'
######################################################################

# RelevantInfo$geneID<-as.double(RelevantInfo$geneID)
# candidate_list<-left_join(candidate_list,RelevantInfo,by='geneID')
# candidate_list= as.vector(candidate_list$gene)
# 
# geneList=factor(as.integer(bg_genes %in% candidate_list))
# names(geneList)= bg_genes

######################################################################
######################################################################
#########    Use this Block for a List of Genes Names ################
######################################################################
######################################################################
#This is used to check a series of gene names, such as those from the
#network communities

######################################################################
# candidate_list= as.vector(HeritableList_MemorySeq$gene)
# RelevantName<-'AllHeritableGenes'
# nodeCount<-10
######################################################################
# 
# geneList=factor(as.integer(bg_genes %in% candidate_list))
# names(geneList)= bg_genes


######################################################################
######################################################################
#########     Use this Block to Create GO Table       ################
######################################################################
######################################################################
#Develops a comprehensive list of GO terms from the ENSEMBL Mart

db= useMart('ENSEMBL_MART_ENSEMBL',dataset='cgpicr_gene_ensembl', host="https://apr2020.archive.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)
gene_2_GO_ENSEMBL=unstack(go_ids[,c(1,2)])
GO_List<-cbind(RelevantInfo$gene,RelevantInfo$GO_Terms)
GO_List<-data.frame(GO_List)
colnames(GO_List)<-c('gene','GO_Terms')
gene_2_GO_List<-as.list(strsplit(GO_List$GO_Terms, ", "))
names(gene_2_GO_List)<-GO_List$gene


combinedGO<-mapply(c,gene_2_GO_List,gene_2_GO_ENSEMBL)
combinedGO<-sapply(combinedGO,unique)
combinedGO<-mapply(na.omit,combinedGO)



######################################################################
######################################################################
#########     Run GO Analysis                         ################
######################################################################
######################################################################
#Inputs include GO gene list from 1 of the 3 input methods above, GO tables,
#and number of nodes,

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = combinedGO,nodeSize=nodeCount)



######################################################################
######################################################################
#########     Search Genes with GO Term               ################
######################################################################
######################################################################
#If interested, this allows user to search for genes containing a certain GO term

rna.pp.terms <- genesInTerm(GOdata)[["GO:0043066"]] # get genes associated with term
Heritable_genesWithGO<-intersect(rna.pp.terms,as.list(sapply(HeritableList_MemorySeq,as.character)))




######################################################################
######################################################################
#########     Perform GO Enrichment Analysis          ################
######################################################################
######################################################################
#There are six possible statistical tools that could be used for analyzing the 
#GO enrichment and these include Fisher variants of classic, weighted, and elim.
#Fisher's tests compare the expected number to the observed number. The weighted
#and elim are more conservative and take into consideration the observed hierarchy
#and connected nodes. The Kolmogorov-smirnov tests require a distribution of 
#gene p-values, such as those generated from differential gene expression analysis

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher.weight=runTest(GOdata, algorithm='weight01', statistic='fisher')
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
#resultKS.weight=runTest(GOdata, algorithm='weight01', statistic='ks')
#resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")


######################################################################
######################################################################
#########     Use this Block for a List of Genes      ################
######################################################################
######################################################################
#Execute based on classicFisher

Method<-'classicFisher'
allRes2 <- GenTable(GOdata, classicFisher = resultFisher,
                    orderBy = Method, ranksOf = "classicFisher", topNodes = 50)
allRes2


p.adj=round(p.adjust(allRes2$classicFisher,method="BH"),digits = 4)
all_res_final=cbind(allRes2,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
results.table.p= all_res_final[which(all_res_final$classicFisher<=0.05),]
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
write.table(all_res_final[1:50,],paste0("summary_topGO_analysis_",RelevantName,"_",Method,".csv"),sep=",",quote=FALSE,row.names=FALSE)

pdf(file=paste0('topGOPlot_fullnames_',RelevantName,"_",Method,'.pdf'), height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(resultFisher.classic), useInfo = "all", sigForAll=TRUE, firstSigNodes=5,.NO.CHAR=50)
dev.off()


######################################################################
######################################################################
#########     Use this Block for a List of Genes      ################
######################################################################
######################################################################
#Execute based on weightedFisher

Method<-'weightedFisher'
allRes2 <- GenTable(GOdata, classicFisher = resultFisher, 
                    weightedFisher=resultFisher.weight, 
                    elimFisher= resultFisher.elim, 
                    #classicKS= resultKS, 
                    #weightedKS= resultKS.weight, 
                    #elimKS = resultKS.elim, 
                    orderBy = Method, ranksOf = "classicFisher", topNodes = 50)
allRes2


p.adj=round(p.adjust(allRes2$weightedFisher,method="BH"),digits = 4)
all_res_final=cbind(allRes2,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
results.table.p= all_res_final[which(all_res_final$weightedFisher<=0.05),]
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
write.table(all_res_final[1:50,],paste0("summary_topGO_analysis_",RelevantName,"_",Method,".csv"),sep=",",quote=FALSE,row.names=FALSE)

pdf(file=paste0('topGOPlot_fullnames_',RelevantName,"_",Method,'.pdf'), height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(resultFisher.weight), useInfo = "all", sigForAll=TRUE, firstSigNodes=5,.NO.CHAR=50)
dev.off()


######################################################################
######################################################################
#########     Use this Block for a List of Genes      ################
######################################################################
######################################################################
#Execute based on elimFisher

Method<-'elimFisher'
allRes2 <- GenTable(GOdata, classicFisher = resultFisher, 
                      weightedFisher=resultFisher.weight, 
                      elimFisher= resultFisher.elim, 
                      #classicKS= resultKS, 
                      #weightedKS= resultKS.weight, 
                      #elimKS = resultKS.elim, 
                      orderBy = Method, ranksOf = "classicFisher", topNodes = 50)
allRes2


p.adj=round(p.adjust(allRes2$elimFisher,method="BH"),digits = 4)
all_res_final=cbind(allRes2,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
results.table.p= all_res_final[which(all_res_final$elimFisher<=0.05),]
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
write.table(all_res_final[1:50,],paste0("summary_topGO_analysis_",RelevantName,"_",Method,".csv"),sep=",",quote=FALSE,row.names=FALSE)

pdf(file=paste0('topGOPlot_fullnames_',RelevantName,"_",Method,'.pdf'), height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(resultFisher.elim), useInfo = "all", sigForAll=TRUE, firstSigNodes=5,.NO.CHAR=50)
dev.off()



######################################################################
######################################################################
######### Screen up/downregulated genes for GO Terms  ################
######################################################################
######################################################################
#If using differntially expressed gene data, can look at a specific GO Term
#and identify which genes with that GO term were significant and which were up
#or down regulated

PosApoptosis.genes <- genesInTerm(GOdata)[["GO:0043065"]]
sig.PosApoptosis.genes<-PosApoptosis.genes[PosApoptosis.genes %in% HeritableList$gene]
Up_PosApoptosis.genes<-filter(candidate_list_temp,gene %in% sig.PosApoptosis.genes & log2FoldChange>0)
Down_PosApoptosis.genes<-filter(candidate_list_temp,gene %in% sig.PosApoptosis.genes & log2FoldChange<0)
