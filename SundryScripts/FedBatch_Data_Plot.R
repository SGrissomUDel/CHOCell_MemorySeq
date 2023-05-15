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
}
loadPackages()

my_data <- read_excel("StressCondition_FedBatch_Data.xlsx", sheet = "Productivity_Data")
Titer <- summarySE(my_data, measurevar="Titer", groupvars=c("Group","Day","Salt"))
Titer_Low<-filter(Titer,Salt=='Low')
Titer_High<-filter(Titer,Salt=='High')
my_data_VCD <- read_excel("StressCondition_FedBatch_Data.xlsx", sheet = "Growth_Data")
VCD <- summarySE(my_data_VCD, measurevar="VCD", groupvars=c("Group","Day","Salt"))
VCD_Low<-filter(VCD,Salt=='Low')
VCD_High<-filter(VCD,Salt=='High')

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

pd <- position_dodge(.3)

Titer_Low_Plot<-ggplot(Titer_Low, aes(x=Day, y=Titer, group=Group, color=Group)) +
  geom_errorbar(
    aes(ymin = Titer - se, ymax = Titer + se),
    width = 2,
    size = 3,
    position = pd
  ) +
  ylim(0,1.8) +
  geom_line(position = pd, size=5,linetype='dashed') +
  geom_point(position = pd, size = 10)+
  theme_classic() + 
  theme(
    axis.title.x = element_text(color = "black", size = 114, face = "bold"),
    axis.title.y = element_text(color = "black", size = 114, face = "bold"),
    axis.text.x = element_text(color = "black", size = 88, face = "bold"),
    axis.text.y = element_text(color = "black", size = 88, face = "bold")
  ) +
  labs(x = expression(Day), y = expression(Titer~(g/L)))+
  theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
        legend.key.height = unit(4.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=72), #change legend title font size
        legend.text = element_text(size=72),
        strip.text = element_text(size = 72))+ #change legend text font size
  scale_color_manual(
    limits = c('Control','Amm Stress', 'Lac Stress','Amm/Lac Stress'),
    labels  =c('Control', 'Ammonia Stress', 'Lactate Stress', 'Ammonia/Lactate Stress'),
    values = c("black",
               "blue",
               "red",
               "purple"))
  
  save_plot("Titer_data_LowSalt.png", Titer_Low_Plot, ncol = 1, base_height=24, base_width=50,limitsize=FALSE)

  Titer_High_Plot<-ggplot(Titer_High, aes(x=Day, y=Titer, group=Group, color=Group)) +
    geom_errorbar(
      aes(ymin = Titer - se, ymax = Titer + se),
      width = 2,
      size = 3,
      position = pd
    ) +
    ylim(0,1.8) +
    geom_line(position = pd, size=5,linetype='dashed') +
    geom_point(position = pd, size = 10)+
    theme_classic() + 
    theme(
      axis.title.x = element_text(color = "black", size = 114, face = "bold"),
      axis.title.y = element_text(color = "black", size = 114, face = "bold"),
      axis.text.x = element_text(color = "black", size = 88, face = "bold"),
      axis.text.y = element_text(color = "black", size = 88, face = "bold")
    ) +
    labs(x = expression(Day), y = expression(Titer~(g/L)))+
    theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
          legend.key.height = unit(4.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_text(size=72), #change legend title font size
          legend.text = element_text(size=72),
          strip.text = element_text(size = 72))+ #change legend text font size
    scale_color_manual(
      limits = c('Control','Amm Stress', 'Lac Stress','Amm/Lac Stress'),
      labels  =c('Control', 'Ammonia Stress', 'Lactate Stress', 'Ammonia/Lactate Stress'),
      values = c("black",
                 "blue",
                 "red",
                 "purple"))
  
  save_plot("Titer_data_HighSalt.png", Titer_High_Plot, ncol = 1, base_height=24, base_width=50,limitsize=FALSE)
  
  
  VCD_High_Plot<-ggplot(VCD_High, aes(x=Day, y=VCD, group=Group, color=Group)) +
    geom_errorbar(
      aes(ymin = VCD - se, ymax = VCD + se),
      width = 2,
      size = 3,
      position = pd
    ) +
    ylim(0,1.2e7) +
    geom_line(position = pd, size=5,linetype='dashed') +
    geom_point(position = pd, size = 10)+
    theme_classic() + 
    theme(
      axis.title.x = element_text(color = "black", size = 114, face = "bold"),
      axis.title.y = element_text(color = "black", size = 114, face = "bold"),
      axis.text.x = element_text(color = "black", size = 88, face = "bold"),
      axis.text.y = element_text(color = "black", size = 88, face = "bold")
    ) +
    labs(x = expression(Day), y = expression(VCD (cells/mL)))+
    theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
          legend.key.height = unit(4.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_text(size=72), #change legend title font size
          legend.text = element_text(size=72),
          strip.text = element_text(size = 72))+ #change legend text font size
    scale_color_manual(
      limits = c('Control','Amm Stress', 'Lac Stress','Amm/Lac Stress'),
      labels  =c('Control', 'Ammonia Stress', 'Lactate Stress', 'Ammonia/Lactate Stress'),
      values = c("black",
                 "blue",
                 "red",
                 "purple"))
  
  save_plot("VCD_data_HighSalt.png", VCD_High_Plot, ncol = 1, base_height=24, base_width=50,limitsize=FALSE)

  VCD_Low_Plot<-ggplot(VCD_Low, aes(x=Day, y=VCD, group=Group, color=Group)) +
    geom_errorbar(
      aes(ymin = VCD - se, ymax = VCD + se),
      width = 2,
      size = 3,
      position = pd
    ) +
    ylim(0,1.2e7) +
    geom_line(position = pd, size=5,linetype='dashed') +
    geom_point(position = pd, size = 10)+
    theme_classic() + 
    theme(
      axis.title.x = element_text(color = "black", size = 114, face = "bold"),
      axis.title.y = element_text(color = "black", size = 114, face = "bold"),
      axis.text.x = element_text(color = "black", size = 88, face = "bold"),
      axis.text.y = element_text(color = "black", size = 88, face = "bold")
    ) +
    labs(x = expression(Day), y = expression(VCD (cells/mL)))+
    theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
          legend.key.height = unit(4.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_text(size=72), #change legend title font size
          legend.text = element_text(size=72),
          strip.text = element_text(size = 72))+ #change legend text font size
    scale_color_manual(
      limits = c('Control','Amm Stress', 'Lac Stress','Amm/Lac Stress'),
      labels  =c('Control', 'Ammonia Stress (10 mM)', 'Lactate Stress (15 mM)', 'Ammonia/Lactate Stress'),
      values = c("black",
                 "blue",
                 "red",
                 "purple"))
  
  save_plot("VCD_data_LowSalt.png", VCD_Low_Plot, ncol = 1, base_height=24, base_width=50,limitsize=FALSE)
  
  
  my_data_summary <- read_excel("StressCondition_FedBatch_Data.xlsx", sheet = "Summary_Data")
  SummaryData_Titer <- summarySE(my_data_summary, measurevar="Titer", groupvars=c("Group","Salt"))
  SummaryData_IVCD <- summarySE(my_data_summary, measurevar="IVCD", groupvars=c("Group","Salt"))
  SummaryData_TiterIVCD<-left_join(SummaryData_Titer,SummaryData_IVCD,by='Group')
  SummaryData_SP <- summarySE(my_data_summary, measurevar="SpecificProductivity", groupvars=c("Group","Salt"))

                               # Change ordering manually
  SummaryData_TiterIVCD<-SummaryData_TiterIVCD %>% arrange(factor(Group, levels = c("Control", "Amm", "Lac", "Amm/Lac", "Osmo", "Osmo/Amm", "Osmo/Lac", "Osmo/Amm/Lac")))
  SummaryData_TiterIVCD$Group<-factor(SummaryData_TiterIVCD$Group, levels = c("Control", "Amm", "Lac", "Amm/Lac", "Osmo", "Osmo/Amm", "Osmo/Lac", "Osmo/Amm/Lac"))
  
  SummaryData_SP<-SummaryData_SP %>% arrange(factor(Group, levels = c("Control", "Amm", "Lac", "Amm/Lac", "Osmo", "Osmo/Amm", "Osmo/Lac", "Osmo/Amm/Lac")))
  SummaryData_SP$Group<-factor(SummaryData_SP$Group, levels = c("Control", "Amm", "Lac", "Amm/Lac", "Osmo", "Osmo/Amm", "Osmo/Lac", "Osmo/Amm/Lac"))
  
  SummaryTiterIVCD<-ggplot(SummaryData_TiterIVCD)+
    geom_bar( aes(x=Group, y=IVCD), stat="identity", colour='black',fill="darkgray", alpha=0.7) +
    geom_errorbar( aes(x=Group, ymin=IVCD-se.y, ymax=IVCD+se.y), width=0.1, colour="black", alpha=0.9, size=5)+
  theme_classic() + 
    geom_line(aes(x=Group, y=Titer*max(IVCD)),stat='identity', size=5,colour='blue',group=1) +
    geom_point(aes(x=Group, y=Titer*max(IVCD)), size = 15,colour='blue')+
    geom_errorbar( aes(x=Group, ymin=Titer*max(IVCD)-se.x*max(IVCD), ymax=Titer*max(IVCD)+se.x*max(IVCD)), width=0.1, colour="darkblue", alpha=0.9, size=5)+
    scale_y_continuous(sec.axis = sec_axis(~./max(SummaryData_TiterIVCD$IVCD)))+
    theme(
      axis.title.x = element_text(color = "black", size = 114, face = "bold"),
      axis.title.y = element_text(color = "black", size = 114, face = "bold"),
      axis.text.x = element_text(color = "black", size = 100, face = "bold",angle = 45,hjust=1),
      axis.text.y = element_text(color = "black", size = 88, face = "bold")
    ) +
    labs(x = expression(Stress~Condition), y = expression(IVCD (cells*days/mL)))+
    scale_y_continuous(
      
      # Features of the first axis
      name = "IVCD (cells*day/mL)",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~./max(SummaryData_TiterIVCD$IVCD),name="Titer (g/L)")
    )+
    theme(legend.key.size = unit(4.5, 'cm'), #change legend key size
          legend.key.height = unit(4.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_text(size=72), #change legend title font size
          legend.text = element_text(size=72),
          strip.text = element_text(size = 72))+ #change legend text font size
    theme( axis.line.y.right = element_line(color = "blue"), 
           axis.ticks.y.right = element_line(color = "blue"))
  save_plot("Summary_TiterIVCD.png", SummaryTiterIVCD, ncol = 1, base_height=35, base_width=50,limitsize=FALSE)
  
  SummarySP<-ggplot(SummaryData_SP)+
    geom_bar( aes(x=Group, y=SpecificProductivity), stat="identity", colour='black',fill="darkgray", alpha=0.7) +
    geom_errorbar( aes(x=Group, ymin=SpecificProductivity-se, ymax=SpecificProductivity+se), width=0.1, colour="black", alpha=0.9, size=5)+
    theme_classic() + 
  theme(
      axis.title.x = element_text(color = "black", size = 114, face = "bold"),
      axis.title.y = element_text(color = "black", size = 114, face = "bold"),
      axis.text.x = element_text(color = "black", size = 100, face = "bold",angle = 45,hjust=1),
      axis.text.y = element_text(color = "black", size = 88, face = "bold")
    ) +
    labs(x = expression(Stress~Condition), y = expression(Specific~Productivity~(pcd)))
   
save_plot("Summary_SP.png", SummarySP, ncol = 1, base_height=35, base_width=50,limitsize=FALSE)

  