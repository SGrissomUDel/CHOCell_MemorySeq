library(ggplot2)
library(dplyr)
####################################################################################
####################################################################################
################ Sorting Data Points & Characteristic Functions ####################
####################################################################################
####################################################################################
maxY <- 10
logBase <- 10
foldchangeThresh <- 0.58
pValThresh <- .05

getColor <- function(log2FoldChange,padj) {
  
  ifelse( log2FoldChange < -foldchangeThresh & -log(padj, base = logBase) > -log(pValThresh, base = logBase),
          "blue",
          ifelse( log2FoldChange > foldchangeThresh & -log(padj, base = logBase) > -log(pValThresh, base = logBase),
                  "red",
                  "gray")
  )
}

getColorLabelName <- function(log2FoldChange,padj) {
  
  ifelse(log2FoldChange < -foldchangeThresh & -log(padj, base = logBase) > -log(pValThresh, base = logBase),
          "Down-Regulated",
          ifelse(log2FoldChange > foldchangeThresh & -log(padj, base = logBase) > -log(pValThresh, base = logBase),
                  "Up-Regulated",
                  "Insignificant")
  )
}

getShapeLabelName <- function(padj) {
  ifelse( (-log(padj, base = logBase) < maxY), 'On Plot', 'Off Plot')
}

getShape <- function(padj) {
  ifelse( (-log(padj, base = logBase) < maxY), 16, 2)
}

getDispY <- function(padj) {
  ifelse( (-log(padj, base = logBase) < maxY), -log(padj, base = logBase), maxY )
}

getSize <- function(padj) {
  1
}

####################################################################################
####################################################################################
########################## Creating Volcano Plot Function ##########################
####################################################################################
####################################################################################

createVolcanoPlot <- function(dataSet) {
  
  ggplot(data = dataSet) +
    geom_point(aes(x=log2FoldChange, 
                   y=getDispY(padj),
                   color = getColorLabelName(log2FoldChange,padj),
                   ),
    shape = getShape(dataSet$padj),
    size = getSize(dataSet$padj)
    
    ) +
    geom_vline(xintercept=c(-foldchangeThresh,foldchangeThresh), color = 'red', linetype = 'dashed') +
    geom_hline(yintercept=(-log(pValThresh, base = logBase)), color = 'red', linetype = 'dashed') +
    ylab(expression('-Log'[10]*' P-values')) +
    xlab(expression('Log'[2]*' Fold Change')) +
    xlim(-6.5,6.5) +
    ylim(0,maxY) +
    labs(
      color = 'Regulation',
      shape = 'Shapes',
    ) +
    scale_color_manual(values = c('blue','gray','red')) +
    theme_classic()
  
}

createVolcanoPlot(All_Every)