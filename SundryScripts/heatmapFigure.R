#
# Installing Packages
#

install.packages('devtools')
library(devtools)

install_github('jokergoo/ComplexHeatmap')
library(ComplexHeatmap)

library(circlize)

##
## (Old) Finding heritable gene list
## Do not run if you already have these variables defined
##

#HeritableList_MemorySeq <- read.delim('C:\\Users\\zdixo\\Desktop\\RAS\\HeritableGenesCHO.txt');
#HeritableList <- HeritableList_MemorySeq

#
# Setting Up Biological Processes
#

GOdata <- new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = combinedGO,nodeSize=nodeCount)

##############################
#### Getting an Average ######
##############################

Avg_ControlZScore <- rowMeans2(ZScore_Heritable, rows = NULL, cols = c(1:2), na.rm = TRUE)
Avg_HighAmmoniaZScore <- rowMeans2(ZScore_Heritable, rows = NULL, cols = c(3:5), na.rm = TRUE)
Avg_HighLactateZScore <- rowMeans2(ZScore_Heritable, rows = NULL, cols = c(6:8), na.rm = TRUE)
Avg_HighOsmolalityZScore <- rowMeans2(ZScore_Heritable, rows = NULL, cols = c(9:11), na.rm = TRUE)
Avg_HighEverythingZScore <- rowMeans2(ZScore_Heritable, rows = NULL, cols = c(12:13), na.rm = TRUE)

ZScore_Heritable_Avg <- cbind(ZScore_Heritable[,0], Avg_ControlZScore, Avg_HighAmmoniaZScore, Avg_HighLactateZScore, Avg_HighOsmolalityZScore, Avg_HighEverythingZScore)
colnames(ZScore_Heritable_Avg) <- c("Control","High Ammonia","High Lactate","High Osmolality,","High Everything")

#######################################################
########### Function to make heatmaps #################
#######################################################

makeHeatmap <- function(goTerm,rowTitle,height,width) {
  filterByGOandZScore <- ZScore_Heritable_Avg[row.names(ZScore_Heritable) %in% c(intersect(genesInTerm(GOdata)[[goTerm]],HeritableList)),];
  
  #view(filterByGOandZScore);
  
  Heatmap(filterByGOandZScore,
          # col = redgreen(75),
          col = colorRamp2(c(-4,-2,0,2,4), 
                           c(
                             rgb(255, 150, 150, maxColorValue = 255),
                             rgb(255, 0, 0, maxColorValue = 255),
                             rgb(0, 0, 0, maxColorValue = 255),
                             rgb(0, 255, 0, maxColorValue = 255),
                             rgb(150, 255, 150, maxColorValue = 255)
                             )
                           ),
          
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          show_row_names = FALSE,
          show_heatmap_legend = FALSE,
          column_order = colnames(filterByGOandZScore),
          
          #row_order = rownames(filterZScoreMetabolic[!(is.na(filterZScoreMetabolic))])
          
          row_title = rowTitle,
          row_title_rot = 0,
          row_title_gp = grid::gpar(fontsize = 15, fontface = "plain"),
          
          # Control font of bottom column titles
          
          heatmap_height = unit(1.3 * height,"cm"),
          heatmap_width = unit(width,"cm"),
          
  )
  
  
}

#
# List of Needed Information for these heatmaps **Order Sensitive when making maps**
#

goTermList <- c("GO:0043065","GO:0043066","GO:0008152","GO:0030968","GO:0033554","GO:0007049","GO:0006457","GO:0016049","GO:0016192","GO:0007010");
goTermNameList <- c("Upregulated \n Apoptosis","Downregulated \n Apoptosis","Metabolism","ER UPR","Stress \n Response","Cell Cycle","Protein Folding","Cell Growth","Vesicular \n Trafficking","Cytoskeleton \n Organization");

######################################################
### Constructing the Heat map drawing and Legend #####
######################################################

htlist = makeHeatmap(goTermList[1],goTermNameList[1],2,10) %v% makeHeatmap(goTermList[2],goTermNameList[2],2,10) %v%
         makeHeatmap(goTermList[3],goTermNameList[3],1.5,10) %v% makeHeatmap(goTermList[4],goTermNameList[4],1,10) %v%
         makeHeatmap(goTermList[5],goTermNameList[5],2,10) %v% makeHeatmap(goTermList[6],goTermNameList[6],1.5,10) %v%
         makeHeatmap(goTermList[7],goTermNameList[7],0.5,10) %v% makeHeatmap(goTermList[8],goTermNameList[8],1,10) %v%
         makeHeatmap(goTermList[9],goTermNameList[9],2,10) %v% makeHeatmap(goTermList[10],goTermNameList[10],4,10)

lgd = Legend(col_fun = colorRamp2(c(-4,-2,0,2,4), 
                                  c(
                                    rgb(255, 150, 150, maxColorValue = 255),
                                    rgb(255, 0, 0, maxColorValue = 255),
                                    rgb(0, 0, 0, maxColorValue = 255),
                                    rgb(0, 255, 0, maxColorValue = 255),
                                    rgb(150, 255, 150, maxColorValue = 255)
                                  )),
             legend_height = unit(4,"cm"),
             legend_width = unit(5,"cm"),
             
             title = "Color Key",
             title_position = "topcenter",
             title_gp = gpar(fontsize = 15, fontface = "plain"),
             
             labels_gp = gpar(fontsize = 12.5, fontface = "plain")
             )

##########################################
##### Drawing the Heat maps and Legend ###
##########################################

draw(htlist);

draw(lgd, x = unit(33, "cm"), y = unit(11.5, "cm"), just = c("left","bottom"));

jpeg("HeritableHeatMap_GOTerms")
grid.newpage()
grid.draw(htlist)
dev.off()