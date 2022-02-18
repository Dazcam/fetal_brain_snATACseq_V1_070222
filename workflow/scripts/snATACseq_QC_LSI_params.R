#--------------------------------------------------------------------------------------
#
#    ArchR - QC to test LSI parameters for batch correction
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Requirements  ----------------------------------------------------------------------

# Required on Hawk before opening R
# module load libgit2/1.1.0
# module load R/4.0.3
# This script was written for Hawk - need to alter lib paths to run elsewhere

## Info  ------------------------------------------------------------------------------

#  Run on archR project that has been run past basic QC stage 

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(cowplot)
library(rmarkdown)
library(argparser)
library(plyr)
library(gtools)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snATACseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "parameter", help = "No Iterative LSI parameter specified")
p <- add_argument(p, "data_dir", help = "No input data directory specified")
p <- add_argument(p, "archR_out_dir", help = "No ArchR output directory specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
args <- parse_args(p)
print(args)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
REGION <- args$region
PARAMETER <- args$parameter
DATA_DIR <- args$data_dir
OUT_DIR <- args$archR_out_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file

addArchRThreads(threads = 24) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")

# Default QC parameters (for markdown doc)
TSS <- 4
FRAGS <- 2500

# Initial QC parameters
#FRAGS_THRESH <- c(1000, 1500, 2000, 2500, 3000) # Default 2500
#TSS_THRESH <- c(3, 4, 5, 6, 7) # Default 4

# Iterative LSI parameters
ITERATIONS <- c(2, 4, 6, 8, 10) # Default 2
RESOLUTIONS <- c(0.1, 0.2, 0.4, 0.8, 1, 2) # Default 0.2
VAR_FEATURES <- c(10000, 15000, 20000, 25000, 30000) # Default 25000
MAX_CLUSTERS <- c(2, 4, 6, 8, 10) # Default 6
N.STARTS <- c(6, 8, 10, 12) # Default 10
DIMENSIONS <- c(20, 25, 30, 35) # Default 30
SAMPLE_CELLS <- c(5000, 7500, 10000, 12500) # Default 10000

# Set the Variables for LSI
if (PARAMETER == 'Iteration') {
  
  PARAMETERS <- ITERATIONS
  RESOLUTION <- 0.2
  VAR_FEATURE <- 25000
  MAX_CLUSTER <- 6
  N.START <- 10
  DIMENSION <- 30
  SAMPLE_CELL <- 10000

} else if (PARAMETER == 'Resolution') {
    
  ITERATION <- 2
  PARAMETERS <- RESOLUTIONS
  VAR_FEATURE <- 25000
  MAX_CLUSTER <- 6
  N.START <- 10
  DIMENSION <- 30
  SAMPLE_CELL <- 10000
  
} else if (PARAMETER == 'Variable_Features') {
  
  ITERATION <- 2
  RESOLUTION <- 0.2
  PARAMETERS <- VAR_FEATURES 
  MAX_CLUSTER <- 6
  N.START <- 10
  DIMENSION <- 30
  SAMPLE_CELL <- 10000
  
} else if (PARAMETER == 'Max_Clusters') {
  
  ITERATION <- 2
  RESOLUTION <- 0.2
  VAR_FEATURE <- 25000
  PARAMETERS <- MAX_CLUSTERS
  N.START <- 10
  DIMENSION <- 30
  SAMPLE_CELL <- 10000
  
} else if (PARAMETER == 'N_starts') {
  
  ITERATION <- 2
  RESOLUTION <- 0.2
  VAR_FEATURE <- 25000
  MAX_CLUSTER <- 6
  PARAMETERS <- N.STARTS
  DIMENSION <- 30
  SAMPLE_CELL <- 10000
  
} else if (PARAMETER == 'Dimensions') {
  
  ITERATION <- 2
  RESOLUTION <- 0.2
  VAR_FEATURE <- 25000
  MAX_CLUSTER <- 6
  N.START <- 10
  PARAMETERS <- DIMENSIONS
  SAMPLE_CELL <- 10000
  
} else {
  
  ITERATION <- 2
  RESOLUTION <- 0.2
  VAR_FEATURE <- 25000
  MAX_CLUSTER <- 6
  N.START <- 10
  DIMENSION <- 30
  PARAMETERS <- SAMPLE_CELLS
  
}
  
##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)

# Loop to extract sample IDs 
if (REGION == "FC") {
  
  SAMPLES <- c("14510_PFC_ATAC", "14611_PFC_ATAC", "14993_PFC_ATAC")
  SAMPLE_IDs <- SAMPLES %>% str_remove("P") %>% str_remove("14")  
  
  
} else {
  
  SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
  SAMPLE_IDs <- SAMPLES %>% str_remove("W") %>% str_remove("14") 
  
} 

##  Inital ArchR QC -------------------------------------------------------------------
# ArchR does some QC when loading the files in so need to load the pre-QC info
# Pre-filter
cat('\nLoading pre-filtered data ... \n')
for (SAMPLE in 1:length(SAMPLE_IDs)) {
  
  # Subset IDs
  sampleID <- substr(SAMPLE_IDs[SAMPLE], 1, 6)
  donorID <- substr(SAMPLE_IDs[SAMPLE], 1, 3)
  
  # Load Pre-filtered data
  preQC_df <- readRDS(paste0(OUT_DIR, "/QualityControl/", SAMPLE_IDs[SAMPLE], "/", 
                             SAMPLE_IDs[SAMPLE], "-Pre-Filter-Metadata.rds"))
  
  assign(paste0("counts_df_", donorID), 
         data.frame("Sample" = sampleID,
                    "Cells_Pass_Filter" = sum(preQC_df$Keep),
                    "Cells_dropped" = sum(preQC_df$Keep == 0),
                    "Total_Frags" = sum(preQC_df$nFrags),
                    "Median_Frags" = median(preQC_df$nFrags[preQC_df$Keep==1]),
                    "Median_TSS_Enrichment" = median(preQC_df$TSSEnrichment[preQC_df$Keep==1])))
  
}

# Counts df
counts_df <- rbind(counts_df_510, counts_df_611, counts_df_993)


for (i in 1:length(PARAMETERS)) {
  
 # Looop to set parameter to test 
   if (PARAMETER == 'Iteration') {
    
    ITERATION <- PARAMETERS[i]
    
  } else if (PARAMETER == 'Resolution') {
    
    RESOLUTION <- PARAMETERS[i]
    
  } else if (PARAMETER == 'Variable_Features') {
    
    VAR_FEATURE <- PARAMETERS[i]
    
  } else if (PARAMETER == 'Max_Clusters') {
    
    MAX_CLUSTER <- PARAMETERS[i]
    
  } else if (PARAMETER == 'N_starts') {
    
    N.START <- PARAMETERS[i]
    
  } else if (PARAMETER == 'Dimensions') {
    
    DIMENSION <- PARAMETERS[i]
    
  } else {
    
    SAMPLE_CELL <- PARAMETERS[i]
    
  }
  
  cat(paste0('\nRunning LSI changing ', PARAMETER, ' param to ', PARAMETERS[i], ' ... \n'))
  
  # At the moment overwriting each LSI, Cluster and UMAP run
  # Would like to save all of them in ArchR object
  # Issue is subsetting ArchR  object by $ need another method - line 113
  # Currently only an issue with the Clustering name
  LSI_NAME <- paste0("LSI_", PARAMETER, "_", PARAMETERS[i])
  #CLUSTERS_NAME <- paste0("Clusters_iters", PARAM)
  UMAP_NAME <- paste0("UMAP_", PARAMETER, "_", PARAMETERS[i])

  cat('Running LSI ... \n')
  archR.2 <- addIterativeLSI(
    ArchRProj = archR,
    useMatrix = "TileMatrix", 
    name = LSI_NAME, 
    iterations = ITERATION, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = RESOLUTION, 
      sampleCells = SAMPLE_CELL, 
      n.start = N.START
    ), 
    varFeatures = VAR_FEATURE, 
    dimsToUse = 1:DIMENSION
  )
  
  cat('Assigning clusters ... \n')
  archR.2 <- addClusters(
      input = archR.2,
      reducedDims = LSI_NAME,
      method = "Seurat",
      name = "Clusters",
      resolution = 0.8,
      force = TRUE
    )
  
  cat('Creating UMAP ... \n')
  archR.2 <- addUMAP(
    ArchRProj = archR.2, 
    reducedDims = LSI_NAME, 
    name = UMAP_NAME, 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
  )
  
  ##  Integrating snATACseq with snRNAseq - Cptr 8  -----------------------------------------
  # Load Seurat RNA data
  if (REGION == 'FC') {
    
    cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
    seurat.obj <- readRDS("../resources/R_objects/seurat.pfc.final.rds")
    seurat.obj$cellIDs <- gsub('FC-', '', seurat.obj$cellIDs)
    
  } else {
    
    cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
    seurat.obj <- readRDS("../resources/R_objects/seurat.wge.final.rds")
    seurat.obj$cellIDs <- gsub('GE-', '', seurat.obj$cellIDs)
    
  }
  
  #  Run unconstrained integration
  cat('\nRunning unconstrained integration ... \n')
  archR.2 <- addGeneIntegrationMatrix(
    ArchRProj = archR.2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = LSI_NAME,
    seRNA = seurat.obj,
    addToArrow = FALSE,
    groupRNA = "cellIDs",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
  )
  
  ## Reporting
  # UMAPs
  cat('Plotting UMAP by cluster and by sample ... \n')
  clusters_UMAP <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                 name = 'Clusters', embedding = UMAP_NAME) +
    NoLegend() + ggtitle(paste0('Clusters - ', PARAMETER, ': ', PARAMETERS[i]))
  
  clusters_UMAP_BySample <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                          name = "Sample", embedding = UMAP_NAME) +
    NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')
  
  
  # Cell count df - per cluster
  cat('Creating cell count per sample confusion matrix ... \n')
  clusters_cnts <- as.data.frame(t(as.data.frame(as.vector((table(archR$Clusters))))))
  rownames(clusters_cnts) <- NULL
  colnames(clusters_cnts) <- names(table(archR$Clusters))
  
  # Confusion matrices - cell counts per donor
  cat('Creating confusion matrix for cell counts per donor ... \n')
  cM_LSI <- confusionMatrix(paste0(archR$Clusters), paste0(archR$Sample))
  colnames(cM_LSI) <- c("611_FC", "510_FC", "993_FC")
  cM_LSI <- cM_LSI[ mixedsort(row.names(cM_LSI)), ]
  rownames(cM_LSI) <- factor(rownames(cM_LSI), 
                             levels = rownames(cM_LSI))
  clust_CM_LSI <- pheatmap::pheatmap(
    mat = cM_LSI, 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
    cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
    treeheight_col = 0,
    angle_col = 0,
    number_color = 'black'
    )
  clust_CM_LSI
  
  # Confusion matrices - cell assignments after unconstrained integration
  cat('Creating confusion matrix for cell assignments after unconstrained integration ... \n')
  cM_geneExp <- as.matrix(confusionMatrix(archR.2$Clusters, archR.2$predictedGroup_Un))
  cM_geneExp <- cM_geneExp[ mixedsort(row.names(cM_geneExp)), ]
  clust_CM_geneExp <- pheatmap::pheatmap(
    mat = as.matrix(cM_geneExp), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
    cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
    treeheight_col = 0,
    treeheight_row = 0,
    number_color = 'black',
    fontsize_col = 8,
    legend = FALSE
  )

  # Unconstrained integration cluster assignments
  preClust <- colnames(cM_geneExp)[apply(cM_geneExp, 1 , which.max)]
  integration_df <- t(as.data.frame(cbind(preClust, rownames(cM_geneExp)))) #Assignments
  rownames(integration_df) <- c("RNA", "ATAC")
  colnames(integration_df) <- NULL 


  # Stacked barplots
  cat('Creating stacked barplots ... \n')
  cnts_per_donor <- as.data.frame(as.matrix(cM_LSI)) %>%
    rownames_to_column("Cluster") 
  cnts_per_donor$Cluster <- as.factor(cnts_per_donor$Cluster)
  cnts_per_donor_melt <- reshape2::melt(cnts_per_donor, id = 'Cluster') 
  cnts_per_donor_melt$Cluster <- factor(cnts_per_donor_melt$Cluster, 
                                        levels = rownames(cM_LSI))
  
  # Get the levels for type in the required order - https://stackoverflow.com/questions/22231124
  cnts_per_donor_melt$variable = factor(cnts_per_donor_melt$variable, 
                                        levels = c("510_FC",  
                                                   "611_FC", 
                                                   "993_FC"))
  cnts_per_donor_melt = arrange(cnts_per_donor_melt, Cluster, desc(variable))
  
  # Calculate percentages
  cnts_per_donor_melt = plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, percent = value/sum(value) * 100)
  
  # Format the labels and calculate their positions
  cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, pos = (cumsum(value) - 0.5 * value))
  cnts_per_donor_melt$label = paste0(sprintf("%.0f", df$percent), "%")
  
  # Plot - Note this could also be shown with bars filling plot
  plot_stacked_pct <- ggplot(cnts_per_donor_melt, aes(x = factor(Cluster), y = percent, fill = variable)) +
    geom_bar(position = position_stack(), stat = "identity") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2) +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1, fill = NA),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour = "#000000", size = 12)) +
  xlab(NULL) + ylab(NULL)
  
  plot_stacked_cnt <- ggplot(cnts_per_donor_melt, aes(x = factor(Cluster), y = value, fill = variable)) +
    geom_bar(position = position_stack(), stat = "identity") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2) +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1, fill = NA),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour = "#000000", size = 12)) +
    xlab(NULL) + ylab(NULL)
  
  cat('Creating group plot ... \n')
  group_plot <- plot_grid(clusters_UMAP, clusters_UMAP_BySample, plot_stacked_pct, plot_stacked_cnt,
                          clust_CM_LSI$gtable, clust_CM_geneExp$gtable,
            ncol = 2, align = 'hv', axis = 'rl')
  
  assign(paste0(PARAMETER ,'_QC_grp_plt_', PARAMETERS[i]), group_plot)
  assign(paste0(PARAMETER ,'_integration_df', PARAMETERS[i]),integration_df) 

}

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

