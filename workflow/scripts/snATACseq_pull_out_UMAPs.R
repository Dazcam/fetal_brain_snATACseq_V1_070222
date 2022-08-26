#--------------------------------------------------------------------------------------
#
#    ArchR - Pull out UMAPs
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

## Info  ------------------------------------------------------------------------------

#  Run the analysis up until cluster QC cell removal 

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(clustree)
library(cowplot)
library(rmarkdown)
library(argparser)
library(plyr)
library(gtools)

## Parse region / set region variable -------------------------------------------------
for (REGION in c('FC', 'GE')) {

  # Assign Marker genes for plots
  if (REGION == 'FC') {
    
    MARKER_GENES <-  c('SLC17A7', 'GAD1', 'GLI3', 'SPI1')
    
  } else {
    
    MARKER_GENES <-  c('GAD2', 'LHX6', 'SCGN', 'ZFHX3', 'GLI3')
    
  }
  
  ##  Load ArchR project  -------------------------------------------------------------------
  archR <- loadArchRProject(paste0('../results/ARCHR/', REGION))
  
  # Plot cluster UMAP - for Integrated LSI clusters
  clusters_UMAP <- plotEmbedding(ArchRProj = archR, colorBy = "cellColData", 
                                 name = "Clusters_broad", embedding = "UMAP") 
  
  # Create gene specific UMAPs using imputation
  archR.2 <- addImputeWeights(archR)
  genes_UMAP <- plotEmbedding(
    ArchRProj = archR, 
    colorBy = "GeneScoreMatrix", 
    name = MARKER_GENES, 
    embedding = 'UMAP',
    imputeWeights = getImputeWeights(archR.2)
  )
  
  # Save RDS objects
  saveRDS(genes_UMAP, paste0('../results/archR_data_processing/rds_files/', REGION, '_UMAP_genes.rds'))
  
  saveRDS(clusters_UMAP, paste0('../results/archR_data_processing/rds_files/', REGION, '_UMAP_clusters_broad.rds'))

}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------