#--------------------------------------------------------------------------------------
#
#    ArchR - Batch correction
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

#  Run Harmony batch correction on scATACseq data 

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


## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snATACseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
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
DATA_DIR <- args$data_dir
OUT_DIR <- args$archR_out_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file

addArchRThreads(threads = 8) # Set Hawk to 12 cores so 0.75 of total
addArchRGenome("hg38")
#setwd(OUT_DIR) # Required or saves all files to ~/


##  Load ArchR project  -------------------------------------------------------------------
cat('\nLoading ArchR project ... \n')
archR <- loadArchRProject(path = OUT_DIR)

## Batch effect correction - correcting for Sample based batch effects  ---------------
# Batch correct the LSI reduction using harmony save as new reduction named 'Harmony'
cat(paste0('\nRunning batch correction for ', REGION, ' ... \n'))
archR.2 <- addHarmony(
  ArchRProj = archR,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Re-cluster using batch corrected LSI reduction data as input - save as 'Clusters_harmony'
cat('\nRe-clustering ... \n')
archR.2 <- addClusters(
  input = archR.2,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_harmony",
  resolution = 0.8,
  force = TRUE
)

# Plot UMAP
archR.2 <- addUMAP(
  ArchRProj = archR.2,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)


## Batch effects - reporting  -------------------------------------------------------------
# Cluster counts - after batch corrected Iterative LSI based clustering
cat('\nCreate tables and confusion matrices ... \n')
clusters_cnts_harmony <- as.data.frame(t(as.data.frame(as.vector((table(archR.2$Clusters_harmony))))))
rownames(clusters_cnts_harmony) <- NULL
colnames(clusters_cnts_harmony) <- names(table(archR.2$Clusters_harmony))

# Confusion matrix - cell counts per donor
cM_harmony <- confusionMatrix(paste0(archR.2$Clusters_harmony), 
                              paste0(archR.2$Sample))
clust_CM_harmony <- pheatmap::pheatmap(
  mat = as.matrix(cM_harmony), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_harmony

# Plot UMAPs
cat('\nPlotting UMAPs ... \n')
clusters_UMAP_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                   name = "Clusters_harmony", embedding = "UMAPHarmony")
clusters_UMAP_BySample_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                            name = "Sample", embedding = "UMAPHarmony")
cluster_plot_har <- ggAlignPlots(clusters_UMAP_har, clusters_UMAP_BySample_har, type = "h")

# Confusion matrix to compare LSI based and batch corrected based clusters
cM_harmony_compare <- confusionMatrix(paste0(archR.2$Clusters),
                                      paste0(archR.2$Clusters_harmony))
clust_CM_harmony_compare <- pheatmap::pheatmap(
  mat = as.matrix(cM_harmony_compare),
  color = paletteContinuous("whiteBlue"),
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_harmony_compare

# Cluster tree to compare LSI based and batch corrected based clusters
clusttree_harmony_df <- as.data.frame(getCellColData(archR.2,
                                                     select = c("Clusters",
                                                                "Clusters_harmony")))
colnames(clusttree_harmony_df) <- c("K1", "K2")
clustTree_harmony_plot <- clustree(clusttree_harmony_df, prefix = "K", prop_filter = 0.01)

## Save ArchR project  ----------------------------------------------------------------
cat('\nSaving project ... \n')
saveArchRProject(ArchRProj = archR.2, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)


## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
