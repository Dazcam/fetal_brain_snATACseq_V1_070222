#--------------------------------------------------------------------------------------
#
#    ArchR - Label cells using marker genes on Gene score matrix
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#   Use this script instead of integration scripts 
##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(clustree)
library(cowplot)
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
MARKER_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"

addArchRThreads(threads = 8) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")


##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)


##  If statement to specify cluster ID to use  --------------------------------------------
# Load Seurat RNA data
if (REGION == 'FC') {
  
  clust_ID <- 'Clusters_reclust' # Due to cluster QC being run on FC
  UMAP_ID<- 'UMAP_reclust'
  
} else {
  
  clust_ID <- 'Clusters'
  UMAP_ID<- 'UMAP'
  
}

# Collate marker genes for Heatmap
MARKER_GENES <- readxl::read_excel(paste0(MARKER_DIR, "snATACseq_cluster_ID_markers.xlsx")) %>%
  unlist(use.names=FALSE) 
MARKER_GENES <- unique(MARKER_GENES[!is.na(MARKER_GENES)])


##  Gene Scores and Marker Genes  -----------------------------------------------------
markersGS <- getMarkerFeatures(
  ArchRProj = archR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = clustID,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Create df of top markers
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
gene_markers <- as.data.frame(unlist(markerList))

# View genes in marker list
markerList[[1]][["name"]]

# Crete vector of genes of interest
markerGenes  <- c(MARKER_GENES)

# Create list of heatmaps of top marker genes
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

# Plot gene expression heatmap
cat('Create gene expression heatmap ... \n')
geneExp_plot <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", 
                                     annotation_legend_side = "bot")
dev.off()

# Plot UMAP - for Integrated LSI clusters
cat('Create UMAPs ... \n')
clusters_reclust_UMAP <- plotEmbedding(ArchRProj = archR, colorBy = "cellColData", 
                                       name = clust_ID, embedding = UMAP_ID) +
  NoLegend() + ggtitle('Clusters')
clusters_reclust_UMAP_BySample <- plotEmbedding(ArchRProj = archR, colorBy = "cellColData", 
                                                name = "Sample", embedding = UMAP_ID) +
  NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')

cluster_plot <- ggAlignPlots(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, type = "h")

cat('Creating group plot ... \n')
group_plot <- plot_grid(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, 
                        ncol = 2, align = 'hv', axis = 'rl')

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
