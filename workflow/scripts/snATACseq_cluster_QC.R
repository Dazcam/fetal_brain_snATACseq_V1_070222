#--------------------------------------------------------------------------------------
#
#    ArchR - Cluster QC
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - run cluster QC 1

#   Remove clusters that have < 10 cells from 2 out of 3 donors this is to be checked 
#   manually for now but should probs be automated. May need to be run multiple times
#   after re-clustering

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

addArchRThreads(threads = 8) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")


##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)



## Clusters QC ------------------------------------------------------------------------
# Retain clusters that have <= 10 cells from 2 out of 3 donors
# As matrix needed to convert sparse matrix to dense matrix
cat('\nCheck for clusters that have <= 10 cells from 2 out of 3 donors ... \n')
donor_cell_cnts <- as.data.frame(as.matrix(confusionMatrix(paste0(archR$Clusters), paste0(archR$Sample)))) %>%
  rownames_to_column(var = 'Cluster')

cluster_list <- vector()

# Loop to extract cluster IDs for clusters that have > 10 cells from 2 out of 3 donors
for (line in 1:dim(donor_cell_cnts)[1]) {
  
  donor_510 <- donor_cell_cnts[line, 2] <= 10
  donor_611 <- donor_cell_cnts[line, 3] <= 10
  donor_993 <- donor_cell_cnts[line, 4] <= 10
  
  if (donor_510 + donor_611 + donor_993 >= 2) {cluster_list <- c(cluster_list, donor_cell_cnts[line, 1])}
  
}

cell_list <- as.data.frame(getCellColData(archR, select = c("donor", "Clusters")))
cells_to_keep <- rownames(cell_list %>% filter(!Clusters %in% cluster_list))
cells_num_deleted <- length(archR$cellNames) - length(cells_to_keep)

# Remove cells
archR
cat(paste0('\nRemoving the following clusters:', cluster_list))
cat(paste0('Cells deleted: ', cells_num_deleted, '\n'))
archR.2 <- archR[cells_to_keep, ]
archR

## Re-cluster after cluster removal 1 -------------------------------------------------
cat('\nRe-clustering cells ... \n')
#  Dimensionality reduction  
archR.2<- addIterativeLSI(
  ArchRProj = archR.2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI_reclust", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

##  Clustering  
archR.2<- addClusters(
  input = archR.2,
  reducedDims = "IterativeLSI_reclust",
  method = "Seurat",
  name = "Clusters_reclust",
  resolution = 0.8
)

##  Visualisation  
archR.2<- addUMAP(
  ArchRProj = archR.2, 
  reducedDims = "IterativeLSI_reclust", 
  name = "UMAP_reclust", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

## Re-clustering 1 - reporting  -----------------------------------------------------------
# Re-cluster counts - after Iterative LSI based clustering
cat('\nCreating tables and plots ... \n')
re_cluster_cnts <- as.data.frame(t(as.data.frame(as.vector((table(archR.2$Clusters_reclust))))))
rownames(re_cluster_cnts) <- NULL
colnames(re_cluster_cnts) <- names(table(archR.2$Clusters_reclust))

# Confusion matrix - cell counts per donor
cM_LSI_reClust <- confusionMatrix(paste0(archR.2$Clusters_reclust), paste0(archR.2$Sample))
re_clust_CM_LSI <- pheatmap::pheatmap(
  mat = as.matrix(cM_LSI_reClust), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
re_clust_CM_LSI

# Plot UMAP - for Integrated LSI clusters
clusters_reclust_UMAP <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                       name = "Clusters_reclust", embedding = "UMAP_reclust")
clusters_reclust_UMAP_BySample <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                                name = "Sample", embedding = "UMAP_reclust")
cluster_plot <- ggAlignPlots(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, type = "h")


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

