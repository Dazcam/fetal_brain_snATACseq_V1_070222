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

#   Remove clusters that do not have >= 30 cells in at least 2 donors. This hs to be checked 
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
library(plyr)
library(gtools)

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

# Loop to extract sample IDs 
if (REGION == "FC") {
  
  # Without 993 AGGR
  #SAMPLES <- c("14510_PFC_ATAC", "14611_PFC_ATAC", "14993_PFC_ATAC")
  #SAMPLE_IDs <- SAMPLES %>% str_remove("P") %>% str_remove("14")  

  # With 993 AGGR
  SAMPLES <- c("14510_PFC_ATAC", "14611_PFC_ATAC", "14993_PFC_ATAC_AGGR")
  SAMPLE_IDs <- SAMPLES %>% str_remove("P") %>% str_remove("14")  %>% str_remove("_AGGR")  

} else if (REGION == "Cer") {

  SAMPLES <- c("14510_Cerebellum_ATAC", "14611_Cerebellum_ATAC", "14993_Cerebellum_ATAC")
  SAMPLE_IDs <- SAMPLES %>% str_remove("ebellum") %>% str_remove("14")  
  
} else {
  
  SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
  SAMPLE_IDs <- SAMPLES %>% str_remove("W") %>% str_remove("14") 
  
} 

LEVELS <- SAMPLE_IDs %>% str_remove("_ATAC") # For stacked barplots
LEVELS

##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)

## Clusters QC ------------------------------------------------------------------------
# Retain clusters that have >= 30 cells in at least 2 donors 
# As matrix needed to convert sparse matrix to dense matrix
cat('\nCheck for clusters that do not have >= 30 cells in at least 2 donors ... \n')
donor_cell_cnts <- as.data.frame(as.matrix(confusionMatrix(paste0(archR$Clusters), paste0(archR$Sample)))) %>%
  rownames_to_column(var = 'Cluster')
donor_cell_cnts

cluster_list <- vector()

# Loop to extract cluster IDs for clusters that have >= 30 cells in at least 2 donors 
cat('Do the following clusters have >= 30 cells in at least 2 donors?')
for (line in 1:dim(donor_cell_cnts)[1]) {
  
  donor_A <- donor_cell_cnts[line, 2] >= 30
  donor_B <- donor_cell_cnts[line, 3] >= 30
  donor_C <- donor_cell_cnts[line, 4] >= 30
  
  # If statement for clusters to retain
  if (donor_A + donor_B + donor_C >= 2) {cluster_list <- c(cluster_list, donor_cell_cnts[line, 1])}
   
  cat(paste0('Cluster ', donor_cell_cnts[line, 1], ': ', donor_A + donor_B + donor_C >= 2), '\n')

}

cat('The clusters to be retained are: ', cluster_list, '\n')

cell_list <- as.data.frame(getCellColData(archR, select = c("donor", "Clusters")))
cells_to_keep <- rownames(cell_list %>% filter(Clusters %in% cluster_list))
cells_num_deleted <- length(archR$cellNames) - length(cells_to_keep)

# Remove cells
archR
cat(paste0('\nRetaining the following clusters:', cluster_list))
cat(paste0('\nCells deleted: ', cells_num_deleted, '\n'))
archR.2 <- archR[cells_to_keep, ]
archR.2

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
cat('Creating confusion matrix for cell counts per donor ... \n')
cM_LSI_reClust <- confusionMatrix(paste0(archR.2$Clusters_reclust), paste0(archR.2$Sample))
colnames(cM_LSI_reClust) <- colnames(cM_LSI_reClust) %>% str_remove("_ATAC")
cM_LSI_reClust <- cM_LSI_reClust[ gtools::mixedsort(row.names(cM_LSI_reClust)), ]
rownames(cM_LSI_reClust) <- factor(rownames(cM_LSI_reClust),
                           levels = rownames(cM_LSI_reClust))

cM_LSI_reClust

re_clust_CM_LSI <- pheatmap::pheatmap(
  mat = as.matrix(cM_LSI_reClust), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
  treeheight_col = 0,
  treeheight_row = 0,
  angle_col = 0,
  number_color = 'black'

)

# Plot UMAP - for Integrated LSI clusters
cat('Create UMAPs ... \n')
clusters_reclust_UMAP <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                       name = "Clusters_reclust", embedding = "UMAP_reclust") +
  NoLegend() + ggtitle('Clusters')
clusters_reclust_UMAP_BySample <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                                name = "Sample", embedding = "UMAP_reclust") +
    NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')

cluster_plot <- ggAlignPlots(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, type = "h")


# Stacked barplots
cat('Creating stacked barplots ... \n')
cnts_per_donor <- as.data.frame(as.matrix(cM_LSI_reClust)) %>% 
  tibble::rownames_to_column("Cluster")
cat('\nA\n')
cnts_per_donor
cnts_per_donor$Cluster <- as.factor(cnts_per_donor$Cluster)
cnts_per_donor_melt <- reshape2::melt(cnts_per_donor, id = 'Cluster')
cnts_per_donor_melt$Cluster <- factor(cnts_per_donor_melt$Cluster,
                                      levels = sort(cnts_per_donor$Cluster))


cat('\nB\n')
cnts_per_donor_melt

# Get the levels for type in the required order - https://stackoverflow.com/questions/22231124
cnts_per_donor_melt$variable = factor(cnts_per_donor_melt$variable,
                                        levels = LEVELS)
cat('\nB2\n')
cnts_per_donor_melt

cnts_per_donor_melt = arrange(cnts_per_donor_melt, Cluster, desc(variable))

cat('\nC\n')
cnts_per_donor_melt

# Calculate percentages
cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, percent = value/sum(value) * 100)

cat('\nD\n')
cnts_per_donor_melt

# Format the labels and calculate their positions
cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, pos = (cumsum(value) - 0.5 * value))
cnts_per_donor_melt$label <- paste0(sprintf("%.0f", cnts_per_donor_melt$percent), "%")

cat('\nE\n')
cnts_per_donor_melt

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

cat('Creating group plot ... \n')
  group_plot <- plot_grid(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, plot_stacked_pct,
                          re_clust_CM_LSI$gtable, ncol = 2, align = 'hv', axis = 'rl')

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

