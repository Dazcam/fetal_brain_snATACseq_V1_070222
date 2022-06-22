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
library(cowplot)
library(rmarkdown)
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
p <- add_argument(p, "pre_or_post_clust_QC", help = "State if running pre- or post- clust QC")
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
CLUST_QC_STATUS <- args$pre_or_post_clust_QC

addArchRThreads(threads = 8) # Set Hawk to 12 cores so 0.75 of total
addArchRGenome("hg38")
#setwd(OUT_DIR) # Required or saves all files to ~/

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


##  Load ArchR project  -------------------------------------------------------------------
cat('\nLoading ArchR project ... \n')
archR <- loadArchRProject(path = OUT_DIR)

##  Specify if scripts is run pre- or post clustQC  
if (CLUST_QC_STATUS == 'PRE')	{
  
  LSI_ID <- 'IterativeLSI'
  clust_ID <- 'Clusters' # Due to cluster QC being run on FC
  UMAP_ID<- 'UMAP_reclust'
  HARMONY_ID <- 'Harmony'  
  CLUSTERS_BATCH_CORRECTED_ID <- 'Clusters_harmony' 
  UMAP_BATCH_CORRECTED_ID <- 'UMAPHarmony'  

} else {
  
  LSI_ID <- 'IterativeLSI_reclust'
  clust_ID <- 'Clusters_reclust'
  UMAP_ID<- 'UMAP_reclust'
  HARMONY_ID <- 'Harmony_reclust'
  CLUSTERS_BATCH_CORRECTED_ID <- 'Clusters_harmony_reclust'  
  UMAP_BATCH_CORRECTED_ID <- 'UMAPHarmony_reclust'

}

# Assign Marker genes for plots
if (REGION == 'FC') {

  MARKER_GENES <-  c('SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'GLI3',
                     'TNC', 'C3', 'SPI1', 'MEF2C')

} else {

  MARKER_GENES <-  c('GAD1', 'GAD2', 'SLC32A1', 'GLI3', 'SLC17A7',
                     'TNC', 'PROX1', 'SCGN', 'LHX6', 'NXPH1',
                     'MEIS2','ZFHX3', 'SPI1', 'LHX8', 'ISL1', 'GBX2')

}

## Batch effect correction - correcting for Sample based batch effects  ---------------
# Batch correct the LSI reduction using harmony save as new reduction named 'Harmony'
cat(paste0('\nRunning batch correction for ', REGION, ' ... \n'))
archR.2 <- addHarmony(
  ArchRProj = archR,
  reducedDims = LSI_ID,
  name = HARMONY_ID,
  groupBy = "Sample",
  force = TRUE
)

# Re-cluster using batch corrected LSI reduction data as input - save as 'Clusters_harmony'
cat('\nRe-clustering ... \n')
archR.2 <- addClusters(
  input = archR.2,
  reducedDims = HARMONY_ID,
  method = "Seurat",
  name = CLUSTERS_BATCH_CORRECTED_ID,
  resolution = 0.8,
  force = TRUE
)

# Plot UMAP
archR.2 <- addUMAP(
  ArchRProj = archR.2,
  reducedDims = HARMONY_ID,
  name = UMAP_BATCH_CORRECTED_ID,
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)


## Batch effects - reporting  -------------------------------------------------------------
# Cluster counts - after Iterative LSI based clustering
cat('\nCreating tables and plots for Iterative LSI based clustering ... \n')
clusters_cnts_harmony <- as.data.frame(t(as.data.frame(as.vector(table(unname(unlist(getCellColData(archR.2)[CLUSTERS_BATCH_CORRECTED_ID])))))))
rownames(clusters_cnts_harmony) <- NULL
colnames(clusters_cnts_harmony) <- names(table(unname(unlist(getCellColData(archR.2)[CLUSTERS_BATCH_CORRECTED_ID]))))

# Confusion matrix - cell counts per donor
cat('Creating confusion matrix for cell counts per donor ... \n')
cM_harmony <- confusionMatrix(paste0(unname(unlist(getCellColData(archR.2)[CLUSTERS_BATCH_CORRECTED_ID]))),
                              paste0(archR.2$Sample))
colnames(cM_harmony) <- colnames(cM_harmony) %>% str_remove("_ATAC")
cM_harmony <- cM_harmony[ gtools::mixedsort(row.names(cM_harmony)), ]
rownames(cM_harmony) <- factor(rownames(cM_harmony),
                           levels = rownames(cM_harmony))

clust_cM_harmony <- pheatmap::pheatmap(
  mat = cM_harmony,
  color = paletteContinuous("whiteBlue"),
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
  treeheight_col = 0,
  treeheight_row = 0,
  angle_col = 0,
  number_color = 'black'
  )
print(clust_cM_harmony)


# Plot UMAPs
cat('\nPlotting UMAPs ... \n')
clusters_UMAP_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                   name = CLUSTERS_BATCH_CORRECTED_ID, 
                                   embedding = UMAP_BATCH_CORRECTED_ID) +
  Seurat::NoLegend() + ggtitle('Clusters')

clusters_UMAP_BySample_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                            name = "Sample", embedding = UMAP_BATCH_CORRECTED_ID) +
  Seurat::NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')

# Stacked barplots
cat('Creating stacked barplots ... \n')
cnts_per_donor <- as.data.frame(as.matrix(cM_harmony)) %>%
  rownames_to_column("Cluster")
cnts_per_donor$Cluster <- as.factor(cnts_per_donor$Cluster)
cnts_per_donor_melt <- reshape2::melt(cnts_per_donor, id = 'Cluster')
cnts_per_donor_melt$Cluster <- factor(cnts_per_donor_melt$Cluster,
                                      levels = rownames(cM_harmony))

# Get the levels for type in the required order - https://stackoverflow.com/questions/22231124
cnts_per_donor_melt$variable = factor(cnts_per_donor_melt$variable,
                                        levels = LEVELS)
cnts_per_donor_melt = arrange(cnts_per_donor_melt, Cluster, desc(variable))

# Calculate percentages
cnts_per_donor_melt = plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, percent = value/sum(value) * 100)

# Format the labels and calculate their positions
cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, pos = (cumsum(value) - 0.5 * value))
cnts_per_donor_melt$label = paste0(sprintf("%.0f", cnts_per_donor_melt$percent), "%")

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
  group_plot <- plot_grid(clusters_UMAP_har, clusters_UMAP_BySample_har, plot_stacked_pct,
                          clust_cM_harmony$gtable, ncol = 2, align = 'hv', axis = 'rl')


# Confusion matrix to compare LSI based and batch corrected based clusters
cM_harmony_compare <- confusionMatrix(paste0(unname(unlist(getCellColData(archR.2)[clust_ID]))),
                                      paste0(unname(unlist(getCellColData(archR.2)[CLUSTERS_BATCH_CORRECTED_ID]))))
clust_cM_harmony_compare <- pheatmap::pheatmap(
  mat = as.matrix(cM_harmony_compare),
  color = paletteContinuous("whiteBlue"),
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
  treeheight_col = 0,
  treeheight_row = 0,
  angle_col = 0,
  number_color = 'black'
)
print(clust_cM_harmony_compare)

# Gene specific UMAPs using imputation
# Note that I'm not saving these imputation weights in Save ArchR sectiion below
archR.3 <- addImputeWeights(archR.2)

genes_UMAP <- plotEmbedding(
  ArchRProj = archR.3, 
  colorBy = "GeneScoreMatrix", 
  name = MARKER_GENES, 
  embedding = UMAP_BATCH_CORRECTED_ID,
  imputeWeights = getImputeWeights(archR.3)
)

all_genes_UMAP <- lapply(genes_UMAP, function(x){
  
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})


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
