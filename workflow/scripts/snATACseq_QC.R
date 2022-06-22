#--------------------------------------------------------------------------------------
#
#    ArchR - Initial QC
#
#--------------------------------------------------------------------------------------

## need to check GE doublet numbers - may be table issue lines 223-228

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
FRAGS_THRESH <- 3000
TSS_THRESH <- 4 
MAX_CLUSTERS <- 6
VAR_FEATURES <- 25000
N_START <- 10

addArchRThreads(threads = 24) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")
#setwd(OUT_DIR) # Required or saves all files to ~/

# Create ArchR output directry
cat('\nCreate output directory for Arch R project  ... \n')
dir.create(OUT_DIR, recursive = TRUE) # Required ArchR doesn't create this for you

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

# Assign Marker genes for plots
if (REGION == 'FC') {

  MARKER_GENES <-  c('SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'GLI3',
                     'TNC', 'C3', 'SPI1', 'MEF2C')

} else {

  MARKER_GENES <-  c('GAD1', 'GAD2', 'SLC32A1', 'GLI3', 'SLC17A7',
                     'TNC', 'PROX1', 'SCGN', 'LHX6', 'NXPH1',
                     'MEIS2','ZFHX3', 'SPI1', 'LHX8', 'ISL1', 'GBX2')

}

# Load Seurat RNA data for unconstrained integration
cat(paste0('\nLoading Seurat snRNAseq data for ', REGION, ' ... \n'))
if (REGION == 'FC') {

  cat(paste0('\nLoading Seurat object for and region specific variables for', REGION, ' ... \n'))
  seurat.obj <- readRDS("../resources/R_objects/seurat.pfc.final.rds")
  seurat.obj$cellIDs <- gsub('FC-', '', seurat.obj$cellIDs)

} else {

  cat(paste0('\nLoading Seurat object for and region specific variables for', REGION, ' ... \n'))
  seurat.obj <- readRDS("../resources/R_objects/seurat.wge.final.rds")
  seurat.obj$cellIDs <- gsub('GE-', '', seurat.obj$cellIDs)

}

##  Load snATACseq data - Cptr 1.5  ---------------------------------------------------
cat('\nCreating Arrow files ... \n')
ArrowFiles <- createArrowFiles(
  inputFiles = c(paste0(DATA_DIR, SAMPLES[1], "/outs/fragments.tsv.gz"),
                 paste0(DATA_DIR, SAMPLES[2], "/outs/fragments.tsv.gz"),
                 paste0(DATA_DIR, SAMPLES[3], "/outs/fragments.tsv.gz")),
  sampleNames = SAMPLE_IDs,
  minTSS = TSS_THRESH, # Dont set this too high because you can always increase later
  minFrags = FRAGS_THRESH, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = paste0(OUT_DIR, "/QualityControl"),
  
)

##  Doublets  - Cptr 2  ---------------------------------------------------------------
cat('\nCalculating Doublet scores ... \n')
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search 
  # with doublet projection.
  LSIMethod = 1,
  outDir = paste0(OUT_DIR, "/QualityControl")
)


##  Create Arrow project  - Cptr 3  ---------------------------------------------------
cat('\nCreate output directory for Arch R project  ... \n')
dir.create(OUT_DIR, recursive = TRUE) # Required ArchR doesn't create this for you

cat('\nCreating ArchR project ... \n')
archR <- ArchRProject(ArrowFiles = ArrowFiles, 
                           outputDirectory = OUT_DIR,
                           copyArrows = TRUE # This is recommened so that if you modify 
                           # the Arrow files you have an original copy for later usage.
)

##  Save and load Arrow project  - Cptr 3.5  ------------------------------------------
cat('\nSaving ArchR project ... \n')
saveArchRProject(ArchRProj = archR, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)

# Load project
# archR <- loadArchRProject(path = "")

##  Add coldata  ----------------------------------------------------------------------
archR$donor <- word(archR$Sample, 1, sep = "_")


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
  preQC_df$log10nFrags <- log10(preQC_df$nFrags)
  
  # TSS Plot
  preQC_tss_uFrag_plot <- ggPoint(
    x = preQC_df[,"log10nFrags"], 
    y = preQC_df[,"TSSEnrichment"], 
    title = SAMPLES[SAMPLE],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(preQC_df[,"log10nFrags"], probs = 0.99)),
    ylim = c(0, quantile(preQC_df[,"TSSEnrichment"], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
  
  # Assign frag plots and counts_df
  assign(paste0("preQC_tss_uFrag_plot_", donorID), preQC_tss_uFrag_plot)
  assign(paste0("counts_df_", donorID), 
         data.frame("Sample" = sampleID,
                    "Cells_Pass_Filter" = sum(preQC_df$Keep),
                    "Cells_dropped" = sum(preQC_df$Keep == 0),
                    "Total_Frags" = sum(preQC_df$nFrags),
                    "Median_Frags" = median(preQC_df$nFrags[preQC_df$Keep==1]),
                    "Median_TSS_Enrichment" = median(preQC_df$TSSEnrichment[preQC_df$Keep==1])))
  
}

## Initial QC reporting  --------------------------------------------------------------
# Pre-filter tss-frag plot
cat('\nCreating pre-filter plots ... \n')
preQC_tss_uFrag_plot <- ggAlignPlots(preQC_tss_uFrag_plot_510, preQC_tss_uFrag_plot_611, 
                                     preQC_tss_uFrag_plot_993, type = "h")
# Counts df
counts_df <- rbind(counts_df_510, counts_df_611, counts_df_993)

## PostQC
archR.meta <- as.data.frame(getCellColData(archR))
archR.meta$log10nFrags <- log10(archR.meta$nFrags)

tss_uFrag_plot <- ggPoint(
  x = archR.meta[,"log10nFrags"], 
  y = archR.meta[,"TSSEnrichment"], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(archR.meta[,"log10nFrags"], probs = 0.99)),
  ylim = c(0, quantile(archR.meta[,"TSSEnrichment"], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")


fragSize_plot <- plotFragmentSizes(ArchRProj = archR)
fragSize_plot 

tss_plot <- plotTSSEnrichment(ArchRProj = archR)
tss_plot

tss_uFrag_plot

ridge_plot <- plotGroups(
  ArchRProj =  archR, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

## Filter doublets  -------------------------------------------------------------------
cat('\nFiltering doublets ... \n')
archR.2 <- filterDoublets(archR)
doublet_df <- cbind(as.data.frame(table(archR$Sample)), as.data.frame(table(archR.2$Sample)))
doublet_df[3] <- NULL
doublet_df$cells_removed <- 100 - doublet_df[3] / doublet_df[2] * 100
colnames(doublet_df) <- c("Sample", "Pre_DoubRem", "Post_DoubRem", "pc_cells_removed")
doublet_df

##  Dimensionality reduction  ---------------------------------------------------------
cat('\nRunning dimensionality reduction - pre-batch correction ... \n')
archR.2 <- addIterativeLSI(
  ArchRProj = archR.2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = MAX_CLUSTERS
  ), 
  varFeatures = VAR_FEATURES, 
  dimsToUse = 1:30
  
)

##  Clustering  -----------------------------------------------------------------------
cat('\nClustering cells  ... \n')
archR.2 <- addClusters(
  input = archR.2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

##  Visualisation  --------------------------------------------------------------------
cat('\nCreating UMAP ... \n')
archR.2 <- addUMAP(
  ArchRProj = archR.2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

## Clustering - reporting  ------------------------------------------------------------
# Cluster counts - after Iterative LSI based clustering
cat('\nCreating tables and plots for Iterative LSI based clustering ... \n')
clusters_cnts <- as.data.frame(t(as.data.frame(as.vector((table(archR.2$Clusters))))))
rownames(clusters_cnts) <- NULL
colnames(clusters_cnts) <- names(table(archR.2$Clusters))

# Confusion matrix - cell counts per donor
cat('Creating confusion matrix for cell counts per donor ... \n')
cM_LSI <- confusionMatrix(paste0(archR.2$Clusters), paste0(archR.2$Sample))
colnames(cM_LSI) <- colnames(cM_LSI) %>% str_remove("_ATAC")
cM_LSI <- cM_LSI[ gtools::mixedsort(row.names(cM_LSI)), ]
rownames(cM_LSI) <- factor(rownames(cM_LSI),
                           levels = rownames(cM_LSI))

clust_CM_LSI <- pheatmap::pheatmap(
  mat = cM_LSI,
  color = paletteContinuous("whiteBlue"),
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
  treeheight_col = 0,
  treeheight_row = 0,
  angle_col = 0,
  number_color = 'black'
  )
clust_CM_LSI

#Old
#cM_LSI <- confusionMatrix(paste0(archR.2$Clusters), paste0(archR.2$Sample))
#clust_CM_LSI <- pheatmap::pheatmap(
#  mat = as.matrix(cM_LSI), 
#  color = paletteContinuous("whiteBlue"), 
#  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
#  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
#  treeheight_col = 0,
#  treeheight_row = 0,
#  angle_col = 0,
#  number_color = 'black'
#)
#clust_CM_LSI


# Plot UMAP - for Integrated LSI clusters
cat('Create UMAPs ... \n')
clusters_UMAP <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                               name = "Clusters", embedding = "UMAP") +
  NoLegend() + ggtitle('Clusters')
clusters_UMAP_BySample <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                        name = "Sample", embedding = "UMAP") +
    NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')


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
  group_plot <- plot_grid(clusters_UMAP, clusters_UMAP_BySample, plot_stacked_pct,
                          clust_CM_LSI$gtable, ncol = 2, align = 'hv', axis = 'rl')


# Gene specific UMAPs using imputation
# Note that I'm not saving these imputation weights in Save ArchR section below
# They are for visual cell IDing only at this stage
archR.3 <- addImputeWeights(archR.2)

genes_UMAP <- plotEmbedding(
  ArchRProj = archR.3, 
  colorBy = "GeneScoreMatrix", 
  name = MARKER_GENES, 
  embedding = 'UMAP',
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

#  Run unconstrained integration  -----------------------------------------------------
cat('\nRunning unconstrained integration ... \n')
archR.3 <- addGeneIntegrationMatrix(
  ArchRProj = archR.3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seurat.obj,
  addToArrow = FALSE,
  groupRNA = "cellIDs",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

## Unconstrained integration - reporting  ---------------------------------------------
# Confusion matrix - unconstrained cell mappings 
cat('\nCreating confusion matrix ... \n')
cM_geneExp <- as.matrix(confusionMatrix(unname(unlist(getCellColData(archR.3)['Clusters'])), 
                        archR.3$predictedGroup_Un))
clust_CM_geneExp <- pheatmap::pheatmap(
  mat = as.matrix(cM_geneExp), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
  cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
  treeheight_col = 0,
  treeheight_row = 0,
  number_color = 'black'
)


# Get df of top cellID matches from RNA for each ATAC cluster
preClust <- colnames(cM_geneExp)[apply(cM_geneExp, 1 , which.max)]
integration_df <- t(as.data.frame(cbind(preClust, rownames(cM_geneExp)))) #Assignments
rownames(integration_df) <- c("RNA", "ATAC")
colnames(integration_df) <- NULL 


# Plot RNA and ATAC UMAPs for comparison
cat('\nCreating UMAP ... \n')
clusters_UMAP <- plotEmbedding(ArchRProj = archR.3, colorBy = "cellColData", 
                               name = 'Clusters', 
                               embedding = 'UMAP') +
  NoLegend() + ggtitle('Clusters')

# Prepare cell groupings for constrained integration
# Only cell-types in preClust need to be included 
cM_unconstrained <- as.matrix(confusionMatrix(unname(unlist(getCellColData(archR.3)['Clusters'])),
                              archR.3$predictedGroup_Un))
preClust <- colnames(cM_unconstrained)[apply(cM_unconstrained, 1 , which.max)]
cM_unconstrained2 <- cbind(preClust, rownames(cM_unconstrained))
unique(unique(archR.3$predictedGroup_Un))

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
