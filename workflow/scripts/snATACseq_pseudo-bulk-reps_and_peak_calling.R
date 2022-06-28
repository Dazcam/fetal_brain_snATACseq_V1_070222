#--------------------------------------------------------------------------------------
#
#    ArchR - Create psuedo-bulk replicates and peak calling
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - reate psuedo-bulk replicates and peak calling

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
library(Seurat)
library(plyr)
library(gtools)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snATACseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "data_dir", help = "No input data directory specified")
p <- add_argument(p, "archR_out_dir", help = "No ArchR output directory specified")
p <- add_argument(p, "peaks_dir", help = "No peaks output directory specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
p <- add_argument(p, "macs2_path", help = "No path to macs2 binary specified")
args <- parse_args(p)
print(args)


##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
REGION <- args$region
DATA_DIR <- args$data_dir
OUT_DIR <- args$archR_out_dir
PEAKS_DIR <- args$peaks_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file
MACS2_PATH <- args$macs2_path
FDR_THRESHOLD <- 0.05
EXTEND_BPs <- 250


addArchRThreads(threads = 8) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")


##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)


## Set broad catagories for cell IDs  -------------------------------------------------
cat(paste0('\nAssign cell IDs to clusters for ', REGION, ' ... \n'))
if (REGION == 'FC') {

  # Reclassify cell IDs into broader catagories
  newLabel <- c("FC-ExN", "FC-ExN", "FC-ExN", "FC-InN", "FC-InN", 
                "FC-InN", "FC-ExN", "FC-Undef", "FC-InN", "FC-RG", 
                "FC-RG", "FC-RG", "FC-MG", "FC-RG", "FC-RG", 
                "FC-RG")
  oldLabel <- c("C1", "C2", "C3", "C4", "C5", 
                "C6", "C7", "C8", "C9", "C10", 
                "C11", "C12", "C13", "C14", "C15", 
                "C16")
  archR$Clusters_broad <- mapLabels(archR$Clusters, newLabels = newLabel, oldLabels = oldLabel)

} else { 

  # Reclassify cell IDs	into broader catagories
  newLabel <- c("LGE-InN", "GE-RG", "GE-RG", "CGE-InN", "MGE-InN", 
                "MGE-InN", "GE-Undef", "LGE-InN", "CGE-InN", "MGE-InN", 
                "LGE-InN")
  oldLabel <- c("C1", "C2", "C3", "C4", "C5", 
                "C6", "C7", "C8", "C9", "C10", 
                "C11")
  archR$Clusters_broad <- mapLabels(archR$Clusters, newLabels = newLabel, oldLabels = oldLabel)

}

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


## Pseudo-bulk replicates - chtr 9  ---------------------------------------------------
cat(paste0('\nCreate pseudo-bulk replicates for ', REGION, ' ... \n'))
archR.2 <- addGroupCoverages(ArchRProj = archR, groupBy = "Clusters_broad", force = TRUE)

##  Peak Calling  - cptr 10  ----------------------------------------------------------
# Set macs2 path - note that you need to set the default python env to python 2
pathToMacs2 <- MACS2_PATH

# Call peaks
cat(paste0('\nCalling peaks for ', REGION, ' ... \n'))
archR.2 <- addReproduciblePeakSet(
  ArchRProj = archR.2, 
  groupBy = "Clusters_broad", 
  pathToMacs2 = MACS2_PATH,
  cutOff = FDR_THRESHOLD, 
  extendSummits = EXTEND_BPs)

## Peak Calling - Reporting  ----------------------------------------------------------

# Create peak cnt table and bed files for LDSC 
# Print peak calling parameters
# Had to cobble code from ArchR repo to generate this - this is printed to screen
# During peak calling but only partially reproduced in log

cat(paste0('\nCreate tables and plots for report ', REGION, ' ... \n'))
coverageParams <- archR.2@projectMetadata$GroupCoverages[["Clusters_broad"]]$Params
coverage_metadata <- archR.2@projectMetadata$GroupCoverages[["Clusters_broad"]]$coverageMetadata
maxPeaks_default <- 150000
peaksPerCell_default <- 500

tableGroups <- table(getCellColData(archR.2, "Clusters_broad", drop = TRUE))
peakCallParams_summary_df <- lapply(seq_along(coverageParams$cellGroups), function(y){
  x <- coverageParams$cellGroups[[y]]
  uniq <- unique(unlist(x))
  n <- lapply(x, length) %>% unlist %>% sum
  nmin <- lapply(x, length) %>% unlist %>% min
  nmax <- lapply(x, length) %>% unlist %>% max
  data.frame(
    Group=names(coverageParams$cellGroups)[y], 
    nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
    nCellsUsed=length(uniq), 
    nReplicates=length(x), 
    nMin=nmin, 
    nMax=nmax, 
    maxPeaks = min(maxPeaks_default, length(uniq) * peaksPerCell_default)
  )
}) %>% Reduce("rbind",.)

# Plot peak call summary
peak_call_summary <- metadata(archR.2@peakSet)$PeakCallSummary
peak_call_summary_plot <- ggplot(peak_call_summary, 
                                 aes(fill=Var1, y=Freq, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  viridis::scale_fill_viridis(discrete = T) +
  ggtitle(paste0("Peak call summary for ", REGION, ' at FDR < ',
                 FDR_THRESHOLD, ' with ', EXTEND_BPs, 'bp extension')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Annotation")) +
  xlab("") +
  ylab(expression("No. of Peaks"~(x10^{"3"})))

## Create bed files for LDSC  ---------------------------------------------------------
# ArchR peak calling output subs in dots for dashes in cell-type part of rds filename 
dir.create('../results/archR_data_processing/peaks/')
cell_types <- gsub("\\-", "\\.", unique(archR$Clusters_broad)) 



cat(paste0('\nCreating bed files for ', REGION, ' ... \n'))
for (CELL_TYPE in cell_types) {
  
  # Load reproducible peak set for each cell-type
  PEAKS <- readRDS(paste0(OUT_DIR, '/PeakCalls/', 
                          CELL_TYPE, '-reproduciblePeaks.gr.rds'))
  
  # Convert to bed file
  PEAKS_DF <- data.frame(seqnames=seqnames(PEAKS),
                         starts=start(PEAKS)-1,
                         ends=end(PEAKS),
                         names=c(rep(".", length(PEAKS))),
                         scores=c(rep(".", length(PEAKS))),
                         strands=strand(PEAKS))
  
  # Write df to bed file - https://www.biostars.org/p/89341/ 
  write.table(PEAKS_DF, 
              file=paste0(PEAKS_DIR, CELL_TYPE, '.hg38.bed'),
              quote=F, 
              sep="\t", 
              row.names=F, 
              col.names=F)
  
  # Assign Granges object 
  assign(paste0(CELL_TYPE, '_peaks'), PEAKS)
  
}

## Create peak count table
cat(paste0('\nCreate peak count table for ', REGION, ' ... \n'))
for (CELL_TYPE in cell_types) {
  
  if (exists("PEAK_CNT_DF")) {
    
    PEAKS <- get(paste0(CELL_TYPE, '_peaks'))
    PEAK_CNT <- cbind(CELL_TYPE, dim(Repitools::annoGR2DF(PEAKS))[1])
    PEAK_CNT_DF <- rbind(PEAK_CNT_DF, PEAK_CNT)
    
  } else {
    
    PEAKS <- get(paste0(CELL_TYPE, '_peaks'))
    PEAK_CNT_DF <- cbind(CELL_TYPE, dim(Repitools::annoGR2DF(PEAKS))[1])
    colnames(PEAK_CNT_DF) <- c("Cell Type", "Peak count")
    
  }
  
  
}

cat('\nPeak counts are: ... \n')
PEAK_CNT_DF <- base::as.data.frame(PEAK_CNT_DF)
PEAK_CNT_DF

## Create UMAP of broad clusters
cat(paste0('\nCreate broad cluster UMAP for ', REGION, ' ... \n'))
UMAP_broad <- plotEmbedding(archR.2, colorBy = "cellColData", name = "Clusters_broad") +
  ggtitle('Clusters')  

## Clustering - reporting  ------------------------------------------------------------
# Cluster counts - after Iterative LSI based clustering
cat('\nCreating tables and plots for Iterative LSI based clustering ... \n')
clusters_cnts <- as.data.frame(t(as.data.frame(as.vector((table(archR.2$Clusters_broad))))))
rownames(clusters_cnts) <- NULL
colnames(clusters_cnts) <- names(table(archR.2$Clusters_broad))

# Confusion matrix - cell counts per donor
cat('Creating confusion matrix for cell counts per donor ... \n')
cM_LSI <- confusionMatrix(paste0(archR.2$Clusters_broad), paste0(archR.2$Sample))
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

# Plot UMAP - for Integrated LSI clusters
cat('Create UMAPs ... \n')
clusters_UMAP <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                               name = "Clusters_broad", embedding = "UMAP") +
  ggtitle('Clusters')
clusters_UMAP_BySample <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                        name = "Sample", embedding = "UMAP") +
  ggtitle('By Donor. R: 510, B: 611, G: 993')


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



## Save ArchR project  ----------------------------------------------------------------
cat('\n\nSaving project ... \n')
saveArchRProject(ArchRProj = archR.2, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)

getAvailableMatrices(archR.2)

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

