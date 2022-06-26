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

# Set GE cluster to remove (expresses SPI1 and SLC17A7)
GE_CLUST_TO_RM <- 'C7'
GE_CLUST_TO_RM_HARMONY <- 'C11'

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

##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)

## Clusters QC ----------------------------------------------------------------------------
for (CLUSTER_TYPE in c('Clusters', 'Clusters_harmony')) {

  if (REGION == 'GE') {
  
    if (CLUSTER_TYPE == 'Clusters') {
    
        CLUSTER_RM <- GE_CLUST_TO_RM  
        
    } else { 
        
        CLUSTER_RM <- GE_CLUST_TO_RM_HARMONY
        
    }
        
    ## Remove of cells expressing SPI1 and SLC17A7 from GE  -------------------------------
    cat('\nObject before removing GE cluster expressing SPI1 and SLC17A7 ...\n')
    print(archR)
    cat('\nRemoving cluster ', CLUSTER_RM, 'from ', CLUSTER_TYPE, 'LSI in the GE ...\n')
    cell_cnt_pre_GE_rm <- length(archR$cellNames)
    SPI1_SLC17_cells  <- BiocGenerics::which(unname(unlist(getCellColData(archR)[CLUSTER_TYPE])) != CLUSTER_RM)
    cells_to_keep <- archR$cellNames[SPI1_SLC17_cells]
    archR <- archR[cells_to_keep, ]

    cat('\nObject after removing GE cluster expressing SPI1 and SLC17A7 ...\n')
    print(archR)
    cell_cnt_post_GE_rm <- length(archR$cellNames)
    
    cat('\nAssign cell cnts before/after removing GE cluster expressing SPI1 and SLC17A7 ')
    assign(paste0(REGION, '_', CLUSTER_TYPE, '_cell_cnt_pre_GE_rm'),
           cell_cnt_pre_GE_rm)
    assign(paste0(REGION, '_', CLUSTER_TYPE, '_cell_cnt_post_GE_rm'),
           cell_cnt_post_GE_rm)

  }

  ## Retain clusters that have >= 30 cells in at least 2 donors  --------------------------
  # As matrix needed to convert sparse matrix to dense matrix
  cat('\nCheck for clusters that do not have >= 30 cells in at least 2 donors ... \n\n')
  donor_cell_cnts <- as.data.frame(as.matrix(confusionMatrix(paste0(unname(unlist(getCellColData(archR)[CLUSTER_TYPE]))) , 
                                                              paste0(archR$Sample)))) %>%
  rownames_to_column(var = 'Cluster')
  print(donor_cell_cnts)

  cluster_list <- vector()

  # Loop to extract cluster IDs for clusters that have >= 30 cells in at least 2 donors 
  cat('\nDo the following clusters have >= 30 cells in at least 2 donors?\n\n')
  for (line in 1:dim(donor_cell_cnts)[1]) {

    donor_A <- donor_cell_cnts[line, 2] >= 30
    donor_B <- donor_cell_cnts[line, 3] >= 30
    donor_C <- donor_cell_cnts[line, 4] >= 30

    # If statement for clusters to retain
    if (donor_A + donor_B + donor_C >= 2) {cluster_list <- c(cluster_list, donor_cell_cnts[line, 1])}

    cat(paste0('Cluster ', donor_cell_cnts[line, 1], ': ', donor_A + donor_B + donor_C >= 2), '\n')

  }

  cat('\nThe clusters to be retained are: ', cluster_list, '\n')

  cell_list <- as.data.frame(getCellColData(archR, select = c("donor", CLUSTER_TYPE)))
  cells_to_keep <- rownames(cell_list %>% filter(get(CLUSTER_TYPE) %in% cluster_list))
  cells_num_deleted <- length(archR$cellNames) - length(cells_to_keep)

  # Remove cells
  cat('\nObject before retaining clusters that have >= 30 cells in at least 2 donors in', REGION,' ...\n')
  print(archR)
   
  cat(paste0('\nRetaining the following clusters for ', REGION, ':', cluster_list))
  cat(paste0('\nCells deleted: ', cells_num_deleted, '\n'))
  archR.2 <- archR[cells_to_keep, ]

  cat('\nObject after retaining clusters that have >= 30 cells in at least 2 donors in', REGION,' ...\n')
  print(archR.2)

    
  ## Re-cluster after cluster removal 1 -------------------------------------------------
  cat('\nRe-clustering cells ... \n')
  #  Dimensionality reduction  
  archR.2 <- addIterativeLSI(
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

  ## Clustering  
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

  ## Batch effect correction - correcting for Sample based batch effects  ---------------
  # Batch correct the LSI reduction using harmony save as new reduction named 'Harmony'
  cat(paste0('\nRunning batch correction for ', REGION, ' ... \n'))
  archR.2 <- addHarmony(
    ArchRProj = archR.2,
    reducedDims = "IterativeLSI_reclust",
    name = 'Harmony_reclust',
    groupBy = "Sample",
    force = TRUE
  )

  # Re-cluster using batch corrected LSI reduction data as input - save as 'Clusters_harmony'
  cat('\nRe-clustering ... \n')
  archR.2 <- addClusters(
    input = archR.2,
    reducedDims = 'Harmony_reclust',
    method = "Seurat",
    name = 'Clusters_harmony_reclust',
    resolution = 0.8,
    force = TRUE
  )

  # Plot UMAP
  archR.2 <- addUMAP(
    ArchRProj = archR.2,
    reducedDims = 'Harmony_reclust',
    name = 'UMAPHarmony_reclust',
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
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

  print(cM_LSI_reClust)

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

  cluster_plot_reclust <- ggAlignPlots(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, type = "h")


  # Stacked barplots
  cat('Creating stacked barplots ... \n')
  cnts_per_donor <- as.data.frame(as.matrix(cM_LSI_reClust)) %>% 
  tibble::rownames_to_column("Cluster")
  cnts_per_donor$Cluster <- as.factor(cnts_per_donor$Cluster)
  cnts_per_donor_melt <- reshape2::melt(cnts_per_donor, id = 'Cluster')
  cnts_per_donor_melt$Cluster <- factor(cnts_per_donor_melt$Cluster,
                                      levels = sort(cnts_per_donor$Cluster))

  # Get the levels for type in the required order - https://stackoverflow.com/questions/22231124
  cnts_per_donor_melt$variable = factor(cnts_per_donor_melt$variable,
                                        levels = LEVELS)
  cnts_per_donor_melt = arrange(cnts_per_donor_melt, Cluster, desc(variable))

  # Calculate percentages
  cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, percent = value/sum(value) * 100)

  # Format the labels and calculate their positions
  cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, pos = (cumsum(value) - 0.5 * value))
  cnts_per_donor_melt$label <- paste0(sprintf("%.0f", cnts_per_donor_melt$percent), "%")

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
  group_plot_reclust <- plot_grid(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, plot_stacked_pct,
                          re_clust_CM_LSI$gtable, ncol = 2, align = 'hv', axis = 'rl')
  
  # Assign grioup plot and cluster cell count table
  cat('Assigning group plots and cluster cell count table ... \n')
  assign(paste0(REGION, '_re_cluster_cnts'), re_cluster_cnts)
  assign(paste0(REGION, '_group_plot_reclust'), group_plot_reclust)
  assign(paste0(REGION, '_clusters_reclust_UMAP'), clusters_reclust_UMAP)

  
                          
  ## Batch effects - reporting  -------------------------------------------------------------
  # Cluster counts - after Iterative LSI based clustering
  cat('\nCreating tables and plots for Iterative LSI based clustering ... \n')
  clusters_cnts_harmony <- as.data.frame(t(as.data.frame(as.vector(table(unname(unlist(getCellColData(archR.2)['Clusters_harmony_reclust'])))))))
  rownames(clusters_cnts_harmony) <- NULL
  colnames(clusters_cnts_harmony) <- names(table(unname(unlist(getCellColData(archR.2)['Clusters_harmony_reclust']))))

  # Confusion matrix - cell counts per donor
  cat('Creating confusion matrix for cell counts per donor ... \n')
  cM_harmony <- confusionMatrix(paste0(unname(unlist(getCellColData(archR.2)['Clusters_harmony_reclust']))),
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
                                     name = 'Clusters_harmony_reclust', 
                                     embedding = 'UMAPHarmony_reclust') +
    Seurat::NoLegend() + ggtitle('Clusters')

  clusters_UMAP_BySample_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                              name = "Sample", embedding = 'UMAPHarmony_reclust') +
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
  group_plot_reclust_harmony <- plot_grid(clusters_UMAP_har, clusters_UMAP_BySample_har, plot_stacked_pct,
                                          clust_cM_harmony$gtable, ncol = 2, align = 'hv', axis = 'rl')


  # Confusion matrix to compare LSI based and batch corrected based clusters
  cM_harmony_compare <- confusionMatrix(paste0(unname(unlist(getCellColData(archR.2)['Clusters_reclust']))),
                                      paste0(unname(unlist(getCellColData(archR.2)['Clusters_harmony_reclust']))))
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
  
  # Assign group plot and cluster cell count table
  cat('Assigning group plots and cluster cell count table ... \n')
  assign(paste0(REGION, '_clusters_cnts_harmony'),  clusters_cnts_harmony)
  assign(paste0(REGION, '_group_plot_reclust_harmony'), group_plot_reclust_harmony)
  assign(paste0(REGION, '_clusters_UMAP_har'), clusters_UMAP_har)
  assign(paste0(REGION, '_clust_cM_harmony_compare'), clust_cM_harmony_compare)
  
  # Assign cell counts
  assign(paste0(REGION, '_', CLUSTER_TYPE, '_cell_cnt_pre_filter'),
         length(archR$cellNames))
  assign(paste0(REGION, '_', CLUSTER_TYPE, '_clusters_retained'), 
         cluster_list)
  assign(paste0(REGION, '_', CLUSTER_TYPE, '_cnts_cells_to_rm'), 
         length(archR$cellNames) - length(cells_to_keep))
  assign(paste0(REGION, '_', CLUSTER_TYPE, '_cell_cnt_post_filter'), 
         length(archR.2$cellNames))
  
  ## In this section we are generating UMAPs and unconstrained integration plots ------
  #  For reclust and harmony reclusts clusters to assess which is best
  
  for (TYPE in c('Clusters_reclust', 'Clusters_harmony_reclust')) {
  
    if (TYPE == 'Clusters_reclust') {
      
      UMAP_ID <- 'UMAP_reclust'
      LSI_ID <- 'IterativeLSI_reclust'      
    
    } else {
    
      UMAP_ID <- 'UMAPHarmony_reclust'
      LSI_ID <- 'Harmony_reclust'      
      
    }

    # Gene specific UMAPs using imputation  -----------------------------------------------
    # Note that I'm not saving these imputation weights in Save ArchR sectiion below
    cat('\nAdd imputation weights for', TYPE, '... \n')
    archR.3 <- addImputeWeights(archR.2)
    
    cat('\nRunning grouped UMAPs for', TYPE, '... \n')
    genes_UMAP <- plotEmbedding(
    ArchRProj = archR.3, 
    colorBy = "GeneScoreMatrix", 
    name = MARKER_GENES, 
    embedding = UMAP_ID,
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
    cat('\nRunning unconstrained integration for', TYPE, '... \n')
    archR.3 <- addGeneIntegrationMatrix(
    ArchRProj = archR.3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = LSI_ID,
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
    cM_geneExp <- as.matrix(confusionMatrix(unname(unlist(getCellColData(archR.3)[TYPE])), 
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

    # Prepare cell groupings for constrained integration
    # Only cell-types in preClust need to be included 
    cM_unconstrained <- as.matrix(confusionMatrix(unname(unlist(getCellColData(archR.3)[TYPE])),
                                archR.3$predictedGroup_Un))
    preClust <- colnames(cM_unconstrained)[apply(cM_unconstrained, 1 , which.max)]
    cM_unconstrained2 <- cbind(preClust, rownames(cM_unconstrained))
    unique(unique(archR.3$predictedGroup_Un))
    
    # Assign group UMAPs and 
    assign(paste0(REGION, '_', TYPE, '_all_genes_UMAP'), all_genes_UMAP) 
    assign(paste0(REGION, '_', TYPE, '_integration_df'), integration_df) 
    assign(paste0(REGION, '_', TYPE, '_clust_CM_geneExp'), clust_CM_geneExp) 
    
    }
                        
}

## Save ArchR project  ----------------------------------------------------------------
#cat('\nSaving project ... \n')
#saveArchRProject(ArchRProj = archR.2, 
#                 outputDirectory = OUT_DIR, 
 #                load = FALSE)

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
