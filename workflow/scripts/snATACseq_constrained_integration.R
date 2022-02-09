#--------------------------------------------------------------------------------------
#
#    ArchR - Constrained integration
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - run contrained integration

#   Here we map the RNA-seq cluster IDs to our ATAC-seq clusters. It is a 2-step process.
#   First unconstrained integration is run which is a broad pass at mapping the RNA-seq 
#   cluster IDs to the ATAC-seq clusters. Then the constrained integration is run, this 
#   time by adding supervised cell groupings as IDed in the unconstrained analysis 
#   i.e. InNs (InN-1, InN-2) and ExNs (ExN-1, ExN-2) etc. are clumped. This enables more 
#   accurate cell mapping between the modalites.

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


##  Integrating snATACseq with snRNAseq - Cptr 8  -----------------------------------------
# Load Seurat RNA data
cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
if (REGION == 'FC') {
  
  cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
  seurat.obj <- readRDS("../resources/R_objects/seurat.pfc.final.rds")
  seurat.obj$cellIDs <- gsub('FC-', '', seurat.obj$cellIDs)
  
} else {
  
  cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
  seurat.obj <- readRDS("../resources/R_objects/seurat.wge.final.rds")
  seurat.obj$cellIDs <- gsub('GE-', '', seurat.obj$cellIDs)
  
}


## Prepare cell groupings for constrained integration  --------------------------------
# Only cell-types in preClust need to be included 
cM_unconstrained <- as.matrix(confusionMatrix(archR$Clusters_harmony, archR$predictedGroup_Un))
preClust <- colnames(cM_unconstrained)[apply(cM_unconstrained, 1 , which.max)]
cM_unconstrained2 <- cbind(preClust, rownames(cM_unconstrained))
cM_unconstrained2
unique(preClust)

# Note tha the unconstrained markdown reports had to be consulted to fill this section in
if (REGION == 'Cer') {
  
  cat(paste0('\nPrepare cell mappings for ', REGION, ' ... \n')) # Can I automate this section?
  
  # Group the unique RNA cluster IDs assigned to ATAC clusters into broad catagories
  cExN <- "ExN-1"
  cInN <- "InN-1"
  cRG <- "RG-2"
  cMG <- "MG"
  cExN_Pro <- "ExN-Pro-1"
  cInN_Pro <- "InN-Pro-2"
  cN_undef <- "N-undef-1"

  
  # Pull out the ATAC cluster IDs that map to RNA clusters in unconstrained run
  clustExN <- rownames(cM_unconstrained)[grep(cExN, preClust)]
  clustInN <- rownames(cM_unconstrained)[grep(cInN, preClust)]
  clustRG <- rownames(cM_unconstrained)[grep(cRG, preClust)]
  clustMG <- rownames(cM_unconstrained)[grep(cMG, preClust)]
  clustExN_Pro <- rownames(cM_unconstrained)[grep(cExN_Pro, preClust)]
  clustInN_Pro <- rownames(cM_unconstrained)[grep(cInN_Pro, preClust)]
  clustN_undef <- rownames(cM_unconstrained)[grep(cN_undef, preClust)]
  
  
  # Get unique cells IDs in these categories from the seurat RNA object
  rnaExN <- colnames(seurat.obj)[grep(cExN, seurat.obj$cellIDs)]
  rnaInN <- colnames(seurat.obj)[grep(cInN, seurat.obj$cellIDs)]
  rnaRG <- colnames(seurat.obj)[grep(cRG, seurat.obj$cellIDs)]
  rnaMG <- colnames(seurat.obj)[grep(cMG, seurat.obj$cellIDs)]
  rnaExN_Pro <- colnames(seurat.obj)[grep(cExN_Pro, seurat.obj$cellIDs)]
  rnaInN_Pro <- colnames(seurat.obj)[grep(cInN_Pro, seurat.obj$cellIDs)]
  rnaN_undef <- colnames(seurat.obj)[grep(cN_undef, seurat.obj$cellIDs)]
  
  # Prepare group list for constrained run
  groupList <- SimpleList(
    ExN = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustExN],
      RNA = rnaExN
    ),
    InN = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustInN],
      RNA = rnaInN
    ),    
    RG = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustRG],
      RNA = rnaRG
    ), 
    MG = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustMG],
      RNA = rnaMG
    ), 
    ExN_Pro = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustExN_Pro],
      RNA = rnaExN_Pro
    ),
    InN_Pro = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustInN_Pro],
      RNA = rnaInN_Pro
    ), 
    N_undef = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustN_undef],
      RNA = rnaN_undef
    ) 
    
  )
  
} else if (REGION == 'FC') {
  
  cat(paste0('\nPrepare cell mappings for ', REGION, ' ... \n')) # Can I automate this section?
  
  # Group the unique RNA cluster IDs assigned to ATAC clusters into broad catagories
  cExN <- paste0(c("ExN-2", "ExN-3", "ExN-4", "ExN-5"), collapse="|")
  cInN <- paste0(c("InN-1", "InN-2", "InN-3"), collapse="|")
  cRG <- paste0(c("RG-1", "RG-2"), collapse="|")
  cMG <- "MG"
  cN_undef <- "N-undef"    

  # Pull out the ATAC cluster IDs that map to RNA clusters in unconstrained run
  clustExN <- rownames(cM_unconstrained)[grep(cExN, preClust)]
  clustInN <- rownames(cM_unconstrained)[grep(cInN, preClust)]
  clustRG <- rownames(cM_unconstrained)[grep(cRG, preClust)]
  clustMG <- rownames(cM_unconstrained)[grep(cMG, preClust)]
  clustN_undef <- rownames(cM_unconstrained)[grep(cN_undef, preClust)]

  # Get unique cells IDs in these categories from the seurat RNA object
  rnaExN <- colnames(seurat.obj)[grep(cExN, seurat.obj$cellIDs)]
  rnaInN <- colnames(seurat.obj)[grep(cInN, seurat.obj$cellIDs)]
  rnaRG <- colnames(seurat.obj)[grep(cRG, seurat.obj$cellIDs)]
  rnaMG <- colnames(seurat.obj)[grep(cMG, seurat.obj$cellIDs)]
  rnaN_undef <- colnames(seurat.obj)[grep(cN_undef, seurat.obj$cellIDs)]

  # Prepare group list for constrained run
  groupList <- SimpleList(
    ExN = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustExN],
      RNA = rnaExN
    ),
    InN = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustInN],
      RNA = rnaInN
    ),    
    RG = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustRG],
      RNA = rnaRG
    ), 
    MG = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustMG],
      RNA = rnaMG
    ),
   N_undef = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustN_undef],
      RNA = rnaN_undef
    )

    
  )
  
} else {
  
  cat(paste0('\nPrepare cell mappings for ', REGION, ' ... \n')) # Can I automate this section?
  
  # Group the unique RNA cluster IDs assigned to ATAC clusters into broad catagories
  cInN <- paste0(c("InN-1", "InN-3", "InN-6", "InN-7"), collapse="|")
  cRG <- paste0(c("RG-1", "RG-2", "RG-3"), collapse="|")
  
  # Pull out the ATAC cluster IDs that map to RNA clusters in unconstrained run
  clustInN <- rownames(cM_unconstrained)[grep(cInN, preClust)]
  clustRG <- rownames(cM_unconstrained)[grep(cRG, preClust)]
  
  # Get unique cells IDs in these categories from the seurat RNA object
  rnaInN <- colnames(seurat.obj)[grep(cInN, seurat.obj$cellIDs)]
  rnaRG <- colnames(seurat.obj)[grep(cRG, seurat.obj$cellIDs)]
  
  # Prepare group list for constrained run
  groupList <- SimpleList(
    InN = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustInN],
      RNA = rnaInN
    ),    
    RG = SimpleList(
      ATAC = archR$cellNames[archR$Clusters_harmony %in% clustRG],
      RNA = rnaRG
    )
    
  )
  
}


# Run constrained integration
cat('\nRun constrained integration ... \n')
archR.2 <- addGeneIntegrationMatrix(
  ArchRProj = archR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seurat.obj,
  addToArrow = FALSE, 
  groupList = groupList,
  groupRNA = "cellIDs",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
  dimsToUse = 1:26,
  k.score = 26,
  npcs = 26,
  dims = 1:26,
  force = TRUE
  #num.cc = 26
)

## Constrained integration - reporting  -----------------------------------------------
# Compare constrained and unconstrained
pal <- paletteDiscrete(values = seurat.obj$cellIDs)
pal

unconstrained_plot <- plotEmbedding(
  archR.2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal,
  
)

constrained_plot <- plotEmbedding(
  archR.2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal,
  
)

integration_compare_plot <- plot_grid(unconstrained_plot, 
                                      constrained_plot, labels = 'AUTO')

# Add integration matrix to ArchR project

#     Add the linked gene expression data to each of the Arrow files.
#     'GroupList' constrains integration and column names to nameCell, 
#     nameGroup, and nameScore for each metadata column we will add 
#     to cellColData. 

cat('\nAdd integration matrix to ArchR project ... \n')
archR.2 <- addGeneIntegrationMatrix(
  ArchRProj = archR.2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seurat.obj,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "cellIDs",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  dimsToUse = 1:26,
  k.score = 26,
  npcs = 26,
  dims = 1:26
)

# Check available matrices
cat('\nCheck avaliable matrices saved in ArchR project ... \n')
getAvailableMatrices(archR.2)

# Add imputed weights 
cat('\nAdd impute weights (is this required?) ... \n')
archR.2 <- addImputeWeights(archR.2)

# Get old ATAC-seq clusters IDs
cM_mappings <- confusionMatrix(archR.2$Clusters_harmony, archR.2$predictedGroup)
labelOld <- rownames(cM_mappings)
labelOld

# Identify the cell type from predictedGroup which best defines that each ATAC cluster
labelNew <- colnames(cM_mappings)[apply(cM_mappings, 1, which.max)]
labelNew
unique(labelNew)
atac_clustIDs <- rbind(labelOld, labelNew)

# Map new labels to clusters
cat('\nMap new labels to clusters ... \n')
archR.2$Clusters_RNAmapped <- mapLabels(archR.2$Clusters_harmony, 
                                           newLabels = labelNew, 
                                           oldLabels = labelOld)

## Final integration - reporting  -----------------------------------------------------
# Confusion matrix to compare old and new cluster mappings
cat('\nCreate plots and tables for report ... \n')
cM_oldNew_compare <- confusionMatrix(paste0(archR.2$Clusters_harmony), 
                                     paste0(archR.2$Clusters_RNAmapped))
clust_CM_oldNew_compare <- pheatmap::pheatmap(
  mat = as.matrix(cM_oldNew_compare), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_oldNew_compare

# Cluster tree to compare old and new cluster mappings
clusttree_OldNew_mappings_df <- as.data.frame(getCellColData(archR.2,
                                                             select = c("Clusters_harmony", 
                                                                        "Clusters_RNAmapped")))
colnames(clusttree_OldNew_mappings_df) <- c("K1", "K2")
clustTree_OldNew_mappings_plot <- clustree(clusttree_OldNew_mappings_df, 
                                           prefix = "K", prop_filter = 0.01) 

# Plot UMAPs
clust_UMAP_oldLabel <- plotEmbedding(archR.2, colorBy = "cellColData", name = "Clusters_harmony")
clust_UMAP_newLabel <- plotEmbedding(archR.2, colorBy = "cellColData", name = "Clusters_RNAmapped")
clust_UMAP_newOldLabel_compare <- plot_grid(clust_UMAP_oldLabel, clust_UMAP_newLabel)

# Donor specific counts after integration
rna_mapped_clusters_df <- as.data.frame(getCellColData(archR.2, select = c("donor", "Clusters_RNAmapped")))
postRNAint_donor_cnts <- rna_mapped_clusters_df %>% dplyr::count(donor, Clusters_RNAmapped) %>%
  pivot_wider(names_from = Clusters_RNAmapped, values_from = n)


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
