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
FDR_THRESHOLD <- 0.01
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
  newLabel <- c("FC-ExN", "FC-ExN", "FC-undef", "FC-RG", "FC-ExN", "FC-ExN", "FC-ExN", "FC-InN", "FC-RG", "FC-RG", "FC-MG", "FC-InN", "FC-InN", "FC-InN")
  oldLabel <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14")
  archR$Clusters_broad <- mapLabels(archR$Clusters_reclust, newLabels = newLabel, oldLabels = oldLabel)

} else { 

  # Reclassify cell IDs	into broader catagories
  newLabel <- c("LGE-InN", "LGE-InN", "LGE-InN", "MGE-InN", "MGE-InN", "MGE-InN", "GE-Proj", "GE-Proj", "CGE-InN", "GE-RG", "GE-RG", "GE-RG")
  oldLabel <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")
  archR$Clusters_broad <- mapLabels(archR$Clusters, newLabels = newLabel, oldLabels = oldLabel)

}


## Pseudo-bulk replicates - chtr 9  ---------------------------------------------------

#  The underlying assumption in this process is that the single cells 
#  that are being grouped together are sufficiently similar that we do 
#  not care to understand the differences between them. These cell 
#  groupings are almost always derived from individual clusters or 
#  supersets of clusters that correspond to known cell types. 

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

