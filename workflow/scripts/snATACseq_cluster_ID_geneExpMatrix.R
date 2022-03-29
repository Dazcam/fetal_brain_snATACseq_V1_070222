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
library(ggplot2)
library(openxlsx)
library(scales)

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
MARKER_DIR <- "../resources/sheets/"
MARKER_FDR <- c(0.05) # Default 0.01
MARKER_LOG2FC <- c(0.585) # Default 1.25 / 0.585 = 1.5x FC

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
MARKER_GENES <-  c('SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'GLI3',
                   'TNC', 'PROX1', 'SCGN', 'LHX6', 'NXPH1', 
                   'MEIS2','ZFHX3','C3')


##  Gene Scores and Marker Genes  -----------------------------------------------------
cat('Extracting marker genes for each cluster ... \n')
markersGS <- getMarkerFeatures(
  ArchRProj = archR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = clust_ID,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Loop to extract top marker genes at different FDR and Log2FC thresholds
for (LOG2FC_THRESH in c(MARKER_LOG2FC)) {

  for (FDR_THRESH in c(MARKER_FDR)) {
  
    cat(paste0('Creating top marker lists for ', REGION, '. FDR Thresh: ', 
               FDR_THRESH,  'Log2FC Thresh: ', LOG2FC_THRESH, '\n'))
    markerList <- getMarkers(markersGS, cutOff = paste0("FDR <= ", FDR_THRESH, " & Log2FC >= ", LOG2FC_THRESH))
    #gene_markers <- as.data.frame(unlist(markerList))

    # Crete vector of genes of interest
    markerGenes  <- c(MARKER_GENES)
  
    # Create list of heatmaps of top marker genes
    cat('Creating gene expression heatmap ... \n')
    heatmapGS <- plotMarkerHeatmap(
      seMarker = markersGS, 
      cutOff = paste0("FDR <= ", FDR_THRESH, " & Log2FC >= ", LOG2FC_THRESH), 
      labelMarkers = markerGenes,
      transpose = TRUE
    )
  
    # Plot gene expression heatmap
    cat('Plotting gene expression heatmap ... \n')
    geneExp_plot <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", 
                                         annotation_legend_side = "bot")
    dev.off()
    assign(paste0('geneExp_plot_FDR_', FDR_THRESH, '_Log2FC_', LOG2FC_THRESH), geneExp_plot)
  
    # Create gene exp tables
    write.xlsx(as.list(markerList), paste0(REPORT_DIR, REGION, "_FDR_", 
                                           FDR_THRESH, "_Log2FC_", 
                                           LOG2FC_THRESH, ".xlsx"))

  }

}

# Plot UMAP - for Integrated LSI clusters
cat('Create UMAPs ... \n')
clusters_reclust_UMAP <- plotEmbedding(ArchRProj = archR, colorBy = "cellColData", 
                                       name = clust_ID, embedding = UMAP_ID) +
  Seurat::NoLegend() + ggtitle('Clusters')
clusters_reclust_UMAP_BySample <- plotEmbedding(ArchRProj = archR, colorBy = "cellColData", 
                                                name = "Sample", embedding = UMAP_ID) +
  Seurat::NoLegend() + ggtitle('By Donor. R: 510, B: 611, G: 993')

cat('Creating group plot ... \n')
group_plot <- plot_grid(clusters_reclust_UMAP, clusters_reclust_UMAP_BySample, 
                        ncol = 2, align = 'hv', axis = 'rl')

#####  New section for calcuating average expression and pct expression  ################

# Extract gene score matrix and convert to df
gs_mat <- getMatrixFromProject(ArchRProj = archR,
                               useMatrix = 'GeneScoreMatrix')
gs_mat_df <- as.data.frame(as.matrix(assay(gs_mat)))

# Extract gene data and combine with gene score data
gene_data <- rowData(gs_mat)
gene_scores_df <- cbind(gene_data, gs_mat_df)

# Remove gene location info and transpose so genes are cols and cells are rows
gene_scores_df2 <- as.data.frame(t(as.data.frame(gene_scores_df[,7:length(colnames(gene_scores_df))])))
colnames(gene_scores_df2) <- gene_data$name
rownames(gene_scores_df2) <- colnames(gene_scores_df)[7:length(colnames(gene_scores_df))]

# Add cell IDs for cluster and sample data to gene score df
for (i in rownames(gene_scores_df2)) {
  
  gene_scores_df2[i, "Clusters"] <- archR$Clusters_broad[which(archR$cellNames == i)]
  gene_scores_df2[i, "Samples"] <- archR$Sample[which(archR$cellNames == i)]
  
}

##  Generate data for plot  -----------------------------------------------------------
# Calculate percentage of cells expressing genes of interest per cluster
cat('\nCalculating percentage of cells expressing genes of interest per cluster ... ')
for (GENE in all_of(MARKER_GENES)) {
  
  cat('\n\nRunning gene:', GENE, '\n')
  
  for (CLUSTER in unique(gene_scores_df2$Clusters)) {
    
    cat('\nRunning Cluster:', CLUSTER)
    
    # Number of cells (rows) labelled with cluster ID X with a gene score value > 0 for gene Y
    pct_numerator <- length(gene_scores_df2[(gene_scores_df2[GENE] > 0 & gene_scores_df2$Clusters == CLUSTER), GENE])
    
    # Number of rows containing specified cluster ID
    pct_divisor <- nrow(gene_scores_df2[gene_scores_df2$Clusters == CLUSTER,])
    
    # Percentage
    percentage <- pct_numerator / pct_divisor * 100
    
    if (exists('pct_scores')) {
      
      ENTRY <- as.data.frame(cbind(GENE, CLUSTER, percentage))
      pct_scores <- rbind(pct_scores, ENTRY)
      
    } else {
      
      pct_scores <- as.data.frame(cbind(GENE, CLUSTER, percentage))
      
    }
    
  }
  
}


# Calculate average expression expresssing genes of interest per cluster
mean_scores <- gene_scores_df2 %>% 
  select(all_of(MARKER_GENES), Clusters) %>%
  group_by(Clusters) %>%
  summarize(across(everything(), ~mean(.[. > 0]))) %>%
  tidyr::gather(GENE, AVG_EXP, all_of(MARKER_GENES)) %>%
  rename(CLUSTER = Clusters) %>%
  inner_join(pct_scores) %>%
  mutate(PERCENTAGE = as.double(percentage)) %>%
  select(-percentage)

## Generate dot plot
mid <- mean(mean_scores$AVG_EXP)
av_exp_plot <- ggplot(data = mean_scores, mapping = aes_string(x = 'GENE', y = 'CLUSTER')) +
  geom_point(mapping = aes_string(size = 'PERCENTAGE', color = 'AVG_EXP')) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  scale_size_area(limits = c(0, 100), max_size = 7) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red") +
  labs(x = 'gene_name', y = 'Clusters') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0)) +
  coord_flip()

# Create gene exp tables
write_tsv(mean_scores, paste0(REPORT_DIR, REGION, "_geneScore_matrix_avg_exp_and_pct_exp.tsv"))


################################################################################################

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
