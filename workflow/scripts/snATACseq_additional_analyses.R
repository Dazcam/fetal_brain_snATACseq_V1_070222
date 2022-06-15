#--------------------------------------------------------------------------------------
#
#    ArchR - Motif and peak to gene linkage analysis analysis
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - Running additonal ArchR analyses on peaks

##  Load Packages  --------------------------------------------------------------------
source('envs/archR_env.R')

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead markdown and report file and directory info ... \n")
#p <- add_argument(p, "archR_out_dir", help = "No ArchR output directory specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "rds_dir", help = "No RDS output directory specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
p <- add_argument(p, "hars_bed", help = "No human accelerated region bed file specified")
args <- parse_args(p)
print(args)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
OUT_DIR <- args$archR_out_dir
MARKDOWN_FILE <- args$markdown_file
RDS_DIR <- args$rds_dir
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file
HARS_FILE <- args$hars_bed
addArchRThreads(threads = 12) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")
FC_motifs <- c('NEUROD2', 'NEUROD4', 'NEUROD6', 'ASCL1', 'SPI1', 'PAX6')
GE_motifs <- c('DLX5', 'ASCL1', 'FOXP1', 'ISL1', 'MEIS2', 'SIX3', 'NFIX', 'WT1')

## Create directories  ------------------------------------------------------------------
cat('\nCreating directories that snakemake is not tracking ... \n')
dir.create(RDS_DIR)

## Main loop for analyses
for (REGION in c("FC", "GE")) {
  
  if (REGION == 'FC') {

    reducedDim_ID <- 'IterativeLSI_reclust'

  } else {

    reducedDim_ID <- 'IterativeLSI'
    clust_ID <- 'Clusters'
    UMAP_ID<- 'UMAP'

  }

  ##  Load ArchR project  ---------------------------------------------------------------
  cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
  archR <- loadArchRProject(paste0('../results/ARCHR/', REGION))
  
  cat(paste0('\nAdding peak matrix ... \n'))
  archR <- addPeakMatrix(archR) 
  
  ##  IDing Marker peaks that are unique to individual cluster - chptr 11  --------------
  # Cluster labels
  table(archR$Clusters_broad)
  
  cat(paste0('\nAdding marker peaks for each cell type ... \n'))
  # Add marker peaks - returns summarisedExperiment with 6 assays
  markersPeaks <- getMarkerFeatures(
    ArchRProj = archR, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters_broad",
    bias = c("TSSEnrichment", "log10(nFrags)"), # Correction fof diffs in data quality
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
  markerList
  
  ## Peak coaccessibility  -  Chptr 15.2  ----------------------------------------------
  cat(paste0("Starting coaccessibility analysis for ", REGION, " ... \n"))
  cat("Adding coaccessibile links ... \n")
  archR <- addCoAccessibility(
    ArchRProj = archR,
    reducedDims = reducedDim_ID
  )

  cat("Retrieving coaccesibility df ... \n")
  cA <- getCoAccessibility(
    ArchRProj = archR,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
  )

  cA_df <- cA
  cA_df
  cA_metadata <- metadata(cA)[[1]]
  cA_metadata

  ## Save RDS objects
  cat(paste0('Writing coaccesibility df and metadata'))
  saveRDS(cA_df, paste0(RDS_DIR, '_', REGION, 'cA_peaks_df.rds'))
  saveRDS(cA_metadata, paste0(RDS_DIR, '_', REGION, 'cA_peaks_metadata.rds'))


  ## Motif enrichment  -  Chptr 12  -----------------------------------------------------
  # Add motif anns to archR project - creates binary matrix for presence/absence of motif in peak
  cat(paste0('\nAdding motif annotations  ... \n'))
  archR <- addMotifAnnotations(ArchRProj = archR, motifSet = "cisbp", name = "Motif", force = TRUE)
  
  ##   --------  Motif enrichment in marker peaks  ----------------- 
  cat(paste0('\nAssessing motif enrichment ... \n'))
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = archR,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  enrichMotifs
  
  # Export enrichMotifs and motif p-val table used for plotting
  saveRDS(enrichMotifs, paste0(RDS_DIR, REGION, '_motifs.rds'))
  saveRDS(assays(enrichMotifs)[["mlog10Padj"]], paste0(RDS_DIR, REGION, '_motifs_mlogPs.rds'))

  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE)
  assign(paste0("motif_heatmap_", REGION), heatmapEM)
  
  ## Create motif heatmap RDS files  ----------------------------------------------------
  cat(paste0('\nCreating rds files  ... \n'))
  # Return heatmap matrix
  heatmap_matrix <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE,
                                      returnMatrix = TRUE)
  
  # Remove all extra info from motif names
  new_cols <- as.data.frame(str_split(colnames(heatmap_matrix), 
                                      fixed("_"), simplify = TRUE)) %>%
    pull(V1)
  colnames(heatmap_matrix) <- new_cols
    
  # Note that Complex heatmap also has a function called heatmap for object conversion!!
  motif_heatmap_plot <- as.ggplot(pheatmap::pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                           cellwidth = 12, cellheight = 15, fontsize_row = 12,
                           fontsize_col = 12)) 
  
  # Save motif heatmap RDS files
  saveRDS(motif_heatmap_plot, paste0(RDS_DIR, REGION, '_motif_enrichment_heatmap_plot.rds')) 
  saveRDS(heatmap_matrix, paste0(RDS_DIR, REGION, '_motif_enrichment_heatmap_matrix.rds'))    
  
  ## HARs enrichment  -  Chptr 12.4  -----------------------------------------------------
  # Add human accelerated regions anns to archR project 
  archR <- addPeakAnnotations(ArchRProj = archR, regions = HAR_regions, name = "Hars")
  
  ##   --------  HARs enrichment in marker peaks  ----------------- 
  enrichHarRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = archR,
    peakAnnotation = "Hars",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
    
  # Export enrichMotifs and motif p-val table used for plotting
  saveRDS(enrichHarRegions, paste0(RDS_DIR, REGION, '_HARs.rds'))
  saveRDS(assays(enrichHarRegions)[["mlog10Padj"]], paste0(RDS_DIR, REGION, '_HARs_mlogPs.rds'))
  
  heatmapEM <- plotEnrichHeatmap(enrichHarRegions, transpose = TRUE)
  assign(paste0("HARs_heatmap_", REGION), heatmapEM)
  
  ## Create motif heatmap RDS files  ----------------------------------------------------
  cat(paste0('\nCreating rds files  ... \n'))
  # Return heatmap matrix
  heatmap_HARs_matrix <- plotEnrichHeatmap(enrichHarRegions, transpose = TRUE,
                                           returnMatrix = TRUE)
  
  # Remove all extra info from motif names
  # new_cols <- as.data.frame(str_split(colnames(heatmap_HARs_matrix), 
  #                                     fixed("_"), simplify = TRUE)) %>%
  #   pull(V1)
  # colnames(heatmap_HARs_matrix) <- new_cols
  
  # Note that Complex heatmap also has a function called heatmap for object conversion!!
  heatmap_HARs_plot <- as.ggplot(pheatmap::pheatmap(heatmap_HARs_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                                                    cellwidth = 12, cellheight = 15, fontsize_row = 12,
                                                    fontsize_col = 12)) 
  
  # Save motif heatmap RDS files
  saveRDS(heatmap_HARs_plot, paste0(RDS_DIR, REGION, '_HARs_enrichment_heatmap_plot.rds')) 
  saveRDS(heatmap_HARs_matrix, paste0(RDS_DIR, REGION, '_HARs_enrichment_heatmap_matrix.rds'))   


  ## Motif Footprinting  -  Chptr 14  ---------------------------------------------------
  #  cat(paste0('\nRunning footprinting ... \n'))
  # Loop to pull out cisbp codes for motifs of interest
  REGION_MOTIFS <- get(paste0(REGION, '_motifs'))
  cat(paste0('\nRegion motifs are:\n'))
  cat(REGION_MOTIFS)

  cat(paste0('\nLoading motif file ... \n'))
  motifs <- enrichMotifs

  cat(paste0('\nRetrieving cisbp codes for motifs of interest in ', REGION, ' ... \n'))
  MOTIFS_recode_all <- vector()

  for (MOTIF in REGION_MOTIFS) {
  
    MOTIFS_recode <- grep(MOTIF, rownames(motifs), value = TRUE)
    print(MOTIFS_recode)
    MOTIFS_recode_all <- c(MOTIFS_recode_all, MOTIFS_recode)
    assign(paste0(REGION, '_', 'motifs_recode_all'), MOTIFS_recode_all) # For .Rmd file   
  
  }
  
  cat(paste0('\nAll ', length(REGION_MOTIFS), ' motifs have a cisbp code? ', 
           length(REGION_MOTIFS) == length(MOTIFS_recode_all),'\n'))


  cat(paste0('\nGenerating motif footprints for ', REGION, ' ... \n'))
  motifPositions <- getPositions(archR)
  motifPositions
  
  markerMotifs <- unlist(lapply(MOTIFS_recode_all, function(x) grep(x, names(motifPositions), value = TRUE)))
  markerMotifs

  # Compute footprints
  cat(paste0('Computing footprints ... \n'))
  seFoot <- getFootprints(
      ArchRProj = archR,
      positions = motifPositions[markerMotifs],
      groupBy = "Clusters_broad"
  )

  # Plot footprint - plot = FALSE required to get grob object
  cat(paste0('Plotting ... \n'))
  footprint_grob <-  plotFootprints(
      seFoot = seFoot,
      ArchRProj = archR,
      normMethod = "Subtract",
      plotName = "Footprints_subtract_bias",
      addDOC = FALSE,
      smoothWindow = 5,
      plot = FALSE

  )

  cat(paste0('\nCreating plots and rds files for footprints ... \n'))
  for (MOTIF in MOTIFS_recode_all) {

    motifs <- MOTIF
    motifs
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
    markerMotifs

    ## Convert to grob ggplot object
    cat(paste0('Creating plot footprint plot for ', MOTIF, ' ... \n'))
    footprint_plot <- ggplotify::as.ggplot(footprint_grob[[MOTIF]])

    assign(paste0(REGION, '_', MOTIF, '_footprint_plot'), footprint_plot)
    saveRDS(footprint_plot, paste0(RDS_DIR, REGION, '_', MOTIF, '_footprint_plot.rds'))
    saveRDS(footprint_grob, paste0(RDS_DIR, REGION, '_', MOTIF, '_footprint_grob.rds'))

    }

}

#  ## Peak to gene links  -  Chptr 15.3  -------------------------------------------------
#  cat(paste0("Starting peak to gene linkage analysis for ", REGION, " ... \n"))
#  cat("Adding Peak to gene links ... \n")
#  archR <- addPeak2GeneLinks(
#  ArchRProj = archR,
#  reducedDims = "Harmony"
#  )

#  cat("\nRetrieving Peak to gene links ... \n")
#  p2g <- getPeak2GeneLinks(
#  ArchRProj = archR,
#  corCutOff = 0.45,
#  resolution = 1,
#  returnLoops = FALSE
#  )

#  cat("\nCreating Peak to gene links dataframes ... \n")
#  p2g_df <- base::as.data.frame(p2g)
#  geneIDs_df <- Repitools::annoGR2DF(metadata(p2g)$geneSet)
#  peakIDs_df <- Repitools::annoGR2DF(metadata(p2g)$peakSet)
  
#  markerGenes  <- c("MYTL1", "FUT9", "DLX1", "DLX2")

#  p <- plotBrowserTrack(
#    ArchRProj = archR, 
#    groupBy = "Clusters_broad", 
#    geneSymbol = markerGenes, 
#    upstream = 50000,
#    downstream = 50000,
#    loops = getPeak2GeneLinks(archR)
#)

#  # Plot p2g links
#  tiff(paste0("results/figures/MYTL1_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
#       compression = "lzw", res = 300)
#  grid::grid.newpage()    
#  grid::grid.draw(p$MYTL1)
#  dev.off()

#  tiff(paste0("results/figures/DLX1_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
#       compression = "lzw", res = 300)
#  grid::grid.newpage()
#  grid::grid.draw(p$DLX1)
#  dev.off()

#  tiff(paste0("results/figures/DLX2_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
#       compression = "lzw", res = 300)
#  grid::grid.newpage()
#  grid::grid.draw(p$DLX2)
#  dev.off()
  

 # cat("\nObtaining peak start/stop coordinates and gene IDs for links ... \n")

  # Need to remove dfs between regions or region specific entries will be appended
#  if (exists("peak2gene_df")) { rm(peak2gene_df) }
#  if (exists("peak2gene_final_df")) { rm(peak2gene_final_df) }
  
#  for (LOOP in 1:nrow(p2g_df)) {
    
#    gene_index <- p2g_df[LOOP, 2]
#    peak_index <- p2g_df[LOOP, 1]
    
#    gene_TSS <- geneIDs_df[gene_index, ]$start
#    gene_ID <- geneIDs_df[gene_index, ]$name
#    peak_start <- peakIDs_df[peak_index, ]$start
#    peak_end <- peakIDs_df[peak_index, ]$end
#    chr <- as.vector(peakIDs_df[peak_index, ]$chr)
    
#    if (exists("peak2gene_df")) {
      
#      peak2gene_row <- cbind(chr, peak_start, peak_end, gene_TSS, gene_ID)
#      peak2gene_df <- rbind(peak2gene_df, peak2gene_row)
      
#    } else {
      
#      peak2gene_df <- cbind(chr, peak_start, peak_end, gene_TSS, gene_ID)
      
#    }
    
    
#  }
  
#  # Join peak start/stop coordinates and gene IDs to original table
#  peak2gene_final_df <- cbind(peak2gene_df, p2g_df)
  
  
#  cat("\nWriting table to file ... \n")
#  write_tsv(peak2gene_final_df, paste0('results/peak2gene_table_', REGION, '.tsv'))
  
#}  


## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
