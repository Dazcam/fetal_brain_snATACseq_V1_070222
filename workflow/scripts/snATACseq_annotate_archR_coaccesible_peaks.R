#--------------------------------------------------------------------------------------
#
#     Coaccessibility - 
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# AnnotationDbi - https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
# Granges - https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
# ChIPseeker - https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Info  ------------------------------------------------------------------------------

# 1. Coaccessibility - prioritisd SNPs - 250kp around each SNP any gene promoters 
# 2.  that are coaccessible we want to know what they are either side (i.e. within 500kb)

# cA_df - queryHits and subjectHits cols denote index of the two corrolated peaks
# metadata_df - indexes of queryHits and subjectHits mentioned above apply to this


##  Load Packages  --------------------------------------------------------------------
biocLite('FunciSNP')
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86) #hg38
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(Repitools) # Granges to df
library(diffloop) # Add/rm chr to granges seqlevels



##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/archR_data_processing/rds_files/"
CA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/coaccessibile_peaks/"
REGIONS <- c("FC", "GE")
dir.create(CA_DIR)
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/fine_mapped_SNPs/"

CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Load data  ----------------------------------------------------------------
cat('\nLoading peak and SNP data ... \n')
for (REGION in REGIONS) {
  
  PEAKS <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_df.rds'))
  METADATA <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_metadata.rds'))
  
  assign(paste0(REGION, '_cA_peaks'), PEAKS)
  assign(paste0(REGION, '_cA_metadata'), METADATA)
    
}

SNPS <- read_tsv(paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_SNPs_only.tsv'))
SNPS_WITH_PEAKS <- read_tsv(paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_with_peaks.tsv'))

##  Combine info on cA peaks and metadata dfs  ----------------------------------------
# Get peak indices for query and subject hits in index files - these may be the same - check
cat('\n\nCombine peak info stored in cA peaks and metadata dfs ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  CA_PEAKS <- get(paste0(REGION, '_cA_peaks'))
  CA_METADATA <- get(paste0(REGION, '_cA_metadata'))

  query_peaks_idx <- CA_PEAKS$queryHits
  subject_peaks_idx <- CA_PEAKS$subjectHits
  
  # Pull peak information for indices from metadata
  query_peaks_info <- CA_METADATA[query_peaks_idx] 
  subject_peaks_info <- CA_METADATA[subject_peaks_idx] 
  
  # Remove rownames
  names(query_peaks_info) <- NULL
  names(subject_peaks_info) <- NULL
  
  # Convert peak info to df
  query_peaks_info_no_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(start, end)
  subject_peaks_info_no_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(start, end)
  colnames(query_peaks_info_no_chr) <- c("query_start", "query_end")
  colnames(subject_peaks_info_no_chr) <- c("subject_start", "subject_end")
  
  # Save peak info with chr
  query_peaks_info_with_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(chr, start, end)
  subject_peaks_info_with_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(chr, start, end)
  
  # Convert peaks S4 df to S3
  CA_PEAKS <- as.data.frame(CA_PEAKS)
  CA_PEAKS <- cbind(query_peaks_info_with_chr, subject_peaks_info_with_chr)
  
  # Check if chr cols identical
  identical(CA_PEAKS[,1], CA_PEAKS[,4])
  colnames(CA_PEAKS) <- c("chr", "query_start", "query_end", "chr_ext", 
                          "subject_start", "subject_end")
  CA_PEAKS <- CA_PEAKS %>%
    select(-chr_ext)
  
  
  CA_PEAKS_WITH_CORR <- cbind(query_peaks_info_no_chr, subject_peaks_info_no_chr, CA_PEAKS) 
  
  cat('Total cA peak pairs in ', REGION, ' after conversion: ', nrow(CA_PEAKS), '\n')
  
  assign(paste0(REGION, '_cA_peaks_combined'), CA_PEAKS)
  assign(paste0(REGION, '_cA_peaks_combined_with_corr_results'), CA_PEAKS_WITH_CORR)
  assign(paste0(REGION, '_cA_query_peaks_combined'), query_peaks_info_with_chr)
  assign(paste0(REGION, '_cA_subject_peaks_combined'), subject_peaks_info_with_chr)
  
  #write_tsv(CA_PEAKS, paste0(CA_DIR, REGION, '_cA_peaks_combined.tsv'))
  
}

##  Find cA peak pairs where one peak contains PGC3 SCZ finemapped SNP  ---------------
cat('\n\nFind cA peak pairs where one peak contains PGC3 SCZ finemapped SNP ... \n')
for (REGION in REGIONS) {
  
  # Pull out ranges from the ArchR cA metadata dataframe
  ALL_CA_PEAKS <- get(paste0(REGION, '_cA_metadata'))
  PEAKS_COMBINED <- get(paste0(REGION, '_cA_peaks_combined'))
  names(ALL_CA_PEAKS) <- NULL
  
  # Create GRanges object for all PGC3 SCZ SNPs that overlap snATACseq peak
  # use only hg38 base position - checked if all BPs are unique in SNPs
  cat('\n\nCreate GRanges object for finemapped SNPs ... \n')
  SNPS_add_range <- SNPS %>%
    mutate(end = hg38_base_position) %>%
    rename(start = hg38_base_position) %>%
    relocate(rsid, chr, start, end)
  SNPS_granges <- makeGRangesFromDataFrame(SNPS_add_range, keep.extra.columns = FALSE,
                                           ignore.strand = TRUE)
  
  # Find overlaps using IRanges
  OVERLAPS <- findOverlaps(ALL_CA_PEAKS, SNPS_granges)
  
  # Combine peak and SNP info in single df
  METADATA_OVERLAPS <- as.data.frame(ALL_CA_PEAKS[queryHits(OVERLAPS)]) %>%
    select(seqnames, start, end)
  SNP_OVERLAPS <- as.data.frame(SNPS_granges[subjectHits(OVERLAPS)]) %>%
    select(seqnames, start) %>%
    rename(hg38_base_position = start)
  OVERLAPS_DF <- cbind(METADATA_OVERLAPS, SNP_OVERLAPS) 
  colnames(OVERLAPS_DF) <- c("chr", "start", "end", "seqnames", "hg38_base_position")
  
  # Do chr columns match between SNPs and peaks?
  identical(as.vector(OVERLAPS_DF$chr), as.vector(OVERLAPS_DF$seqnames))
  
  # Add maual check for overlap here
  
  # Rm duplicate column
  OVERLAPS_DF <- OVERLAPS_DF %>%
    select(-seqnames)
  
  
  # Join with rsIDs - increase in entries here as mulitiple SNPs in single peaks
  OVERLAPS_DF_JOIN <- OVERLAPS_DF %>% left_join(SNPS) %>%
    select(chr, start, end, hg38_base_position, rsid)
  
  # Unique SNPs in final df
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata:',
      length(unique(OVERLAPS_DF_JOIN$rsid)), '\n')
 
  QUERY_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'query_start', 'end' = 'query_end'))
  
  SUBJECT_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'subject_start', 'end' = 'subject_end'))
  
  # Unique SNPs after joining with cA peaks
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata query:',
      length(unique(QUERY_JOIN$rsid)), '\n')
  
  assign(paste0(REGION, 'cA_peak_pairs_with_SNP'), QUERY_JOIN)

}

##  Annotate cA peak pairs conatining SNPs  -------------------------------------------
edb <- EnsDb.Hsapiens.v86
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat('\n\nRunning regulatory peak annotation for cA peaks pairs containing SNPs ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  # Load peaks - need to remove chr to match edb encoding for annotatePeaks
  cA_PEAKS <- get(paste0(REGION, 'cA_peak_pairs_with_SNP')) #%>%
   # mutate(chr = gsub("chr", "", chr))
  
  # Split df to run each set of cA peaks separately
  cat('\n\nSplit df to run each set of cA peaks separately ... \n')
  cA_PEAKS_with_SNP <- cA_PEAKS %>% 
    select(-subject_start, -subject_end)
  cA_PEAKS_no_SNP <- cA_PEAKS %>% 
    select(-start, -end) %>%
    rename(start = subject_start,
           end = subject_end)
    
  # Create GRanges object for each peak set
  cat('\n\nCreate GRanges object for each peak set ... \n')
  cA_PEAKS_with_SNP_GR <- makeGRangesFromDataFrame(cA_PEAKS_with_SNP, keep.extra.columns = FALSE,
                                                   ignore.strand = TRUE)
  cA_PEAKS_no_SNP_GR <- makeGRangesFromDataFrame(cA_PEAKS_no_SNP, keep.extra.columns = FALSE,
                                                  ignore.strand = TRUE)
  
  # Annotate peaks
  cat('\n\nAnnotating peaks with SNP ... \n\n\n')
  cA_PEAKS_with_SNP_ANN <- annotatePeak(cA_PEAKS_with_SNP_GR, tssRegion=c(-1000, 100),
                                  TxDb = txdb, annoDb = "org.Hs.eg.db", level = 'gene')
  print(cA_PEAKS_with_SNP_ANN)
  
  cat('\n\nAnnotating peaks without SNP ... \n\n')
  cA_PEAKS_no_SNP_ANN <- annotatePeak(cA_PEAKS_no_SNP_GR, tssRegion=c(-1000, 100),
                                  TxDb = txdb, annoDb = "org.Hs.eg.db", level = 'gene')
  print(cA_PEAKS_no_SNP_ANN)
  
  print(head(cA_PEAKS_with_SNP_ANN@anno))
  print(head(cA_PEAKS_no_SNP_ANN@anno))
  
  cA_PEAKS_with_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_with_SNP_ANN@anno)
  cA_PEAKS_no_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_no_SNP_ANN@anno)
  
  # Add SNPs info back into df
  
  
  
  assign(paste0(REGION, 'cA_peak_pairs_with_SNP_ann'), cA_PEAKS_with_SNP_ANN_DF)
  assign(paste0(REGION, 'cA_peak_pairs_with_SNP_ann'), cA_PEAKS_no_SNP_ANN_DF)
  
  write_tsv(cA_PEAKS_with_SNP_ANN_DF, paste0(CA_DIR, REGION, '_cA_regulatory_anns_with_SNP.tsv'))
  write_tsv(cA_PEAKS_no_SNP_ANN_DF, paste0(CA_DIR, REGION, '_cA_regulatory_anns_no_SNP.tsv'))
  
}

##  Annotate peaks in metadata with regulatory info  ----------------------------------
# cat('\n\nRunning regulatory peak annotation of metadata ... \n')
# for (REGION in REGIONS) {
#   
#   cat('\nRunning ', REGION, ' ...\n')
#   
#   CA_METADATA <- get(paste0(REGION, '_cA_metadata'))
#   
#   # Need to remove chr from granges to match edb encoding
#   CA_METADATA <- diffloop::rmchr(CA_METADATA)
#   
#   # annotatePeak() doesn't like identical rownames
#   mcols(CA_METADATA)$cell_type <- names(CA_METADATA)
#   names(CA_METADATA) <- NULL
#   
#   # Annotate peaks
#   CA_METADATA_ANN <- annotatePeak(CA_METADATA, tssRegion=c(-1000, 100),
#                                   TxDb=edb, annoDb="org.Hs.eg.db", level = 'gene')
#   
#   CA_METADATA_ANN_DF  <- as.data.frame(CA_METADATA_ANN@anno)
#   #  write_tsv(CA_METADATA_ANN_DF, paste0(CA_DIR, REGION, '_cA_metadata_regulatory_anns.tsv'))
#   
# }

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
