#--------------------------------------------------------------------------------------
#
#     Annotate ArchR corrrolated coaccesible peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# AnnotationDbi - https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
# Granges - https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
# ChIPseeker - https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Info  ------------------------------------------------------------------------------

# 1. Read in snATACseq cA peaks and metadata
# 2. Combine all peak information for cA correlated peaks in single df
# 3. Find cA peak pairs where one peak contains PGC3 SCZ finemapped SNP
# 4. Annotate cA peak pairs containing SNPs
# 5. Extract correlated peaks where peak NOT containing SNP are annotated as promoter
# 6. Pull out correlation information for specific peak pairs where one contains a SNP

# cA_df - queryHits and subjectHits cols denote index of the two correlated peaks
# metadata_df - indexes of queryHits and subjectHits mentioned above apply to this


##  Load Packages  --------------------------------------------------------------------
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
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/fine_mapped_SNPs/"
REGIONS <- c("FC", "GE")
dir.create(CA_DIR)


#CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")

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
  CA_PEAKS_NO_CORR <- cbind(query_peaks_info_with_chr, subject_peaks_info_with_chr)
  
  # Check if chr cols identical
  identical(CA_PEAKS_NO_CORR [,1], CA_PEAKS_NO_CORR [,4])
  colnames(CA_PEAKS_NO_CORR ) <- c("chr", "query_start", "query_end", "chr_ext", 
                          "subject_start", "subject_end")
  CA_PEAKS_NO_CORR  <- CA_PEAKS_NO_CORR  %>%
    select(-chr_ext)
  
  CA_PEAKS_WITH_CORR <- cbind(query_peaks_info_no_chr, subject_peaks_info_no_chr, CA_PEAKS) 
  
  cat('Total cA peak pairs in ', REGION, ' after conversion: ', nrow(CA_PEAKS), '\n')
  
  assign(paste0(REGION, '_cA_peaks_combined'), CA_PEAKS_NO_CORR)
  assign(paste0(REGION, '_cA_peaks_combined_with_corr_results'), CA_PEAKS_WITH_CORR)
  assign(paste0(REGION, '_cA_query_peaks_combined'), query_peaks_info_with_chr)
  assign(paste0(REGION, '_cA_subject_peaks_combined'), subject_peaks_info_with_chr)
  
  #write_tsv(CA_PEAKS, paste0(CA_DIR, REGION, '_cA_peaks_combined.tsv'))
  
}

##  Find cA peak pairs where one peak contains PGC3 SCZ finemapped SNP  ---------------
cat('\n\nFind cA peak pairs where one peak contains PGC3 SCZ finemapped SNP ... \n')
for (REGION in REGIONS) {
  
  # Note there may be a fair bit of redundancy here
  # sum(sort(subject_peaks_idx) == sort(query_peaks_idx)) = 153886
  # length(query_peaks_idx) = 153894
  
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
  
  # Combine overlapping peak and SNP info in single df - only base position at this stage 
  METADATA_OVERLAPS <- as.data.frame(ALL_CA_PEAKS[queryHits(OVERLAPS)]) %>%
    select(seqnames, start, end)
  SNP_OVERLAPS <- as.data.frame(SNPS_granges[subjectHits(OVERLAPS)]) %>%
    select(seqnames, start) %>%
    rename(hg38_base_position = start)
  OVERLAPS_DF <- cbind(METADATA_OVERLAPS, SNP_OVERLAPS) 
  colnames(OVERLAPS_DF) <- c("chr", "start", "end", "seqnames", "hg38_base_position")
  
  # Do chr columns match between SNPs and peaks?
  identical(as.vector(OVERLAPS_DF$chr), as.vector(OVERLAPS_DF$seqnames))

  # Rm duplicate chr column
  OVERLAPS_DF <- OVERLAPS_DF %>%
    select(-seqnames)
  
  # Double check overlaps
  for (SNP in 1:nrow(OVERLAPS_DF)) {
    
    test <- ifelse(OVERLAPS_DF$hg38_base_position[SNP] >= OVERLAPS_DF$start[SNP] & OVERLAPS_DF$hg38_base_position[SNP] <= OVERLAPS_DF$end[SNP], TRUE, FALSE) 
    print(test)
    
  }
  
  # Join with rsIDs - increase in entries here as multiple SNPs in single peaks
  OVERLAPS_DF_JOIN <- OVERLAPS_DF %>% left_join(SNPS) %>%
    select(chr, start, end, hg38_base_position, rsid)
  
  # Unique SNPs in final df
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata:',
      length(unique(OVERLAPS_DF_JOIN$rsid)), '\n')
 
  # Pull out peaks from query and subject in combined metadata peaks
  QUERY_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'query_start', 'end' = 'query_end'))
  
  SUBJECT_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'subject_start', 'end' = 'subject_end'))
  
  # Note that the QUERY_JOIN and SUBJECT_JOIN objects are the same - the metadata object must contain
  # the peak information both ways
  identical(QUERY_JOIN, SUBJECT_JOIN %>% rename(subject_start = query_start, subject_end = query_end))

  # Unique SNPs after joining with cA peaks
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata query:',
      length(unique(QUERY_JOIN$rsid)), '\n')
  
  assign(paste0(REGION, '_cA_peak_pairs_with_SNP'), QUERY_JOIN)

}

##  Annotate cA peak pairs containing SNPs  -------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat('\n\nRunning regulatory peak annotation for cA peaks pairs containing SNPs ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  # Load peaks - need to remove chr to match edb encoding for annotatePeaks
  cA_PEAKS <- get(paste0(REGION, '_cA_peak_pairs_with_SNP')) #%>%
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
  
  cA_PEAKS_with_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_with_SNP_ANN@anno) %>%
    mutate(rsid = cA_PEAKS$rsid,
           hg38_base_position = cA_PEAKS$hg38_base_position) %>%
    relocate(seqnames, start, end, hg38_base_position, rsid) %>%
    select(-strand, -width)
  
    
  cA_PEAKS_no_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_no_SNP_ANN@anno)
  
  # Add SNPs info back into df
  
  assign(paste0(REGION, '_cA_peak_pairs_with_SNP_ann'), cA_PEAKS_with_SNP_ANN_DF)
  assign(paste0(REGION, '_cA_peak_pairs_no_SNP_ann'), cA_PEAKS_no_SNP_ANN_DF)
  
  #write_tsv(cA_PEAKS_with_SNP_ANN_DF, paste0(CA_DIR, REGION, '_cA_regulatory_anns_with_SNP.tsv'))
  #write_tsv(cA_PEAKS_no_SNP_ANN_DF, paste0(CA_DIR, REGION, '_cA_regulatory_anns_no_SNP.tsv'))
  
}


# Extract corr peaks where peak NOT containing SNP are annotated as promoter ----------
for (REGION in REGIONS) {
  
  cA_with_SNPs <- get(paste0(REGION, '_cA_peak_pairs_with_SNP_ann'))
  cA_no_SNPs <- get(paste0(REGION, '_cA_peak_pairs_no_SNP_ann'))
  
  # Extract indices for peaks with no SNP annotated as promoter
  cA_no_SNPs_promoter_idx <- which(cA_no_SNPs$annotation == 'Promoter')
  
  # Pull indices out of correlated peak tables
  cA_with_SNPs_promoter_in_correlated_peak <-  cA_with_SNPs[cA_no_SNPs_promoter_idx,]
  cA_no_SNPs_promoter_in_this_peak <- cA_no_SNPs[cA_no_SNPs_promoter_idx,]
  
  # write_tsv(cA_with_SNPs_promoter_in_correlated_peak, 
  #           paste0(CA_DIR, REGION, '_cA_regulatory_anns_with_SNP_promoters_only_in_corrolated_peak.tsv'))
  # write_tsv(cA_no_SNPs_promoter_in_this_peak, 
  #           paste0(CA_DIR, REGION, '_cA_regulatory_anns_no_SNP_promoters_only.tsv'))
  # 
  # Check the indices are the same
  identical(rownames(cA_no_SNPs_promoter_in_this_peak), rownames(cA_with_SNPs_promoter_in_correlated_peak))
  
}

##  Extract correlation information for specific peak pairs where one contains a SNP --



FC_cA_PEAKS_COR <- FC_cA_peaks_combined_with_corr_results %>%
  rename(chr = seqnames)
GE_cA_PEAKS_COR <- GE_cA_peaks_combined_with_corr_results %>%
  rename(chr = seqnames)

fc_SNPs <- data.frame(
  chr = c('chr2', 'chr2', 'chr2'),
  query_start = c(172090681, 172092453, 172092453),
  query_end = c(172091181, 172092953,172092953),
  hg38_base_position = c(172090705, 172092806,172092915),
  rsid = c('rs62183854', 'rs13388257', 'rs13390848'))

ge_SNPs <- data.frame(
  chr = c('chr2'),
  start = c(172091234),
  end = c(172091734),
  hg38_base_position = c(172091721),
  rsid = c('rs62183855'))

FC_cA_PEAKS_COR_OVERLAPS <- fc_SNPs %>% inner_join(FC_cA_PEAKS_COR)
GE_cA_PEAKS_COR_OVERLAPS <- ge_SNPs %>% inner_join(GE_cA_PEAKS_COR)

write_tsv(FC_cA_PEAKS_COR_OVERLAPS, paste0(CA_DIR, 'FC_cA_DLX1_SNP_overlaps.tsv'))
write_tsv(GE_cA_PEAKS_COR_OVERLAPS, paste0(CA_DIR, 'GE_cA_DLX1_SNP_overlaps.tsv'))

GEcA_peak_pairs_with_SNP_ann %>%
  filter(grepl('rs62183855', rsid))

unique(GE_cA_PEAKS_COR_OVERLAPS$rsid)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
