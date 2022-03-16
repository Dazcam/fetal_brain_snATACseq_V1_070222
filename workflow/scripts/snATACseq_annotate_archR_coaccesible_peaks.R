#--------------------------------------------------------------------------------------
#
#     Coaccessibility - 
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# AnnotationDbi - https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
# Granges - https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html

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
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/PGC3_snps/"

CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Load data  ----------------------------------------------------------------
for (REGION in REGIONS) {
  
  PEAKS <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_df.rds'))
  METADATA <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_metadata.rds'))
  
  assign(paste0(REGION, '_cA_peaks'), PEAKS)
  assign(paste0(REGION, '_cA_metadata'), METADATA)
    
}

##  Combine info on cA peaks and metadata dfs  ----------------------------------------
# Get peak indices for query and subject hits in index files
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
  
  # Move names to column as df conversion chokes due to identical rownames
  mcols(query_peaks_info)$query_cell_type <- names(query_peaks_info)
  mcols(subject_peaks_info)$subject_cell_type <- names(subject_peaks_info)
  names(query_peaks_info) <- NULL
  names(subject_peaks_info) <- NULL
  
  # Convert peak info to df
  query_peaks_info_no_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(query_cell_type, start, end)
  subject_peaks_info_no_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(subject_cell_type, start, end)
  colnames(query_peaks_info_no_chr) <- c("query_cell_type", "query_start", "query_end")
  colnames(subject_peaks_info_no_chr) <- c("subject_cell_type", "subject_start", "subject_end")
  
  # Save peak info with chr
  query_peaks_info_with_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(chr, start, end, query_cell_type)
  subject_peaks_info_with_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(chr, start, end, subject_cell_type)
  
  
  # Convert peaks S4 df to S3
  CA_PEAKS <- as.data.frame(CA_PEAKS)
  CA_PEAKS <- cbind(query_peaks_info_no_chr, subject_peaks_info_no_chr, CA_PEAKS)
  
  peaks_diff_cell_type_rm <- sum(CA_PEAKS$query_cell_type == CA_PEAKS$subject_cell_type)
  
  cat('Total cA peak pairs in ', REGION, ' after conversion: ', nrow(CA_PEAKS), '\n')
  cat('Total cA peak pairs in ', REGION, ' after corralation between peaks from different cell types removed: ', peaks_diff_cell_type_rm, '\n')
  
  assign(paste0(REGION, '_cA_peaks_combined'), CA_PEAKS)
  assign(paste0(REGION, '_cA_query_peaks_combined'), query_peaks_info_with_chr)
  assign(paste0(REGION, '_cA_subject_peaks_combined'), subject_peaks_info_with_chr)
  
  write_tsv(CA_PEAKS, paste0(CA_DIR, REGION, '_cA_peaks_combined.tsv'))
  
}

##  Read in PGC3 index SNPs one is not an rsID just removed it  ------------------------
edb <- EnsDb.Hsapiens.v86

cat('\n\nRunning regulatory peak annotation of metadata ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  CA_METADATA <- get(paste0(REGION, '_cA_metadata'))
  
  # Need to remove chr from granges to match edb encoding
  CA_METADATA <- diffloop::rmchr(CA_METADATA)

  # annotatePeak() doesn't like identical rownames
  mcols(CA_METADATA)$cell_type <- names(CA_METADATA)
  names(CA_METADATA) <- NULL

  # Annotate peaks
  CA_METADATA_ANN <- annotatePeak(CA_METADATA, tssRegion=c(-1000, 1000),
                                            TxDb=edb, annoDb="org.Hs.eg.db")
  
  CA_METADATA_ANN_DF  <- as.data.frame(CA_METADATA_ANN@anno)
  write_tsv(CA_METADATA_ANN_DF, paste0(CA_DIR, REGION, '_cA_metadata_regulatory_anns.tsv'))

}


FC_cA_metadata_chr_rm_genes <- plotDistToTSS(FC_cA_metadata_chr_rm)
vennpie(FC_cA_metadata_chr_rm_ann)
upsetplot(FC_cA_metadata_chr_rm_ann, vennpie = TRUE)


diffloop::addchr(FC_cA_metadata)

addGeneIDs(annotatedPeak[1:6], orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"))
