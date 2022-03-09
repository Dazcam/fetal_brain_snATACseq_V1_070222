#--------------------------------------------------------------------------------------
#
#     Map PGC3 SCZ index SNPs and proxies with R2 > 0.8 to snATACseq peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# LDlinkR - requires a personal access token
# Github - https://github.com/CBIIT/LDlinkR
# Vingette: https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html

## Info  ------------------------------------------------------------------------------

# 1. Download PGC3 SCZ GWAS index SNPs
# 2. Extract proxy SNPs of PGC3 index SNPs with r2 > 0.8 using LDlinkR - hg19
#    - 9 SNPs failed 
#    - Not encoded as an rsID in xslx file: 8:4180090_T_A 
#    - Not in 1000G reference panel: rs2494638, rs12514660, rs2914025, rs55858161
#    - Not a not a biallelic variant: rs62152282, rs35026989, rs11241041, rs650520
# 3. Map SNPs to hg38 using BioMart
# 4. Check for overlap of PGC3 index and proxy SNPs in snATACseq peaks (500 bp ext)

##  Load Packages  --------------------------------------------------------------------
library(LDlinkR)  
library(readxl)
library(tidyverse)
library(Repitools)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/PGC3_snps/"
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/archR_data_processing/peaks/"
REGIONS <- c("FC", "GE")
CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Create directories  ----------------------------------------------------------------
dir.create(SNP_DIR)
setwd(SNP_DIR) # Required as LDproxy_batch outputs to cwd

##  Read in PGC3 index SNPs one is not an rsID just removed it  ------------------------
snps <- read_excel(paste0(IN_DIR, 'PGC3_Sup_table_3_combined_discovery_replication_loci_jul21.xls')) %>%
  dplyr::select(`top-index`) %>%
  filter(!grepl('8:4180090_T_A', `top-index`)) %>%
  base::as.data.frame(snps)

# Extract proxy SNPs in LD - can't specify threshold of 0.8 for this so need to go
# With default of 0.01. LD taken from 1000Gs phase 3 - CEU population
snps_in_LD <- LDlinkR::LDproxy_batch(snp = snps, pop = "CEU", r2d = "r2", token = TOKEN)

## Remove SNPs that failed proxy batch call from snp list -----------------------------
error_snps <- c('rs2494638', 'rs12514660', 'rs2914025', 'rs55858161',
                'rs62152282', 'rs35026989', 'rs11241041', 'rs650520')
snps_no_errors <- setdiff(snps, error_snps)

##  Get proxy SNPs with R2 > 0.8 for all index SNPs. ----------------------------------
for (SNP in snps$`top-index`) {
  
  skip_to_next <- FALSE
  cat(paste0('\nLoading proxy SNPs for ', SNP, '... \n'))
  
  tryCatch(
    
    
    proxies <<- read_delim(paste0(SNP_DIR, SNP, ".txt"), delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE),
    
    error = function(e) {
      
      cat(paste0("No file for: ", SNP, "\n"))
      skip_to_next <<- TRUE })
  
  if (skip_to_next) { next }     
  
  proxies_r2_0.8 <- proxies %>%
    filter(R2 > 0.8)
  
  # Generate stats
  index_in_proxies <- SNP %in% proxies$Coord
  index_in_proxies_r2_0.8 <- SNP %in% proxies_r2_0.8$Coord
  snps_in_proxies <- nrow(proxies)
  snps_in_proxies_r2_0.8 <- nrow(proxies_r2_0.8)
  cat(paste0('\n', nrow(proxies), ' snps in ', SNP, ' proxies table ... \n'))
  cat(paste0('\nChecking if ', SNP, ' index SNP is in proxies table: ', index_in_proxies, '\n'))
  cat(paste0('\nChecking if ', SNP, ' index SNP is in proxies r2 > 0.8 table: ', index_in_proxies_r2_0.8, '\n'))
  
  if (exists("snp_proxies_df")) {
    
    snp_proxies_df <- rbind(snp_proxies_df, proxies_r2_0.8) 
    snp_stats <- cbind(SNP, index_in_proxies, snps_in_proxies, snps_in_proxies_r2_0.8)
    snp_proxies_summary_df <- rbind(snp_proxies_summary_df, snp_stats)
    
  } else {
    
    snp_proxies_df <- proxies_r2_0.8
    snp_proxies_summary_df <- cbind(SNP, index_in_proxies, snps_in_proxies, snps_in_proxies_r2_0.8)
    
  }
  
  
}

# Retain only unique SNPs in list 
cat('\nRetain only unique SNPs  ... \n')
snp_proxies_unique_df <- snp_proxies_df %>% distinct(Coord, .keep_all = TRUE) %>%
  filter(grepl('rs', Coord)) # Retain SNPs with rs start only
scz_snps <- gtools::mixedsort(as.vector(snp_proxies_unique_df$Coord))


cat(paste0(length(scz_snps), ' SNPs retained.\n'))
write_tsv(as.data.frame(scz_snps), 'scz_snps.tsv')
scz_snps <- read_tsv(paste0(SNP_DIR, 'scz_snps.tsv')) %>%
  pull()

# Tested rsnps package - too slow - only good for a hnadful of SNPs
# library(rsnps) # ~83hrs for all ~200K SNPs
# SNPs <- rsnps::ncbi_snp_query(scz_snps)
# rsnps::ncbi_snp_query('rs473')
# 
# start_time <- Sys.time()
# snps <- c(scz_snps[500:510])
# rsnps::ncbi_snp_query(snps)
# end_time <- Sys.time()
# end_time - start_time

## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
## ~15-40 mins per 50K SNPs - note sometimes crashes when BiomaRt in heavy use
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs  ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# Need to batch the query: https://support.bioconductor.org/p/23684/
SNPs_a <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[1:50000], 
                mart = ensembl, 
                uniqueRows=TRUE)
Sys.sleep(1)
SNPs_b <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[50001:100000], 
                mart = ensembl, 
                uniqueRows=TRUE)
Sys.sleep(1)
SNPs_c <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[100001:150000], 
                mart = ensembl, 
                uniqueRows=TRUE)
Sys.sleep(1)
SNPs_d <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[150001:200000], 
                mart = ensembl, 
                uniqueRows=TRUE)
SNPs_e <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[200001:length(scz_snps)], 
                mart = ensembl, 
                uniqueRows=TRUE)
SNPs <- rbind(SNPs_a, SNPs_b, SNPs_c, SNPs_d, SNPs_e)

cat(paste0(nrow(SNPs), ' SNPs retained.\n'))

# Some SNPs include duplicates due to chr patches - I removed these
cat('\nRemoving SNPs on CHR patches ... \n')
snps_no_patches <- SNPs %>%
  filter(!grepl('_', chr_name)) %>%
  dplyr::select(-chrom_end) %>%
  dplyr::rename('snpID' = refsnp_id, 'hg38_base_position' = chrom_start)
cat(paste0(nrow(snps_no_patches), ' SNPs retained. \n'))

# 
write_tsv(as.data.frame(snps_no_patches), '_PGC3_SCZ_r2_0.8_SNPs_hg38.tsv')

## Check for overlap of SNPs in snATACseq peaks of all cell types  --------------------
for (CELL_TYPE in CELL_TYPES) {
  
  cat(paste0('\nLoading peaks for ', CELL_TYPE, ' ... \n'))
  peaks_df <- read_tsv(paste0(PEAK_DIR, CELL_TYPE,'.hg38.bed'), col_names = FALSE)
  colnames(peaks_df) <- c('chr', 'start', 'end', 'name', 'score', 'strand')

  cell_overlaps <- data.frame()
  
  cat(paste0('\nChecking for SNP overlaps in ', CELL_TYPE, ' ... \n'))
  for (i in 1:nrow(snps_no_patches)) {
    
    BASE_POSITION <- snps_no_patches$hg38_base_position[i]
    CHR <- snps_no_patches$chr_name[i]
    cat(paste0('\nSNP: ', 
               snps_no_patches$snpID[i], ', position: ',
               BASE_POSITION, '... \n'))
    
    overlaps <- filter(peaks_df, start <= BASE_POSITION, end >= BASE_POSITION, chr == paste0('chr', CHR))
    
    if (nrow(overlaps) > 0) { 
      
      overlaps <- cbind(overlaps, snps_no_patches[i,])
      print(overlaps) 
      
    }
    
    cell_overlaps <- rbind(cell_overlaps, overlaps)
    
  } 
  
  cat(paste0('\nAll ', nrow(snps_no_patches), ' SNPs checked in ', CELL_TYPE, ' ... \n'))
  
  cat(paste0('\nWriting overlapping SNPs to file ... \n'))
  write_tsv(cell_overlaps, paste0(PEAK_DIR, CELL_TYPE,'_PGC3_SCZ_r2_0.8_SNP_overlaps.tsv'))
  assign(paste0(CELL_TYPE, '_SNP_overlaps'), cell_overlaps)
  
}

cat('Done.\n')

## Need to cross reference this with scripts/R/snATACseq_map_PGC3_snps_to_peaks.R on cluster

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------