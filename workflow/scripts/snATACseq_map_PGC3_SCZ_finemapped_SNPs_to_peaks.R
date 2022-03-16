#--------------------------------------------------------------------------------------
#
#     Map PGC3 SCZ finemapped SNPs to snATACseq peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# LDlinkR - requires a personal access token
# Github - https://github.com/CBIIT/LDlinkR
# Vingette: https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html

## Info  ------------------------------------------------------------------------------

# 1. Download PGC3 SCZ GWAS fine mapped SNPs
# 2. Extract proxy SNPs of PGC3 index SNPs with r2 > 0.8 using LDlinkR - hg19
#    - 9 SNPs failed 
#    - Not encoded as an rsID in xslx file: 8:4180090_T_A 
#    - Not in 1000G reference panel: rs2494638, rs12514660, rs2914025, rs55858161
#    - Not a not a biallelic variant: rs62152282, rs35026989, rs11241041, rs650520
# 3. Map SNPs to hg38 using BioMart
# 4. Check for overlap of PGC3 index and proxy SNPs in snATACseq peaks (500 bp ext)
# 5. Create binary df for whether SNP is in/not in peak for all cell types

##  Load Packages  --------------------------------------------------------------------
library(LDlinkR)  
library(readxl)
library(tidyverse)
library(Repitools)
library(biomaRt)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/PGC3_snps/fine_mapped/"
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/"
REGIONS <- c("FC", "GE")
CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Create directories  ----------------------------------------------------------------
dir.create(SNP_DIR)
dir.create(paste0(PEAK_DIR, 'fine_mapped_SNPs/'))

##  Read in PGC3 index SNPs one is not an rsID just removed it  ------------------------
snps <- read_excel(paste0(IN_DIR, 'PGC3_SCZ_Supplementary_Table_11_FINEMAP_UPDATED.xlsx'), sheet = 'ST11a 95% Credible Sets') %>%
  dplyr::select(rsid) %>%
  filter(!grepl(':|_', rsid)) %>% # 23 SNPs with 1:28690628_T_C  encoding removed for now
  base::as.data.frame(snps) %>%
  distinct(rsid, .keep_all = TRUE) %>%
  arrange(rsid) %>%
  pull()

cat(paste0(length(snps), ' SNPs retained from PGC3 table.\n'))

## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
## ~15-40 mins per 50K SNPs - note sometimes crashes when BiomaRt in heavy use
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs  ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# Need to batch the query: https://support.bioconductor.org/p/23684/
SNPs <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = snps, 
                mart = ensembl, 
                uniqueRows=TRUE)

cat(paste0(nrow(SNPs), ' SNPs retained after lift over to hg38.\n'))

# Some SNPs include duplicates due to chr patches - I removed these
cat('\nRemoving SNPs on CHR patches ... \n')
snps_no_patches <- SNPs %>%
  filter(!grepl('_', chr_name)) %>%
  dplyr::select(-chrom_end) %>%
  dplyr::rename('rsid' = refsnp_id, 'hg38_base_position' = chrom_start)
cat(paste0(nrow(snps_no_patches), ' SNPs retained. \n'))

# Add posterior probability and index SNPs columns 
SNPs_join <- snps_no_patches %>% 
  left_join(snps <- read_excel(paste0(IN_DIR, 'PGC3_SCZ_Supplementary_Table_11_FINEMAP_UPDATED.xlsx'), 
                               sheet = 'ST11a 95% Credible Sets')) %>%
  dplyr::select(rsid, chr_name, hg38_base_position, index_snp, finemap_posterior_probability)

# Write to file
write_tsv(as.data.frame(SNPs_join), 'pgc3_scz_finemapped_SNPs_hg38.tsv')
# snps_no_patches <- read_tsv(paste0(SNP_DIR, 'pgc3_scz_index_snps_and_proxies_hg38.tsv'),
#                             col_types = cols(chr_name = col_character()))
  
## Check for overlap of SNPs in snATACseq peaks of individual cell types  -------------
for (CELL_TYPE in CELL_TYPES[10]) {
  
  cat(paste0('\nLoading peaks for ', CELL_TYPE, ' ... \n'))
  peaks_df <- read_tsv(paste0(PEAK_DIR, CELL_TYPE,'.hg38.ext500bp.bed'), col_names = FALSE)
  colnames(peaks_df) <- c('chr', 'start', 'end', 'name', 'score', 'strand')

  cell_overlaps <- data.frame()
  
  cat(paste0('\nChecking for SNP overlaps in ', CELL_TYPE, ' ... \n'))
  for (i in 1:nrow(snps_no_patches)) {
    
    BASE_POSITION <- snps_no_patches$hg38_base_position[i]
    CHR <- snps_no_patches$chr_name[i]
    cat(paste0('\nSNP: ', 
               snps_no_patches$rsid[i], ', position: ',
               BASE_POSITION, ' ... \n'))
    
    overlaps <- filter(peaks_df, start <= BASE_POSITION, end >= BASE_POSITION, chr == paste0('chr', CHR))
    
    if (nrow(overlaps) > 0) { 
      
      overlaps <- cbind(overlaps, snps_no_patches[i,])
      print(overlaps) 
      
    }
    
    cell_overlaps <- rbind(cell_overlaps, overlaps)
    
  } 
  
  cat(paste0('\nAll ', nrow(snps_no_patches), ' SNPs checked in ', CELL_TYPE, ' ... \n'))
  
  cat(paste0('\nWriting overlapping SNPs to file ... \n'))
  write_tsv(cell_overlaps, paste0(PEAK_DIR, 'fine_mapped_SNPs/', CELL_TYPE,'_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext500bp.tsv'))
  assign(paste0(CELL_TYPE, '_SNP_overlaps'), cell_overlaps)
  
}


## Create binary count table to have all SNP/peak overlaps in one df  -----------------
# Create key of all SNPs overlapping peaks
cat('\nCreating key dataframe for all SNP/peak overlapping ... \n')
all_SNPs_key <- rbind(FC.ExN_SNP_overlaps, FC.InN_SNP_overlaps, FC.RG_SNP_overlaps, 
                      FC.MG_SNP_overlaps, FC.undef_SNP_overlaps, LGE.InN_SNP_overlaps, 
                      MGE.InN_SNP_overlaps, CGE.InN_SNP_overlaps, GE.RG_SNP_overlaps, 
                      GE.Proj_SNP_overlaps)
all_SNPs_key_df <- all_SNPs_key %>%
  arrange(rsid) %>%
  distinct(rsid)

# Create binary df - SNP in/not in peak for each cell
cat('\nCreating binary df for whether SNP is in/not in peak ... \n')
for (CELL_TYPE in CELL_TYPES) {
  
  cat(paste0('\nObtaining binary counts for: ', CELL_TYPE, ' ... \n'))
  # Test vector of rsIDs for cell type
  snp_test <- get(paste0(CELL_TYPE, '_SNP_overlaps')) %>%
    pull(rsid)
  
  if (exists('all_SNPs_binary_df')) {
    
    all_SNPs_binary_df <- all_SNPs_binary_df %>%
      rowwise() %>%
      mutate(!!CELL_TYPE := ifelse(rsid %in% snp_test, 1, 0)) %>%
      ungroup()
  
  } else {
    
    all_SNPs_binary_df <- all_SNPs_key_df %>%
      rowwise() %>%
      mutate(!!CELL_TYPE := ifelse(rsid %in% snp_test, 1, 0)) %>%
      ungroup()
    
  }
  
  
}

all_SNPs_binary_unique_rsids_df <- all_SNPs_binary_df %>%
  dplyr::select(-start, -end, -name, -score, -chr_name, -strand) %>%
  relocate(rsid) %>%
  distinct()

# Write table
cat('\nWriting binary count tables ... \n')
write_tsv(all_SNPs_binary_df, 
          paste0(PEAK_DIR, 
                 'fine_mapped_SNPs/', 
                 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext500bp.tsv'))

write_tsv(all_SNPs_binary_unique_rsids_df, 
          paste0(PEAK_DIR, 
                 'fine_mapped_SNPs/', 
                 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext500bp_unique_rsIDs.tsv'))

cat('Done.\n')



#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------