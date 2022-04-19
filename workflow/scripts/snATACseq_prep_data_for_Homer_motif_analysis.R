#--------------------------------------------------------------------------------------
#
#     Prepare data for Homer motif analysis
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# 1. Extract unique peaks that contained a PGC3 SCZ finemapped SNP
# 2. Extend (and extract) fine mapped PGC3 SCZ SNPs that fell under a peak by 15 bps either side


##  Load Packages  --------------------------------------------------------------------
library(tidyverse)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/"
PEAK_EXTENSION <- c("ext250bp")
REGIONS <- c("FC", "GE")


##  Load data  -------------------------------------------------------------------------
snp_peak_ovrlps <- read_tsv(paste0(PEAK_DIR, 'fine_mapped_SNPs/', 
                                   'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_', 
                                   PEAK_EXTENSION, '_with_peaks.tsv'))

##  Munge data  -------------------------------------------------------------------------
# Extract unique peaks that contained a PGC3 SCZ finemapped SNP
peaks <- snp_peak_ovrlps %>%
  select(chr, start, end) %>%
  distinct()

# Extend fine mapped PGC3 SCZ SNPs that fell under a peak by 15 bps either side
snps <- snp_peak_ovrlps %>%
  select(chr, rsid, hg38_base_position) %>%
  mutate(start = hg38_base_position - 15,
         end = hg38_base_position + 15) %>%
  relocate(chr, start, end, rsid) %>%
  distinct()

##  Save data  -------------------------------------------------------------------------
write_tsv(peaks, paste0(PEAK_DIR, 
                        'fine_mapped_SNPs/all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_', 
                        PEAK_EXTENSION, '_peaks_for_homer.tsv'))
write_tsv(snps, paste0(PEAK_DIR, 
                        'fine_mapped_SNPs/all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_', 
                       PEAK_EXTENSION, '_snps_plus_15bpext_for_homer.tsv'))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------