#--------------------------------------------------------------------------------------
#
#     Compare peaks 
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# LDlinkR - requires a personal access token
# Github - https://github.com/CBIIT/LDlinkR
# Vingette: https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html

## Info  ------------------------------------------------------------------------------



##  Load Packages  --------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(bedr)
#library(Repitools)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/archR_data_processing/peaks/"
ZIFFRA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/public_datasets/ziffra_2021/"
REGIONS <- c("FC", "GE")
CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
ZIFFRA_CELL_TYPES <- c()
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Load public data  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
macs2_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype') 
specific_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST4 Specificpeaks_byCelltype') 

## Load snATACseq peaks - format for bedr chr:start-end
options(scipen = 999) # required to prevent number being abbr. in scientific notation
for (CELL_TYPE in CELL_TYPES) {
  
  PEAKS <- read_delim(paste0(PEAK_DIR, CELL_TYPE, ".hg38.bed"), 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE, col_names = FALSE)  %>%
    unite("peaks", X1:X2, sep = ':') %>%
    unite("peaks", peaks:X3, sep = '-') %>%
    select(peaks) %>%
    pull
  assign(paste0(CELL_TYPE, '_peaks'), PEAKS)
  
}

## Munge public data into list of peaks per cell-type
## Join key and cell type lists
ziffra_peaks_all <- macs2_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(seqnames, start, end, peak_name) %>%
  rename(chr = seqnames)  %>%
  rename_with(~str_remove(., '_MACSpeaks')) %>%
  gather(cell_type, val, -peak_name, -chr, -start, -end) %>%
  filter(val == 1) %>%
  group_split(cell_type) 

ziffra_peaks_specific <- specific_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(seqnames, start, end, peak_name) %>%
  rename_with(~str_remove(., '_Specificpeaks')) %>%
  gather(cell_type, val, -peak_name, -seqnames, -start, -end) %>%
  filter(val == 1) %>%
  group_split(cell_type)





# 
zif_test <- ziffra_peaks_all[1:50,] %>%
  unite("peaks", chr:start, sep = ':') %>%
  unite("peaks", peaks:end, sep = '-') %>%
  select(peaks) %>%
  pull()

zif_test <- bedr.sort.region(zif_test, method = "lexicographical", engine ='R')
zif_test <- bedr.merge.region(zif_test)

# Note efficient - pull my cell type out of loop as it's running the checks for those multiple times
for (i in 1:length(ziffra_peaks_all)) {
  
  ZIFFRA_PEAKS <- ziffra_peaks_all[[i]]
  ZIFFRA_CELL_TYPE <- unique(PEAKS$cell_type)
  
  ZIFFRA_PEAKS <- ZIFFRA_PEAKS %>% 
    unite("peaks", chr:start, sep = ':') %>%
    unite("peaks", peaks:end, sep = '-') %>%
    select(peaks) %>%
    pull()
  
  # Sort and merge overlapping peaks
  cat('\n\n Checking if ', ZIFFRA_CELL_TYPE, ' is valid and sorted\n\n')
  ZIFFRA_PEAKS <- bedr.sort.region(ZIFFRA_PEAKS, method = "lexicographical", engine ='R')
  ZIFFRA_PEAKS <- bedr.merge.region(ZIFFRA_PEAKS)
  
  for (CELL_TYPE in CELL_TYPES) {
    
    cat('\n\n Checking if ', CELL_TYPE, ' is valid and sorted\n\n')
    MY_PEAKS <- get(paste0(CELL_TYPE, "_peaks"))
    MY_PEAKS <- bedr.sort.region(MY_PEAKS) # Sorts the validates
    is.sorted <- is.sorted.region(MY_PEAKS, method = "lexicographical", engine = 'R')
    
    cat('\n\n Running jaccard on  ', CELL_TYPE, ' and ', ZIFFRA_CELL_TYPE,' \n\n')
    jaccard.stats <- jaccard(MY_PEAKS, ZIFFRA_PEAKS)
    
    cat('\n\n Assigning files ... \n\n')
    assign(paste0(CELL_TYPE, '_peaks_srtd'), PEAKS)
    assign(paste0(CELL_TYPE, '_', ZIFFRA_CELL_TYPE, '_jaccard'), jaccard.stats)
    
  }

}




reldist.stats <- reldist(FC.ExN_peaks, zif_test)

