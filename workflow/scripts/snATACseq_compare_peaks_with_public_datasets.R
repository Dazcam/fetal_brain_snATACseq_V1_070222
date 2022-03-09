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

library(Repitools)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/sheets/"
SNP_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/PGC3_snps/"
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

## Join key and cell type lists
ziffra_peaks_all <- macs2_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(peak_name, seqnames, start, end) %>%
  rename_with(~str_remove(., '_MACSpeaks')) %>%
  gather(cell_type, val, -peak_name, -seqnames, -start, -end) %>%
  filter(val == 1) %>%
  group_split(cell_type)

ziffra_peaks_specific <- specific_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(peak_name, seqnames, start, end) %>%
  rename_with(~str_remove(., '_Specificpeaks')) %>%
gather(cell_type, val, -peak_name, -seqnames, -start, -end) %>%
  filter(val == 1) %>%
  group_split(cell_type)
  


