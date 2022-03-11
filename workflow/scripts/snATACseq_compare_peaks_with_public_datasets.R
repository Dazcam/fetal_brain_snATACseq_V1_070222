#--------------------------------------------------------------------------------------
#
#     Compare peaks - Test for similarity between peaks 
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# LDlinkR - requires a personal access token
# Github - https://github.com/CBIIT/LDlinkR
# Vingette: https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html

## Info  ------------------------------------------------------------------------------

#  Testing for similarity between peaks (ext 250bp) in our cell types and peaks from
#  celltypes in Ziffra et al (2021)

##  Load Packages  --------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(bedr)
library(openxlsx)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/archR_data_processing/peaks/"
ZIFFRA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/public_datasets/ziffra_2021/"
OUT_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peak_similarity/"
CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.RG", "FC.MG", "FC.undef", 
                "LGE.InN", "MGE.InN", "CGE.InN", "GE.RG", "GE.Proj")
ZIFFRA_CELL_TYPES = c("earlyEN", "dlEN", "ulEN", "IN_MGE", "RG",
           "Microglia", "EndoMural", "IN_MGE", "IN_MGE", "MGE",
           "IN_CGE", "RG", "Insula")
TOKEN <- "" # Needed for LDlinkR - saved in email

##  Load public data  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
cat('\nLoading public data ... \n')
peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
macs2_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype') 
specific_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST4 Specificpeaks_byCelltype') 

## Load snATACseq peaks - format for bedr chr:start-end  ------------------------------
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
cat('\nLoading snATACseq data and munging into bedr format ... \n')
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

## Set up pairwise comparison df  -----------------------------------------------------
# Create df for pairwise tests
cat('\nCreating pairwise df ... \n')
pairwise_df <- data.frame(bray  = c("FC.ExN", "FC.ExN", "FC.ExN", "FC.InN", "FC.RG", 
                                    "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "MGE.InN",
                                    "CGE.InN", "GE.RG", "GE.Proj"),
                          ziffra = c("earlyEN", "dlEN", "ulEN", "IN_MGE", "RG",
                                     "Microglia", "EndoMural", "IN_MGE", "IN_MGE", "MGE",
                                     "IN_CGE", "RG", "Insula")
)

## Munge public data into list of peaks per cell-type  --------------------------------
## Join key and cell type lists
cat('\nMunging Ziffra peaks into bedr format ... \n')
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

##  Sort snATACseq bed files for Jaccard testing  -------------------------------------
cat('\nSort snATACseq peaks ... \n')
for (CELL_TYPE in CELL_TYPES) {
  
  cat('\n\n Checking if ', CELL_TYPE, ' is valid and sorted\n\n')
  MY_PEAKS <- get(paste0(CELL_TYPE, "_peaks"))
  MY_PEAKS <- bedr.sort.region(MY_PEAKS) # Sorts the validates
  is.sorted <- is.sorted.region(MY_PEAKS, method = "lexicographical", engine = 'R')
  MY_PEAKS <- bedr.merge.region(MY_PEAKS)
  assign(paste0(CELL_TYPE, '_peaks_srtd'), MY_PEAKS)
  
}


##  Sort Ziffra bed files for Jaccard testing  -----------------------------------------
cat('\nSort and merge (where necessary) Ziffra peaks ... \n')
for (i in 1:length(ziffra_peaks_all)) {
  
  ZIFFRA_PEAKS <- ziffra_peaks_all[[i]]
  ZIFFRA_CELL_TYPE <- unique(ZIFFRA_PEAKS$cell_type)
  
  ZIFFRA_PEAKS <- ZIFFRA_PEAKS %>% 
    unite("peaks", chr:start, sep = ':') %>%
    unite("peaks", peaks:end, sep = '-') %>%
    select(peaks) %>%
    pull()
  
  # Sort and merge overlapping peaks
  cat('\n\n Checking if ', ZIFFRA_CELL_TYPE, ' is valid and sorted\n\n')
  ZIFFRA_PEAKS <- bedr.sort.region(ZIFFRA_PEAKS, method = "lexicographical", engine ='R')
  ZIFFRA_PEAKS <- bedr.merge.region(ZIFFRA_PEAKS)
  assign(paste0(ZIFFRA_CELL_TYPE, '_peaks_srtd'), ZIFFRA_PEAKS)
  
}


## Run pairwise tests listed in pairwise df  ------------------------------------------
cat('\nRun pairwise tests ... \n')
jaccard_df <- data.frame()

for (i in CELL_TYPES) {
  
  for (j in ZIFFRA_CELL_TYPES) {
    
    BRAY_cell_type <- i
    ZIFFRA_cell_type <- j
    BRAY_peaks <- get(paste0(i, '_peaks_srtd'))
    ZIFFRA_peaks <- get(paste0(j, '_peaks_srtd'))
    print(length(BRAY_peaks))
    print(length(ZIFFRA_peaks))
    
    cat('\nRunning Jaccard test for ', BRAY_cell_type, ' and ', ZIFFRA_cell_type, '\n')
    jaccard.stats <- jaccard(BRAY_peaks, ZIFFRA_peaks)
    rownames(jaccard.stats) <- NULL
    jaccard.stats <- unlist(unname(as.vector(jaccard.stats)))
    jaccard.stats <- c(BRAY_cell_type, ZIFFRA_cell_type, jaccard.stats)
    jaccard_df <- rbind(jaccard_df, jaccard.stats)
    
  }

}

colnames(jaccard_df) <- c('bray_cell', 'ziffra_cell', "intersection", 
                          "union-intersection", "jaccard", "n_intersections")


jaccard_list <- jaccard_df %>%
  arrange(desc(jaccard)) %>%
  group_split(bray_cell)
  
jaccard_list <- as.list(jaccard_list) # Need this write.xlsx doesn't work otherwise

# Write table 
cat('\nWriting Jaccard pairwise table\n')
dir.create(OUT_DIR)
write.xlsx(jaccard_list, paste0(OUT_DIR, "ziffra_peak_comparison_jaccard_tests.xlsx"))
          

cat('\nDone.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------




