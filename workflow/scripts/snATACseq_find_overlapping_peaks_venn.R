#--------------------------------------------------------------------------------------
#
#    Find overlapping snATACseq peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# https://www.biostars.org/p/430015/
# https://www.biostars.org/p/453725/

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(ChIPpeakAnno)
library(cowplot)
library(VennDiagram)
#library(DiffBind)

##  Set variables  --------------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/"
PUBLIC_DATA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/public_datasets/"
REGIONS <- c('FC.ExN', 'FC.InN', 'FC.RG', 'MGE.InN', 'LGE.InN', 'CGE.InN', 'GE.RG')
BULK_REGIONS <- c('FC.UNION', 'GE.UNION')

## Functions  -------------------------------------------------------------------------
# Create a venn grob object from Chip peak anno output
pairwise_venn <- function(venn_counts, name1, name2) {
  
  dev.off()
  
  VENN_COUNTS <- venn_counts
  
  left <- unname(VENN_COUNTS$vennCounts[3,3])  + unname(VENN_COUNTS$vennCounts[4,3]) 
  right <- unname(VENN_COUNTS$vennCounts[4,3]) + unname(VENN_COUNTS$vennCounts[2,3])
  union <- unname(VENN_COUNTS$vennCounts[4,3])
  
  VENN <- draw.pairwise.venn(area1 = left, 
                             area2 = right, 
                             cross.area = union,
                             category = c(name1, name2),
                             scaled = FALSE,
                             fill = c("#0073C2FF", "#EFC000FF"),
                             fontfamily = "sans",
                             
                             # Numbers
                             cex = 2,
                             
                             # Set names
                             cat.cex = 1.5,
                             cat.fontface = "bold",
                             cat.default.pos = "outer",
                             cat.pos = c(-27, 27),
                             cat.dist = c(0.055, 0.055),
                             cat.fontfamily = "sans")
  
  return(VENN)
  
  
}

threeway_venn <- function(venn_counts, name1, name2, name3) {
  
  graphics.off()
  
  VENN_COUNTS <- as_tibble(venn_counts$vennCounts[])
  
  left <- sum(VENN_COUNTS[5:7, 5], VENN_COUNTS[8, 4])
  right <- sum(VENN_COUNTS[c(3, 4, 7), 6], VENN_COUNTS[8, 4])
  middle <- sum(VENN_COUNTS[c(2, 4, 6), 7], VENN_COUNTS[8, 4])
  intersect_12 <- sum(VENN_COUNTS[7, 4], VENN_COUNTS[8, 4])
  intersect_23 <- sum(VENN_COUNTS[4, 4], VENN_COUNTS[8, 4])
  intersect_13 <- sum(VENN_COUNTS[6, 4], VENN_COUNTS[8, 4])
  intersect_all3 <- as.double(VENN_COUNTS[8, 4])
  
  VENN <- draw.triple.venn(area1 = left, 
                           area2 = right, 
                           area3 = middle,
                           n12 = intersect_12,
                           n23 = intersect_23, 
                           n13 = intersect_13, 
                           n123 = intersect_all3,
                           category = c(name1, name2, name3),
                           fill = c("#0073C2FF", "#EFC000FF", '#21908DFF'),
                           fontfamily = "sans",
                           
                           # Numbers
                           cex = 2,
                           
                           # Set names
                           cat.cex = 2,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = c(-27, 27, 135),
                           cat.dist = c(0.055, 0.055, 0.085),
                           cat.fontfamily = "sans")
  
  return(VENN)
  
  
}



##  Load data  ------------------------------------------------------------------------
for (REGION in REGIONS) {
  
  PEAKS <- read_tsv(paste0(DATA_DIR, REGION, '.hg38.bed'), col_names = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3) %>%
    select(chr, start, end)
  PEAKS_GR <- toGRanges(PEAKS, format="BED", header=FALSE) 
  
  assign(paste0(REGION, '_peaks'), PEAKS)
  assign(paste0(REGION, '_gr'), PEAKS_GR)
  
}

for (REGION in BULK_REGIONS) {
  
  PEAKS <- read_tsv(paste0(DATA_DIR, REGION, '.hg19.bed'), col_names = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3) %>%
    select(chr, start, end) %>%
    filter(!grepl('_random|X|Y', chr))
  PEAKS_GR <- toGRanges(PEAKS, format="BED", header=FALSE) 
  
  assign(paste0(REGION, '_peaks'), PEAKS)
  assign(paste0(REGION, '_gr'), PEAKS_GR)
  
}



## Peak overlaps of FC and GE cell types  ---------------------------------------------
# Find peak overlaps
FC_overlaps <- findOverlapsOfPeaks(FC.ExN_gr, FC.InN_gr, FC.RG_gr, minoverlap = 100)
#GE_overlaps <- findOverlapsOfPeaks(MGE.InN_gr, CGE.InN_gr, LGE.InN_gr, GE.RG_gr, minoverlap = 100)
GE_overlaps_InN <- findOverlapsOfPeaks(MGE.InN_gr, CGE.InN_gr, LGE.InN_gr, minoverlap = 100)

# Create Venns
FC_venn <- makeVennDiagram(FC_overlaps, minoverlap = 100)
#GE_venn <- makeVennDiagram(GE_overlaps, minoverlap = 100) 
GE_InN_venn <- makeVennDiagram(GE_overlaps_InN, minoverlap = 100) 

GE_InN_venn2 <- threeway_venn(GE_InN_venn, 'MGE InN', 'CGE InN', 'LGE InN') 
FC_venn2 <- threeway_venn(FC_venn, 'FC ExN', 'FC InN', 'FC RG') 

plot_grid(FC_venn2, GE_InN_venn2, labels = 'AUTO', label_size = 20)


# GE UNION overlaps - Markenscoff 2020  -----------------------------------------------
# Note that this was done on hg19!!!!
mscoff_peaks <- readxl::read_excel(paste0(PUBLIC_DATA_DIR, 
                                          'markenscoff_2020/markenscoff_2020_supp_table_2.xlsx'), 
                   sheet = 'S2B') %>%
  select(`OCR coordinates (hg19)`) %>%
  dplyr::rename('peaks' = 'OCR coordinates (hg19)') %>%
  separate(col = peaks, sep = ':', c("chr", "region")) %>%
  separate(col = region, sep = '-', c("start", "end")) 
MSCOFF_PEAKS_GR <- toGRanges(mscoff_peaks, format="BED", header=FALSE) 

# Find peak overlaps
MSCOFF_PEAKS_overlaps <- findOverlapsOfPeaks(GE.UNION_gr, MSCOFF_PEAKS_GR, minoverlap = 100)

# Create Venns
MSCOFF_venn <- makeVennDiagram(MSCOFF_PEAKS_overlaps, minoverlap = 100, 
                               fill = c("#CC79A7", "#56B4E9"), 
                               col = c("#D55E00", "#0072B2"), 
                               cat.col = c("#D55E00", "#0072B2"),
                               main = "GE.UNION vs. Markenscoff bulk GE")

MSCOFF_venn2 <- pairwise_venn(MSCOFF_venn, 'GE union', 'Markenscoff-Papadimitriou GE') 


# FC UNION overlaps - Kouakou 2021 and de la Torre Ubieta 2019   ---------------------------
kouakou_peaks <- read_tsv(paste0(PUBLIC_DATA_DIR,'kouakou_2021/kouakou_bulk_fetal_FC_hg19.bed'), FALSE) %>%
  dplyr::rename(chr = X1, start = X2, end = X3)
KOUAKOU_PEAKS_GR <- toGRanges(kouakou_peaks, format="BED", header=FALSE) 

ubieta_peaks <- read_csv(paste0(PUBLIC_DATA_DIR,'de_la_Torre_Ubieta_2019/GSE95023_RAW/GSE95023_readswithinpeaks.csv'), TRUE) %>%
  dplyr::rename(chr = CHR, start = START, end = END) %>%
  dplyr::select(chr, start, end) %>%
  mutate(chr = paste0("chr", chr))
UBIETA_PEAKS_GR <- toGRanges(ubieta_peaks, format="BED", header=FALSE) 

# Find peak overlaps
KOUAKOU_overlaps <- findOverlapsOfPeaks(FC.UNION_gr, KOUAKOU_PEAKS_GR, minoverlap = 100)
UBIETA_overlaps <- findOverlapsOfPeaks(FC.UNION_gr, UBIETA_PEAKS_GR, minoverlap = 100)

# Create Venns
KOUAKOU_venn <- makeVennDiagram(KOUAKOU_overlaps, minoverlap = 100) 
UBIETA_venn <- makeVennDiagram(UBIETA_overlaps, minoverlap = 100)

KOUAKOU_venn2 <- pairwise_venn(KOUAKOU_venn, 'FC union', 'Kouakou FC') 
UBIETA_venn2 <- pairwise_venn(UBIETA_venn, 'FC union', 'de la Torre-Ubieta FC') 


plot_grid(UBIETA_venn2, KOUAKOU_venn2, MSCOFF_venn2, align = c("hv"), 
          labels = 'AUTO', label_size = 16)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
