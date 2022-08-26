#--------------------------------------------------------------------------------------
#
#    snATACseq HARs barplot
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
#library(ChIPpeakAnno)
library(cowplot)

##  Set Variables  --------------------------------------------------------------------
FISHERS_DIR <- '~/Desktop/fetal_brain_snATACseq_070222/results/HARs/'
HARs_DIR <- '~/Desktop/fetal_brain_snATACseq_070222/resources/public_datasets/girskis_2021/'
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/results/peaks/"
CELL_TYPES <- c("FC.ExN", "FC.InN", "FC.MG", "FC.RG", "FC.undef",
                "CGE.InN", "MGE.InN", "LGE.InN", "GE.RG", "GE.Undef")
  
##  Load Data  ------------------------------------------------------------------------
# Bedtools fisher tests
for (CELL_TYPE in CELL_TYPES) {
  
  fishers <- scan(paste0(FISHERS_DIR, CELL_TYPE, '_fishers.txt'), what="", sep="\n", skip = 13) %>%
    strsplit('\t') %>%
    unlist() %>%
    as.data.frame() %>%
    t() %>%
    as_tibble() %>%
    dplyr::rename("left" = 'V1', "right" = 'V2', "two_tailed" = 'V3', "ratio" = 'V4') %>%
    mutate(Cell_type = CELL_TYPE) %>%
    mutate(REGION = ifelse(str_detect(Cell_type, 'FC')  == TRUE, 'FC', 'GE'))
  
  assign(paste0(CELL_TYPE, '_fishers'), fishers)
  
}


# This is now obsolete as we did the overlaps in bedtools

# # HARs - Duplicate entry line 134 and 136 - this makes overlap function choke
# HARs <- read_tsv(paste0(HARs_DIR, 'GSE180714_HARs.bed')) %>%
#   filter(HAR_ID != 'HARsv2_0136')
# HARs_GR <- toGRanges(HARs, format="BED", header = FALSE, ) 
# 
# 
# # Cell specific peaks
# for (CELL_TYPE in CELL_TYPES) {
#   
#   PEAKS <- read_tsv(paste0(PEAK_DIR, CELL_TYPE, '.hg38.bed'), col_names = FALSE) %>%
#     dplyr::rename(chr = X1, start = X2, end = X3) %>%
#     select(chr, start, end)
#   PEAKS_GR <- toGRanges(PEAKS, format="BED", header = FALSE) 
#   
#   assign(paste0(CELL_TYPE, '_peaks'), PEAKS)
#   assign(paste0(CELL_TYPE, '_gr'), PEAKS_GR)
#   
# }


##  Create Fishers plot  --------------------------------------------------------------
fishers_grp_df_FC <- rbind(FC.ExN_fishers, FC.InN_fishers, FC.MG_fishers, FC.RG_fishers, 
                        FC.undef_fishers) %>% 
  mutate(log10p = -log10(as.numeric(two_tailed))) %>% 
  mutate(ratio = as.numeric(ratio)) 

fishers_grp_df_GE <- rbind(CGE.InN_fishers, MGE.InN_fishers, LGE.InN_fishers, 
                           GE.RG_fishers, GE.Undef_fishers) %>% 
  mutate(log10p = -log10(as.numeric(two_tailed))) %>% 
  mutate(ratio = as.numeric(ratio)) 
  


fishers_grp_df_FC$Cell_type <- factor(fishers_grp_df_FC$Cell_type, 
                               levels=c("FC.undef", "FC.RG", "FC.MG", "FC.InN", "FC.ExN"))
fishers_grp_df_GE$Cell_type  <- factor(fishers_grp_df_GE$Cell_type, 
                                       levels=c("GE.Undef", "GE.RG","CGE.InN","LGE.InN", "MGE.InN"))

FC_plot <- ggplot(fishers_grp_df_FC, aes(fill = log10p, y = Cell_type, x = ratio)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab('Ratio') +
  ggtitle("Frontal Cortex") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#075AFF",
                       high = "#FF0000")

GE_plot <- ggplot(fishers_grp_df_GE, aes(fill = log10p, y = Cell_type, x = ratio)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab('Ratio') +
  ggtitle("Ganglionic eminence") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#075AFF",
                       high = "#FF0000")

plot_grid(FC_plot, GE_plot, labels = 'AUTO', label_size = 16)

##  Create overlap tables  ------------------------------------------------------------
# Cell specific peaks - No longer need this as used bedtools to generate these tables
# for (CELL_TYPE in CELL_TYPES) {
#   
#   CELL_TYPE_GR <- get(paste0(CELL_TYPE, '_gr'))
#   
#   GRANGES <- findOverlapsOfPeaks(CELL_TYPE_GR, HARs_GR, minoverlap = 100)
#   assign(paste0(CELL_TYPE, '_overlaps_df'), GRANGES$overlappingPeaks[[1]])
#   
# }

# Need to tidy up plot 
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


         
