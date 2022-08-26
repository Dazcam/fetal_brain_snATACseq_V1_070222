#--------------------------------------------------------------------------------------
#
#    snATACseq LDSR barcharts
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
library(viridis)
library(cowplot)

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snATACseq_070222/results/LDSR/'
GWASs <- c('ASD', 'ADHD', 'BPD', 'MDD', 'SCZ', 'BV', 'INSMNA', 'INTEL', 'NRTCSM')
  
##  Load Data  ------------------------------------------------------------------------
for (GWAS in GWASs) {
  
  # SUMSTATS <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_baseline.v1.2_summary_', GWAS, '.tsv')) %>%
  #   mutate(GWAS = rep(GWAS, nrow(.))) %>%
  #   mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'Frontal Cortex', 'Ganglionic eminence'))
  # assign(paste0(GWAS, '_ldsr'), SUMSTATS)
  
  SUMSTATS_COND <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_baseline.v1.2_summary_conditional_noMHC_', GWAS, '.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) %>%
    mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'Frontal Cortex', 'Ganglionic eminence'))
  assign(paste0(GWAS, '_cond_ldsr'), SUMSTATS_COND)
  
}

ldsr_cond_grp_df <- rbind(ADHD_cond_ldsr, ASD_cond_ldsr, BPD_cond_ldsr, MDD_cond_ldsr, SCZ_cond_ldsr, 
                     BV_cond_ldsr, INSMNA_cond_ldsr, INTEL_cond_ldsr, NRTCSM_cond_ldsr) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  separate(Category, c('cell_type', 'suffix'), sep = '_', remove = FALSE, extra = 'drop') %>%
  mutate(across('GWAS', str_replace, 'ASD', 'AUT')) 

ldsr_cond_grp_df$cell_type <- factor(ldsr_cond_grp_df$cell_type, 
                                     levels=c("FC.ExN", "FC.InN", "FC.MG", "FC.RG", "FC.undef",
                                              "CGE.InN", "MGE.InN", "LGE.InN", "GE.RG", "GE.Undef"))

ldsr_cond_grp_df$GWAS <- factor(ldsr_cond_grp_df$GWAS, 
                                levels=c('ADHD', 'AUT', 'BPD', 'MDD', 'SCZ', 'BV', 'INSMNA', 'INTEL', 'NRTCSM'))

ldsr_plot <- ggplot(ldsr_cond_grp_df, aes(fill = cell_type, y = LDSR, x = GWAS)) + 
  geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity") +
  theme_bw() +
  coord_flip() +
  facet_grid(cols = vars(REGION), rows = vars(GWAS), scales = 'free_y') +
  ylab(expression(-log[10](P))) +
  xlab(NULL) +
  geom_hline(yintercept=-log10(0.05/90), linetype = "dashed", color = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines")) +
  scale_fill_manual(values = c('#6098ab','#f18e2a', '#e1575a', '#75b7b2', '#58a14e',
                               '#edc949', '#b07aa1', '#ff9ca7', '#9c755f', '#bab0ab'))


# Color scheme - to match UMAPs
# scale_fill_manual(values = c("#1F78B4", "#FFD92F", "#E31A1C", "#33A02C", "#FF7F00",
#                              "#009F75", "#88C6ED", "#EA6A47", "#394BA0", "#D54799")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


