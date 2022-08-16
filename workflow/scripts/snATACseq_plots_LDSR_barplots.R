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
GWASs <- c('ASD', 'ADHD', 'BPD', 'MDD', 'SCZ', 'HEIGHT', 'BV', 'INSMNA', 'INTEL', 'NRTCSM')
  
##  Load Data  ------------------------------------------------------------------------
for (GWAS in GWASs) {
  
  SUMSTATS <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_baseline.v1.2_summary_', GWAS, '.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) %>%
    mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'FC', 'GE'))
  assign(paste0(GWAS, '_ldsr'), SUMSTATS)
  
  SUMSTATS_COND <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_baseline.v1.2_summary_conditional_', GWAS, '.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) %>%
    mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'FC', 'GE'))
  assign(paste0(GWAS, '_cond_ldsr'), SUMSTATS_COND)
  
}

ldsr_grp_df <- rbind(ADHD_ldsr, ASD_ldsr, BPD_ldsr, MDD_ldsr, SCZ_ldsr, 
      HEIGHT_ldsr, BV_ldsr, INSMNA_ldsr, INTEL_ldsr, NRTCSM_ldsr) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0))

ldsr_cond_grp_df <- rbind(ADHD_cond_ldsr, ASD_cond_ldsr, BPD_cond_ldsr, MDD_cond_ldsr, SCZ_cond_ldsr, 
                     HEIGHT_cond_ldsr, BV_cond_ldsr, INSMNA_cond_ldsr, INTEL_cond_ldsr, NRTCSM_cond_ldsr) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0))

ldsr_grp_df$Category <- factor(ldsr_grp_df$Category, 
                               levels=c("FC.ExN", "FC.InN", "FC.MG", "FC.RG", "FC.undef",
                                        "CGE.InN", "MGE.InN", "LGE.InN", "GE.RG", "GE.Undef"))
#ldsr_grp_df$Category  <- ldsr_grp_df$Category 
ldsr_grp_df$GWAS <- factor(ldsr_grp_df$GWAS, 
                           levels=c('NRTCSM', 'INTEL', 'INSMNA', 'BV', 'HEIGHT',  
                                     'SCZ', 'MDD', 'BPD', 'ADHD', 'ASD'))

ldsr_cond_grp_df$Category <- factor(ldsr_cond_grp_df$Category, 
                               levels=c("FC.ExN_vs_FC.UNION", "FC.InN_vs_FC.UNION", "FC.MG_vs_FC.UNION", "FC.RG_vs_FC.UNION", "FC.undef_vs_FC.UNION",
                                        "CGE.InN_vs_GE.UNION", "MGE.InN_vs_GE.UNION", "LGE.InN_vs_GE.UNION", "GE.RG_vs_GE.UNION", "GE.Undef_vs_GE.UNION"))
#ldsr_cond_grp_df$Category  <- ldsr_cond_grp_df$Category 
ldsr_cond_grp_df$GWAS <- factor(ldsr_cond_grp_df$GWAS, 
                                levels=c('NRTCSM', 'INTEL', 'INSMNA', 'BV', 'HEIGHT',  
                                'SCZ', 'MDD', 'BPD', 'ADHD', 'ASD'))

for (DF in c('ldsr_grp', 'ldsr_cond_grp')) {
  
 PLOT <-  ggplot(get(paste0(DF, '_df')), aes(fill = Category, y = LDSR, x = GWAS)) + 
    geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity") +
    theme_bw() +
    coord_flip() +
    facet_wrap(~REGION) +
    ylab(expression(-log[10](P))) +
    geom_hline(yintercept=-log10(0.05/100), linetype = "dashed", color = "black") +
    geom_hline(yintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13))
  
  assign(paste0(DF, '_plot'), PLOT)
  
}

plot_grid(ldsr_grp_plot, ldsr_cond_grp_plot)


# Need to tidy up plot 
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


