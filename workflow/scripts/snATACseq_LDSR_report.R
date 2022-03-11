#--------------------------------------------------------------------------------------
#
#    Create tables and plots sLDSC results - ATAC
#
#--------------------------------------------------------------------------------------

## Requirements  ----------------------------------------------------------------------

# Required on Hawk before opening R
# module load libgit2/1.1.0
# module load R/4.0.3

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

## Load packages  ---------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(rmarkdown)
library(argparser)

## Parse cell_type / set cell_type variable -------------------------------------------
cat('Parsing args ... \n')
p <- arg_parser("Create sLDSC Rmarkdown report for snATAC-seq data ... \n")
p <- add_argument(p, "markdown_file", help = "No snATACseq sLDSC Rmarkdown report file specified")
p <- add_argument(p, "out_dir", help = "No Rmarkdown html output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
p <- add_argument(p, "peak_extension", help = "No report filename specified")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
SUMMARY_FILE <- args$summary_file
MARKDOWN_FILE <- args$markdown_file
REPORT_FILE <- args$report_file
OUT_DIR <- args$out_dir
PEAK_EXT <- args$peak_extension

## Load Data  -------------------------------------------------------------------------
for (GWAS in c('SCZ', 'MDD', 'ADHD', 'ASD', 'HEIGHT', 'BPD')) {
  
  # Load LDSC summary
  cat(paste0('\nLoading data for', GWAS, ' ... \n'))
  gwas_df <- read_delim(paste0("../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_", PEAK_EXT, "_summary_", GWAS, ".tsv"), 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
  
  ## Split into seprate dfs  ------------------------------------------------------------
  cat('\nCreating plots and tables ... \n')
#  cer_df <- gwas_df %>% filter(grepl('cer', Category)) %>%
#    mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
#    mutate(neglog10 = -log10(Zp))
  fc_df <- gwas_df %>% filter(grepl('FC', Category)) %>%
    mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
    mutate(neglog10 = -log10(Zp))
  ge_df <- gwas_df %>% filter(grepl('GE', Category)) %>%
    mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
    mutate(neglog10 = -log10(Zp))
  
#  cer_enrich_plot <- ggplot(cer_df, aes(x=Category, y=Enrichment)) +
#    geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
#    scale_y_continuous(limits = c(0, 20), expand = c(0.02, 0)) +
#    theme_bw() + 
#    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
#          panel.grid.major = element_blank(), 
#          panel.grid.minor = element_blank(),
#          panel.border = element_rect(colour = "black", size = 0.8),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
#          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5, angle = 45),
#          axis.text.y  = element_text(colour="#000000", size=12)) +
#    ylab("Enrichment") +
#    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
#                  position=position_dodge(.9)) 
  
  fc_enrich_plot <- ggplot(fc_df, aes(x=Category, y=Enrichment)) +
    geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
    scale_y_continuous(limits = c(0, 25), expand = c(0.02, 0)) +
    theme_bw() + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour="#000000", size=12)) +
    ylab("Enrichment") +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                  position=position_dodge(.9)) 
  
  ge_enrich_plot <- ggplot(ge_df, aes(x=Category, y=Enrichment)) +
    geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
    scale_y_continuous(limits = c(0, 25), expand = c(0.02, 0)) +
    theme_bw() + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour="#000000", size=12)) +
    ylab("Enrichment") +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                  position=position_dodge(.9)) 
  
#  cer_pVal_plot <- ggplot(cer_df, aes(x=Category, y=neglog10)) +
#    geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
#    scale_y_continuous(limits = c(0, 12), expand = c(0.02, 0)) +
#    theme_bw() +
#    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
#          panel.grid.major = element_blank(), 
#          panel.grid.minor = element_blank(),
#          panel.border = element_rect(colour = "black", size = 0.8),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
#          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5, angle = 45),
#          axis.text.y  = element_text(colour="#000000", size=12)) +
#    ylab("-log10(Z-Score P)") 
  
  fc_pVal_plot <- ggplot(fc_df, aes(x=Category, y=neglog10)) +
    geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
    scale_y_continuous(limits = c(0, 25), expand = c(0.02, 0)) +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour="#000000", size=12)) +
    ylab("-log10(Z-Score P)") 
  
  ge_pVal_plot <- ggplot(ge_df, aes(x=Category, y=neglog10)) +
    geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85, angle = 45) +
    scale_y_continuous(limits = c(0, 25), expand = c(0.02, 0)) +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
          axis.text.y  = element_text(colour="#000000", size=12)) +
    ylab("-log10(Z-Score P)") 
  
  # Combine plots
#  cer_plot <- plot_grid(cer_enrich_plot, cer_pVal_plot, labels = c('A', 'B'), label_size = 16)
  fc_plot <- plot_grid(fc_enrich_plot, fc_pVal_plot, labels = c('A', 'B'), label_size = 16)
  ge_plot <- plot_grid(ge_enrich_plot, ge_pVal_plot, labels = c('A', 'B'), label_size = 16)

  # Assign dfs and plots
#  assign(paste0('cer_', GWAS, '_df'), cer_df)
  assign(paste0('fc_', GWAS, '_df'), fc_df)
  assign(paste0('ge_', GWAS, '_df'), ge_df)
#  assign(paste0('cer_', GWAS, '_plot'), cer_plot)
  assign(paste0('fc_', GWAS, '_plot'), fc_plot)
  assign(paste0('ge_', GWAS, '_plot'), ge_plot)

  
}


# Generate markdown document
cat('\nGenerating markdown document ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE,  output_dir = OUT_DIR)
cat('Done.')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


