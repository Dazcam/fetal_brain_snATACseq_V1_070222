#--------------------------------------------------------------------------------------
#
#    Create tables and plots sLDSC conditional results - snATACseq
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
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
SUMMARY_FILE <- args$summary_file
MARKDOWN_FILE <- args$markdown_file
OUT_DIR <- args$out_dir

## Load Data  -------------------------------------------------------------------------
for (CONDITION in c('fc.ExN', 'fc.RG', 'ge.RG')) {
  
  for (GWAS in c('SCZ', 'BPD', 'MDD', 'HEIGHT')) {
    
    # Load LDSC summary
    cat(paste0('\nLoading data for ', GWAS, ' ... \n'))
    gwas_df <- read_delim(paste0("results/snATACseq_LDSC_summary_", GWAS, ".vs.", CONDITION, ".tsv"), 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
    
    ## Split into seprate dfs  ------------------------------------------------------------
    cat('\nCreating plots and tables ... \n')
    df <- gwas_df %>% 
      mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
      mutate(neglog10 = -log10(Zp))
    
    enrich_plot <- ggplot(df, aes(x=Category, y=Enrichment)) +
      geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
      scale_y_continuous(limits = c(-4, 30), expand = c(0.02, 0)) +
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
    
    pVal_plot <- ggplot(df, aes(x=Category, y=neglog10)) +
      geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
      scale_y_continuous(limits = c(0, 16), expand = c(0.02, 0)) +
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
    

    # Combine plots
    combined_plot <- plot_grid(enrich_plot, pVal_plot, labels = c('A', 'B'), label_size = 16)
    
    
    # Assign dfs and plots
    assign(paste0('ldsc_', GWAS, '_vs.', CONDITION, '_df'), df)
    assign(paste0('ldsc_', GWAS, '_vs.', CONDITION, '_plot'), combined_plot)

    
    
  }
  
}


# Generate markdown document
cat('\nGenerating markdown document ... \n')
render(MARKDOWN_FILE, output_dir = OUT_DIR)
cat('Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
