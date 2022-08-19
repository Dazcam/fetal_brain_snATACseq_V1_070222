#--------------------------------------------------------------------------------------
#
#    snATACseq UMAPs group plots
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
library(cowplot)

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snATACseq_070222/results/archR_data_processing/rds_files/'

##  Load Data  ------------------------------------------------------------------------
FC_cluster_umap <- readRDS(paste0(DATA_DIR, 'FC_UMAP_clusters_broad.rds')) +
  ggtitle('Frontal Cortex') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
        axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 13),
        axis.text.x  = element_text(colour = "#000000", size = 13),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  scale_color_manual(values = c("#1F78B4", "#FFD92F", "#E31A1C", "#33A02C", "#FF7F00"), 
         labels = c("FC-ExN", "FC-InN", "FC-MG", "FC-RG", "FC-Undef"))


GE_cluster_umap <- readRDS(paste0(DATA_DIR, 'GE_UMAP_clusters_broad.rds')) +
  ggtitle('Ganglionic eminence') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
        axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 13),
        axis.text.x  = element_text(colour = "#000000", size = 13),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  scale_color_manual(values = c("#009F75", "#88C6ED", "#EA6A47", "#394BA0", "#D54799"), 
                     breaks = c("1-CGE-InN", "4-LGE-InN", "5-MGE-InN", "2-GE-RG", "3-GE-Undef"),
                     labels = c("CGE-InN", "LGE-InN", "MGE-InN", "GE-RG", "GE-Undef")) 



FC_gene_umaps <- readRDS(paste0(DATA_DIR, 'FC_UMAP_genes.rds'))
GE_gene_umaps <- readRDS(paste0(DATA_DIR, 'GE_UMAP_genes.rds'))

dput(RColorBrewer::brewer.pal(11, "BrBG"))

# Create lists
FC_plots <- list("Clusters" = FC_cluster_umap)
GE_plots <- list("Clusters" = GE_cluster_umap)

for (NAME in names(FC_gene_umaps)) {
  
  UMAP_PLOT <- FC_gene_umaps[[NAME]] + ggtitle(NAME) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
          axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
          axis.title.y = element_text(colour = "#000000", size = 13),
          axis.text.x  = element_text(colour = "#000000", size = 13),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 9)) +
    labs(fill = expression(Log[2]('NCounts+1')))  +
    scale_fill_gradientn(colours = c('#FAFAFA', "grey", "blue"))
  
  FC_plots[[NAME]] <- UMAP_PLOT
  
} 


for (NAME in names(GE_gene_umaps)) {
  
  UMAP_PLOT <- GE_gene_umaps[[NAME]] + ggtitle(NAME) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
          axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
          axis.title.y = element_text(colour = "#000000", size = 13),
          axis.text.x  = element_text(colour = "#000000", size = 13),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 9)) +
    labs(fill = expression(Log[2]('NormCounts+1'))) +
    scale_fill_gradientn(colours = c('#FAFAFA', "grey", "blue")) 
  
  GE_plots[[NAME]] <- UMAP_PLOT
  
} 

plot_grid(plotlist = FC_plots, align = 'hvrl', ncol = 2)
plot_grid(plotlist = GE_plots, align = 'hvrl', ncol = 2)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------