# Fetal brain snATACseq analyses (version 1)

Initial approach was to run FC and GE snATACseq data together using ArchR. 

***

Datasets for peak correspondence

+ [Markenscoff-Papadimitriou et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32610082/) - Bulk (e.g. FC, MGE, CGE, LGE)
 
+ [Song et al. 2020](https://pubmed.ncbi.nlm.nih.gov/33057195/) - FACS ExN, InN, RG and IP 
 
+ [Ziffra et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34616060/) - snATACseq from several regions of cortex (just test PFC) and MGE
 
+ [Trevino et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34390642/) - scATACseq from cortex (not subdivided) – do against our FC

***

Scripts

1. [snATACseq_cellRanger.smk](workflow/rules/snATACseq_cellRanger.smk) - Run Cell Ranger `atac-count` on fastQ files 

2. [snATACseq_processing.smk](workflow/rules/snATACseq_processing.smk) - snATACseq QC, clustering, visualisation, peak calling and downstream analysis in ArchR

    + [snATAC_QC.R](workflow/scripts/snATAC_QC.R) - Running 
    + [snATACseq_cluster_QC.R](workflow/scripts/snATACseq_cluster_QC.R) - Remove clusters with fewer than ...
    + [snATACseq_cluster_ID_geneExpMatrix.R](workflow/scripts/snATACseq_cluster_ID_geneExpMatrix.R) - Assign cluster identity using canonical marker genes
    + [snATACseq_unconstrained_integration.R](workflow/scripts/snATACseq_unconstrained_integration.R) - Run unconstrained integration of snATACseq and snRNAseq data
    + [snATACseq_pseudo-bulk-reps_and_peak_calling.R](workflow/scripts/snATACseq_pseudo-bulk-reps_and_peak_calling.R) - Create pseudo peak replicates and call peaks using default params
    + [snATACseq_call_peaks_ext_500bp.R](workflow/scripts/snATACseq_call_peaks_ext_500bp.R) - Call peaks with 500bp peak extension
    + [snATACseq_additional_analyses.R](workflow/scripts/snATACseq_additional_analyses.R) - Motif enrichment analysis, TF footprinting, co-accesibility analysis


+ snATACseq_map_PGC3_SCZ_index_and_proxy_SNPs_to_peaks.R
