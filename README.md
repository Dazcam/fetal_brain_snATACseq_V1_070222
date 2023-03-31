# Fetal brain snATACseq analyses (version 1)

Initial approach was to run FC and GE snATACseq data together using ArchR. Much replicated code here that is optimised in version 2 (but some nice tid bits that may be useful in the future that we didn't use). 

***

Datasets for peak correspondence

+ [Markenscoff-Papadimitriou et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32610082/) - Bulk (e.g. FC, MGE, CGE, LGE)
 
+ [Song et al. 2020](https://pubmed.ncbi.nlm.nih.gov/33057195/) - FACS ExN, InN, RG and IP 
 
+ [Ziffra et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34616060/) - snATACseq from several regions of cortex (just test PFC) and MGE
 
+ [Trevino et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34390642/) - scATACseq from cortex (not subdivided) – do against our FC

***

Main scripts of interest

ArchR

1. [snATACseq_cellRanger.smk](workflow/rules/snATACseq_cellRanger.smk) - Run Cell Ranger `atac-count` on fastQ files 

2. [snATACseq_processing.smk](workflow/rules/snATACseq_processing.smk) - snATACseq QC, clustering, visualisation, peak calling and downstream analysis in ArchR

    + [snATACseq_QC.R](workflow/scripts/snATACseq_QC.R) - Running QCs and clustering up until unconstrained integration
    + [snATACseq_cluster_QC.R](workflow/scripts/snATACseq_cluster_QC.R) - Remove clusters that do not have >= 30 cells in at least 2 donors (decided to drop this as unnecessary)
    + [snATACseq_cluster_ID_geneExpMatrix.R](workflow/scripts/snATACseq_cluster_ID_geneExpMatrix.R) - Assign cluster identity using canonical marker genes
    + [snATACseq_unconstrained_integration.R](workflow/scripts/snATACseq_unconstrained_integration.R) - Run unconstrained integration of snATACseq and snRNAseq data
    + [snATACseq_pseudo-bulk-reps_and_peak_calling.R](workflow/scripts/snATACseq_pseudo-bulk-reps_and_peak_calling.R) - Create pseudo peak replicates and call peaks using default params
    + [snATACseq_call_peaks_ext_500bp.R](workflow/scripts/snATACseq_call_peaks_ext_500bp.R) - Call peaks with 500bp peak extension
    + [snATACseq_additional_analyses.R](workflow/scripts/snATACseq_additional_analyses.R) - Motif enrichment analysis, TF footprinting, co-accesibility analysis
 
 
 sLDSR
 
3.  [snATACseq_remove_MHC_regions_from_peaks.smk](workflow/rules/snATACseq_remove_MHC_regions_from_peaks) - Use bedtools to remove MHC peaks from snATACseq peak files
 
4. [snATACseq_munge_sumstats.smk](workflow/rules/snATACseq_munge_sumstats.smk) - Prepare 10 GWAS sumstats files for sLDSR

5. [snATACseq_LDSR.smk](workflow/rules/snATACseq_LDSR.smk) - Run sLD score regression on individual snATACseq cell type peaks and various GWAS sumstats

6. [snATACseq_LDSR_union_conditional.smk](workflow/rules/snATACseq_LDSR_union_conditional.smk) - Run sLD score regression on individual snATACseq cell type peaks and various GWAS sumstats adding UNION peaksets to baseline model

7. [snATACseq_LDSR_conditional.smk](workflow/rules/snATACseq_LDSR_conditional.smk) - Run sLD score regression on individual snATACseq cell type peaks and various GWAS sumstats adding UNION peaksets to baseline model. Alternative (and better) approach to 7, using key value pairs.

Peak file analyses

8. [snATACseq_HARs_overlaps.smk](workflow/rules/snATACseq_HARs_overlaps.smk) - Find overlaps between human accelerated regions (HARs) and snATACseq peaks.

9. [snATACseq_annotate_archR_coaccesible_peaks.R](scripts/snATACseq_annotate_archR_coaccesible_peaks.R) - Find overlaps of SCZ finemapped SNPs and snATACseq cA peaks.

10. [snATACseq_compare_peaks_with_public_datasets.R](scripts/snATACseq_compare_peaks_with_public_datasets.R) - Test jaccard similarity between our peaks and Ziffra et al. peaks.

11. [snATACseq_find_overlapping_peaks_venn.R](scripts/snATACseq_find_overlapping_peaks_venn.R) - Find overlaps of between our peaks and public dataset peaks (listed above).

12. [snATACseq_map_PGC3_SCZ_index_and_proxy_SNPs_to_peaks.R](scripts/snATACseq_map_PGC3_SCZ_index_and_proxy_SNPs_to_peaks.R) - Find overlaps of SCZ index SNPs (and proxies) and snATACseq peaks.

13. [snATACseq_map_PGC3_SCZ_finemapped_SNPs_to_peaks.R](scripts/snATACseq_map_PGC3_SCZ_finemapped_SNPs_to_peaks.R) - Find overlaps of SCZ finemapped SNPs and snATACseq peaks.



Additional exploratory scripts

1. [snATACseq_garfield.smk](workflow/rules/snATACseq_garfield.smk) - Testing Garfield on mQTL

    + [snATAC_create_input_QTL.R](workflow/scripts/snATACseq_create_input_QTL.sh) - Transform mQTL data into Garfield friendy format integration
    + [snATACseq_garfield_munge_UK10K_SNP_list.sh](workflow/scripts/snATACseq_garfield_munge_UK10K_SNP_list.sh) - Transform uk10k SNP data into bedtools friendy format
    + [snATACseq_garfield_create_annotation_file.sh](workflow/scripts/snATACseq_garfield_create_annotation_file.sh) - Create annotation file for each snATACseq cell type and munge into Garfield friendy format
    + [snATACseq_garfield.sh](workflow/scripts/snATACseq_garfield.sh) - Run Garfield
