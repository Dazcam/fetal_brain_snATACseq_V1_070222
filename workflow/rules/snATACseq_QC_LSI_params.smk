# -------------------------------------------------------------------------------------
#
#
#    Script for snATACseq data. ArchR: QC testing LSI parameters 
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"


# -------------  RULES  ---------------

rule snATAC_seq_QC_LSI_params:
    input:  markdown = "scripts/snATACseq_QC_LSI_params.Rmd",
            qc_html = "../results/snATACseq_QC_{REGION}.html" # For rule order  
    output: "../results/QC_LSI_params/snATACseq_QC_{REGION}_LSI_{LSI_PARAM}.html"
    params: data_dir = "../resources/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/QC_LSI_params/",
            report_file = "snATACseq_QC_{REGION}_LSI_{LSI_PARAM}.html"
    log:    "../results/logs/QC_LSI_params/snATAC_QC_{REGION}_LSI_{LSI_PARAM}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_QC_LSI_params.R {wildcards.REGION} {wildcards.LSI_PARAM} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

            """

