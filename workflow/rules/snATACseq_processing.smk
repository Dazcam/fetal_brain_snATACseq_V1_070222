# -------------------------------------------------------------------------------------
#
#    QC and processing snATAC-seq data
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# ----------  SET VARIABLES  ----------

#ROBJ_DIR = SCRATCH + "results/R_objects/"
#RESULTS_DIR = SCRATCH + "results/reports/ATAC/"
#PEAKS_DIR = SCRATCH + "results/bed_files_for_LDSC_ATAC/"

# -------------  RULES  ---------------

rule snATAC_seq_QC:
    input:  markdown = "scripts/snATACseq_QC.Rmd"
    output: "../results/snATACseq_QC_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/", 
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/",
            report_file = "snATACseq_QC_{REGION}.html"
    log:    "../results/logs/archR_data_processsing/snATAC_QC_{REGION}.log"
    shell:
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_QC.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
            
            """

#rule snATAC_seq_cluster_QC:
#    input:  markdown = MARKDOWN_DIR + "snATACseq_cluster_QC.Rmd",
#            html = RESULTS_DIR + "snATACseq_QC_Cer.html"
#    output: RESULTS_DIR + "snATACseq_cluster_QC_Cer.html"
#    params: data_dir = ATAC_DATA_DIR,
#            archR_out_dir = RESULTS_DIR + "ARCHR/Cer",
#            report_dir = RESULTS_DIR,
#            report_file = "snATACseq_cluster_QC_Cer.html"
#    log:    LOG_DIR + "snATAC_cluster_QC_Cer.log"
#    shell:
#            """
#
#            export R_LIBS_USER=/scratch/c.c1477909/R/library
#            module load libgit2/1.1.0
#            module load pandoc/2.7.3
#            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
#            scripts/R/snATACseq_cluster_QC.R Cer {params.data_dir} \
#            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

#            """

rule snATAC_remove_batch_effects:
    input:  markdown = "scripts/snATACseq_remove_batch_effects.Rmd",
            html = "../results/snATACseq_QC_{REGION}.html", # Needed for rule order
    output: "../results/snATACseq_remove_batch_effects_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/",
            report_file = "snATACseq_remove_batch_effects_{REGION}.html"
    log:    "../results/logs/archR_data_processsing/snATAC_remove_batch_effects_{REGION}.log"
    shell:
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_remove_batch_effects.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
            
            """

rule snATAC_unconstrained_integration:
    input:  markdown = "scripts/snATACseq_unconstrained_integration.Rmd",
            html = "../results/snATACseq_remove_batch_effects_{REGION}.html" # Needed for rule order
    output: "../results/snATACseq_unconstrained_integration_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/",
            report_file = "snATACseq_unconstrained_integration_{REGION}.html"
    log:    "../results/logs/archR_data_processsing/snATAC_unconstrained_integration_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
	    scripts/snATACseq_unconstrained_integration.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

            """

rule snATAC_constrained_integration:
    input:  markdown ="scripts/snATACseq_constrained_integration.Rmd",
            qc_html = "../results/snATACseq_unconstrained_integration_{REGION}.html" # Needed for rule order
    output: "../results/snATACseq_constrained_integration_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/",
            report_file = "snATACseq_constrained_integration_{REGION}.html"
    log:    "../results/logs/archR_data_processsing/snATAC_constrained_integration_{REGION}.log"
    shell:
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_constrained_integration.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

            """

#rule snATAC_pseudo_bulk_reps_and_peak_calling:
#    input:  markdown = MARKDOWN_DIR + "snATACseq_pseudo-bulk-reps_and_peak_calling.Rmd",
#            qc_html = RESULTS_DIR + "snATACseq_constrained_integration_{REGION}.html" # Needed for rule order
#    output: RESULTS_DIR + "snATACseq_pseudo-bulk-reps_and_peak_calling_{REGION}.html"
#    params: data_dir = ATAC_DATA_DIR,
#            archR_out_dir = RESULTS_DIR + "ARCHR/{REGION}",
#            peaks_dir = PEAKS_DIR,
#            report_dir = RESULTS_DIR,
#            report_file = "snATACseq_pseudo-bulk-reps_and_peak_calling_{REGION}.html",
#            macs2 = config['MACS2_PATH']
#    log:    LOG_DIR + "snATAC_pseudo-bulk-reps_and_peak_calling_{REGION}.log"
#    conda:  'envs/macs2.yml'
#    shell:
#            """

#            export R_LIBS_USER=/scratch/c.c1477909/R/library
#            module load libgit2/1.1.0
#            module load pandoc/2.7.3
#            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
#            scripts/R/snATACseq_pseudo-bulk-reps_and_peak_calling.R {wildcards.REGION} {params.data_dir} \
#            {params.archR_out_dir} {params.peaks_dir} {input.markdown} {params.report_dir} {params.report_file} {params.macs2} 2> {log}

#            """

#rule snATAC_additional_analyses:
#    input:  markdown = MARKDOWN_DIR + "snATACseq_additional_analyses.Rmd",
#            html1 = RESULTS_DIR + "snATACseq_pseudo-bulk-reps_and_peak_calling_FC.html", # Needed for rule order
#            html2 = RESULTS_DIR + "snATACseq_pseudo-bulk-reps_and_peak_calling_GE.html",
#    output: RESULTS_DIR + "snATACseq_additional_analyses.html"
#    params: report_dir = RESULTS_DIR,
#            report_file = "snATACseq_additional_analyses.html",
#    log:    LOG_DIR + "snATAC_additional_analyses.log"
#    shell:
#            """

#            export R_LIBS_USER=/scratch/c.c1477909/R/library
#            module load libgit2/1.1.0
#            module load pandoc/2.7.3
#            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
#            scripts/R/snATACseq_additional_analyses.R {input.markdown} \
#            {params.report_dir} {params.report_file} 2> {log}             

#            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
