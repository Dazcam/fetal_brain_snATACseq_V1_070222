# -------------------------------------------------------------------------------------
#
#    QC and processing snATAC-seq data
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule snATAC_seq_QC:
    input:  markdown = "scripts/snATACseq_QC.Rmd"
    output: "../results/rmarkdown_reports/snATACseq_QC_{REGION}.html"
    params: data_dir = "../results/snATACseq_CR-atac_1.2.0/", 
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/rmarkdown_reports/",
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

rule snATAC_seq_cluster_QC:
    input:  markdown = "scripts/snATACseq_cluster_QC.Rmd",
            html = "../results/rmarkdown_reports//snATACseq_QC_FC.html"
    output: "../results/rmarkdown_reports/snATACseq_cluster_QC_FC.html"
    params: data_dir = "../results/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/FC",
            report_dir = "../results/rmarkdown_reports/",
            report_file = "snATACseq_cluster_QC_FC.html"
    log:    "../results/logs/archR_data_processsing/snATAC_cluster_QC_FC.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_cluster_QC.R FC {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

            """

rule snATAC_seq_cluster_ID:
    input:  markdown = "scripts/snATACseq_cluster_ID_geneScoreMatrix.Rmd",
            html = "../results/rmarkdown_reports/snATACseq_cluster_QC_FC.html"
    output: "../results/rmarkdown_reports/snATACseq_cluster_ID_geneScoreMatrix_{REGION}.html"
    params: data_dir = "../results/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/rmarkdown_reports/",
            report_file = "snATACseq_cluster_ID_geneScoreMatrix_{REGION}.html",
            clustID_dir = "../results/archR_data_processing/cluster_ID/"
    log:    "../results/logs/archR_data_processing/snATACseq_cluster_ID_geneScoreMatrix_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_cluster_ID_geneScoreMatrix.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} {params.clustID_dir} 2> {log}

            
            """

rule snATAC_remove_batch_effects:
    input:  markdown = "scripts/snATACseq_remove_batch_effects.Rmd",
            html = "../results/rmarkdown_reports/snATACseq_cluster_QC_FC.html", # Needed for rule order
    output: "../results/rmarkdown_reports/snATACseq_remove_batch_effects_{REGION}.html"
    params: data_dir = "../results/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/rmarkdown_reports/",
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
            html = "../results/rmarkdown_reports/snATACseq_cluster_QC_FC.html" # Needed for rule order
    output: "../results/rmarkdown_reports/snATACseq_unconstrained_integration_{REGION}.html"
    params: data_dir = "../results/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            report_dir = "../results/rmarkdown_reports/",
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

#rule snATAC_constrained_integration:
#    input:  markdown ="scripts/snATACseq_constrained_integration.Rmd",
#            qc_html = "../results/snATACseq_unconstrained_integration_{REGION}.html" # Needed for rule order
#    output: "../results/snATACseq_constrained_integration_{REGION}.html"
#    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
#            archR_out_dir = "../results/ARCHR/{REGION}",
#            report_dir = "../results/",
#            report_file = "snATACseq_constrained_integration_{REGION}.html"
#    log:    "../results/logs/archR_data_processsing/snATAC_constrained_integration_{REGION}.log"
#    shell:
#            """
            
#            export R_LIBS_USER=/scratch/c.c1477909/R/library
#            module load libgit2/1.1.0
#            module load pandoc/2.7.3
#            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
#            scripts/snATACseq_constrained_integration.R {wildcards.REGION} {params.data_dir} \
#            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}

#            """

rule snATAC_pseudo_bulk_reps_and_peak_calling:
    input:  markdown = "scripts/snATACseq_pseudo-bulk-reps_and_peak_calling.Rmd",
            qc_html = "../results/rmarkdown_reports/snATACseq_unconstrained_integration_{REGION}.html" # Needed for rule order
    output: "../results/rmarkdown_reports/snATACseq_pseudo-bulk-reps_and_peak_calling_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            peaks_dir = "../results/peaks/",
            report_dir = "../results/rmarkdown_reports/",
            report_file = "snATACseq_pseudo-bulk-reps_and_peak_calling_{REGION}.html",
            macs2 = config['MACS2_PATH']
    log:    "../results/logs/archR_data_processing/snATAC_pseudo-bulk-reps_and_peak_calling_{REGION}.log"
    conda:  '../envs/macs2.yml'
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_pseudo-bulk-reps_and_peak_calling.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {params.peaks_dir} {input.markdown} {params.report_dir} {params.report_file} {params.macs2} 2> {log}

            """

rule snATACseq_call_peaks_ext_500bp:
    input:  markdown = "scripts/snATACseq_call_peaks_ext_500bp.Rmd",
            qc_html = "../results/rmarkdown_reports/snATACseq_unconstrained_integration_{REGION}.html" # Needed for rule order
    output: "../results/rmarkdown_reports/snATACseq_call_peaks_ext_500bp_{REGION}.html"
    params: data_dir = "/scratch/c.c1477909/snATACseq_CR-atac_1.2.0/",
            archR_out_dir = "../results/ARCHR/{REGION}",
            archR_ext_out_dir = "../results/ARCHR/{REGION}_peaks_ext500bp",
            peaks_dir = "../results/peaks/",
            report_dir = "../results/rmarkdown_reports/",
            report_file = "snATACseq_call_peaks_ext_500bp_{REGION}.html",
            macs2 = config['MACS2_PATH']
    log:    "../results/logs/archR_data_processing/snATACseq_call_peaks_ext_500bp_{REGION}.log"
    conda:  '../envs/macs2.yml'
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/snATACseq_call_peaks_ext_500bp.R {wildcards.REGION} {params.data_dir} {params.archR_out_dir} \
            {params.archR_ext_out_dir} {params.peaks_dir} {input.markdown} {params.report_dir} {params.report_file} {params.macs2} 2> {log}

            """

rule snATAC_additional_analyses:
    input:  markdown = "scripts/snATACseq_additional_analyses.Rmd",
            html1 = "../results/rmarkdown_reports/snATACseq_pseudo-bulk-reps_and_peak_calling_FC.html", # Needed for rule order
            html2 = "../results/rmarkdown_reports/snATACseq_pseudo-bulk-reps_and_peak_calling_GE.html",
    output: "../results/rmarkdown_reports/snATACseq_additional_analyses.html"
    params: rds_dir = "../results/archR_data_processing/rds_files/",
            report_dir = "../results/rmarkdown_reports/",
            report_file = "snATACseq_additional_analyses.html",
            hars_bed_file = config['HARs']
    log:    "../results/logs/archR_data_processing/snATAC_additional_analyses.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_additional_analyses.R {input.markdown} {params.rds_dir} \
            {params.report_dir} {params.report_file} {params.hars_bed_file} 2> {log}             

            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
