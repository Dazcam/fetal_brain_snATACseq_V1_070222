# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq data peak files
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule annot2bed:
    input:   folder = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3"
    params:  file = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.{CHR}.annot.gz"
    output:  "../results/LDSR_annotation_files/baseline.{CHR}_no_head.bed"
    message: "Preparing annotation files for {wildcards.CHR}"
    shell:
             "zcat {params.file} | tail -n +2 | awk -v OFS=\"\t\" '{{print \"chr\"$1, $2-1, $2, $3, $4}}' "
             "| sort -k1,1 -k2,2n > {output}"

rule lift_over:
    input:   mybed = "../results/archR_data_processing/peaks/{CELL_TYPE}.hg38.ext500bp.bed",
             chain_file = "../resources/liftover/hg38ToHg19.over.chain.gz"
    output:  "../results/archR_data_processing/peaks/{CELL_TYPE}.hg19.ext500bp.bed"
    message: "Lifting {input.mybed} to hg19"
    log:     "../results/logs/LDSR/scATACseq.{CELL_TYPE}_liftOver.log"
    params:  "../results/archR_data_processing/peaks/{CELL_TYPE}_hg38_ext500bp_unlifted.bed"
    shell:
             """

             ../resources/liftover/liftOver {input.mybed} {input.chain_file} {output} {params} 2> {log}

             """

rule intersect_mybed:
    input:   annot = rules.annot2bed.output,
             mybed = "../results/archR_data_processing/peaks/{CELL_TYPE}.hg19.ext500bp.bed"
    output:  "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.{CHR}.annot.gz"
    params:  out = "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.{CHR}.annot"
    message: "Creating LDSR annotation files for {wildcards.CHR}"
    shell:
             "module load bedtools; " 
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.out}; "
             "bedtools intersect -a {input.annot} -b {input.mybed} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.out}; "
             "gzip {params.out}"

rule ldsr:
    input:   annot = "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase3 {wildcards.CHR}" 
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.chr{CHR}_ldsc.log"
    shell:
             "python ../resources/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['ATAC_CELL_TYPES'], CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/snATACseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/snATACseq.{CELL_TYPE}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/snATACseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.{GWAS}_baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snATACseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2.results", CELL_TYPE = config["ATAC_CELL_TYPES"], GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snATACseq_LDSC_baseline.v1.2_summary_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR_part_herit/baseline_v1.2/",
             results_dir = "../results/LDSR_part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/atac_celltypes.tsv"
    log:     "../resultslogs/LDSR/snATACseq.{GWAS}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSC_FC.ExN_SCZ_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.ph_dir}snATACseq_LDSC_"$Line"_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """

rule ldsc_mrkdwn_report:
    input:   summary = expand("../results/LDSR_part_herit/baseline_v1.2/snATACseq_LDSC_baseline.v1.2_summary_{GWAS}.tsv", GWAS = config['LDSC_GWAS']),
             markdown = "scripts/snATACseq_LDSC_report.Rmd"
    output:  "../results/snATACseq_LDSC_baseline.v1.2_report.html"
    message: "Creating LDSC Rmarkdown report"
    params:  out_dir = "../results/",
             report_file = "snATACseq_LDSC_baseline.v1.2_report.html"
    log:     "../results/logs/LDSR/snATACseq.LDSC.baseline.v1.2.report.log"
    shell:
             """
             
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3 
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_LDSC_report.R {input.markdown} {params.out_dir} {params.report_file} 2> {log}

             """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

