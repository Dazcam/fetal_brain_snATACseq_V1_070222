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
    output:  "../results/LDSR/annotation_files/baseline.{CHR}_no_head.bed"
    message: "Preparing annotation files for {wildcards.CHR}"
    shell:
             "zcat {params.file} | tail -n +2 | awk -v OFS=\"\t\" '{{print \"chr\"$1, $2-1, $2, $3, $4}}' "
             "| sort -k1,1 -k2,2n > {output}"

rule lift_over:
    input:   mybed = "../results/archR_data_processing/peaks/{CELL_TYPE}.hg38.{EXT}.bed",
             chain_file = "../resources/liftover/hg38ToHg19.over.chain.gz"
    output:  "../results/archR_data_processing/peaks/{CELL_TYPE}.hg19.{EXT}.bed"
    message: "Lifting {input.mybed} to hg19"
    log:     "../results/logs/LDSR/scATACseq.{CELL_TYPE}.hg19.{EXT}_liftOver.log"
    params:  "../results/archR_data_processing/peaks/{CELL_TYPE}.hg38.{EXT}_unlifted.bed"
    shell:
             """

             ../resources/liftover/liftOver {input.mybed} {input.chain_file} {output} {params} 2> {log}

             """

rule intersect_mybed:
    input:   annot = rules.annot2bed.output,
             mybed = "../results/archR_data_processing/peaks/{CELL_TYPE}.hg19.{EXT}.bed"
    output:  "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.{CHR}.annot.gz"
    params:  out = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.{CHR}.annot"
    message: "Creating LDSR annotation bed files for {wildcards.CHR}"
    shell:
             "module load bedtools; " 
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.out}; "
             "bedtools intersect -a {input.annot} -b {input.mybed} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.out}; "
             "gzip {params.out}"

rule ldsr:
    input:   annot = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase3 {wildcards.CHR} {wildcards.EXT}" 
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.{EXT}.chr{CHR}_ldsr.log"
    shell:
             "python ../resources/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['ATAC_CELL_TYPES'], EXT = config['PEAK_EXTENSION'],  CHR = range(1,23))
    output:  "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_{EXT}_{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE}, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.{EXT}.{GWAS}_baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand("../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_{EXT}_{GWAS}_baseline.v1.2.results", CELL_TYPE = config["ATAC_CELL_TYPES"], EXT = config['PEAK_EXTENSION'], GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_{EXT}_summary_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.EXT} {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             results_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/atac_celltypes.tsv"
    log:     "../results/logs/LDSR/snATACseq.{GWAS}_{EXT}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_FC.ExN_ext250bp_SCZ_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.ph_dir}snATACseq_LDSR_"$Line"_{wildcards.EXT}_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """

rule ldsr_mrkdwn_report:
    input:   summary = expand("../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_{EXT}_summary_{GWAS}.tsv", EXT = config['PEAK_EXTENSION'], GWAS = config['LDSC_GWAS']),
             markdown = "scripts/snATACseq_LDSR_report.Rmd"
    output:  "../results/rmarkdown_reports/snATACseq_LDSR_baseline.v1.2_{EXT}_report.html"
    message: "Creating LDSR Rmarkdown report"
    params:  out_dir = "../results/rmarkdown_reports/",
             report_file = "snATACseq_LDSR_baseline.v1.2_{EXT}_report.html"
    log:     "../results/logs/LDSR/snATACseq.LDSR.baseline.v1.2.{EXT}.report.log"
    shell:
             """
             
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3 
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snATACseq_LDSR_report.R {input.markdown} {params.out_dir} {params.report_file} {wildcards.EXT} 2> {log}

             """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

