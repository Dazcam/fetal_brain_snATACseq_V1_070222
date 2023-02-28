# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq UNION files
#
#
# -------------------------------------------------------------------------------------

# Note: Decided just to hard code FC and GE together rather than mess around with wildcard 
# constarints. Could try and integrate this by adding these cell types to 
# ATAC_CELL_TYPES later

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule lift_over:
    input:   FC = "../results/peaks/FC.UNION.hg38.bed",
             GE = "../results/peaks/GE.UNION.hg38.bed"
             chain_file = "../resources/liftover/hg38ToHg19.over.chain.gz"
    output:  FC = "../results/peaks/FC.UNION.hg19.bed",
             GE = "../results/peaks/GE.UNION.hg19.bed",
    message: "Lifting FC.UNION and GE.UNION to hg19"
    log:     "../results/logs/LDSR/scATACseq.UNION.hg19_liftOver.log"
    params:  FC = "../results/archR_data_processing/peaks/FC_UNION.hg38_unlifted.bed",
             GE = "../results/archR_data_processing/peaks/FC_UNION.hg38_unlifted.bed",
    shell:
             """

             ../resources/liftover/liftOver {input.FC} {input.chain_file} {output.FC} {params.FC} 2> {log}
             ../resources/liftover/liftOver {input.GE} {input.chain_file} {output.GE} {params.GE} 2> {log}

             """

rule intersect_mybed:
    input:   annot = "../results/LDSR/annotation_files/baseline.{CHR}_no_head.bed",
             FC = "../results/peaks/FC.UNION.hg19.bed",
             GE = "../results/peaks/GE.UNION.hg19.bed",
    output:  FC = "../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}.annot.gz",
             GE	= "../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}.annot.gz",
    params:  FC = "../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}.annot",
             GE = "../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}.annot",
    message: "Creating LDSR annotation bed files for {wildcards.CHR}"
    shell:
             "module load bedtools; " 
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.FC}; "
             "bedtools intersect -a {input.annot} -b {input.FC} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.FC}; "
             "gzip {params.FC}"
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.GE}; "
             "bedtools intersect -a {input.annot} -b {input.GE} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.GE}; "
             "gzip {params.GE}"
             

rule ldsr:
    input:   annot_FC = "../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}.annot.gz",
             annot_GE = "../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  FC = "../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}.l2.ldscore.gz",
             GE = "../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}.l2.ldscore.gz",
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores_FC = "../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}",
             ldscores_GE = "../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase3 {wildcards.CHR}" 
    log:     "../results/logs/LDSR/snATACseq.UNION.chr{CHR}_ldsr.log"
    shell:
             "python ../resources/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot_FC} --out {params.ldscores_FC} --print-snps {params.snps} 2> {log}"
             "python ../resources/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot_GE} --out {params.ldscores_GE} --print-snps {params.snps} 2> {log}" 


rule partitioned_heritability:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR_FC = expand("../results/LDSR/annotation_files/snATACseq.FC.UNION.{CHR}.annot.gz", CHR = range(1,23)),
             LDSR_GE = expand("../results/LDSR/annotation_files/snATACseq.GE.UNION.{CHR}.annot.gz", CHR = range(1,23))
    output:  FC = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_FC.UNION_{GWAS}_baseline.v1.2.results",
             GE = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_GE.UNION_{GWAS}_baseline.v1.2.results",
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.",
             FC_out = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_FC.UNION_{GWAS}_baseline.v1.2",
             GE_out = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_GE.UNION_{GWAS}_baseline.v1.2",
    message: "Running Prt Hrt conditional with {wildcards.CELL_TYPE}, the union peak sets and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.vs.UNION.{GWAS}_baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{input.LDSR_FC},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.FC_out} --print-coefficients 2> {log}"
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{input.LDSR_GE},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.GE_out} --print-coefficients 2> {log}" 


rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand("../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_{GWAS}_baseline.v1.2.results", CELL_TYPE = config["ATAC_CELL_TYPES"], GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_summary_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             results_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/atac_celltypes.tsv"
    log:     "../results/logs/LDSR/snATACseq.{GWAS}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_FC.ExN_SCZ_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.ph_dir}snATACseq_LDSR_"$Line"_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """

rule ldsr_mrkdwn_report:
    input:   summary = expand("../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_summary_{GWAS}.tsv", GWAS = config['LDSC_GWAS']),
             markdown = "scripts/snATACseq_LDSR_report.Rmd"
    output:  "../results/rmarkdown_reports/snATACseq_LDSR_baseline.v1.2_report.html"
    message: "Creating LDSR Rmarkdown report"
    params:  out_dir = "../results/rmarkdown_reports/",
             report_file = "snATACseq_LDSR_baseline.v1.2_report.html"
    log:     "../results/logs/LDSR/snATACseq.LDSR.baseline.v1.2.report.log"
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
