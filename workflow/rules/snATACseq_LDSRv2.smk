# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq data peak files
#
#
# -------------------------------------------------------------------------------------

# Note that to get this to work you need to run to ldsr rule with LDSR_CELL_TYPES rule 
# in order to get annotation files and LD Scores for UNION peaks that will be used for
# conditional analysis. From partitioned_heritability rule use ATAC_CELL_TYPES

# Also used shell solution for if statement in the partitioned heritability rule as 
# Snakemake does not accept using conda env with run directive for some reason 

# See: https://github.com/snakemake/snakemake/issues/163

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

rule intersect_mybed:
    input:   annot = "../results/LDSR/annotation_files/baseline.{CHR}_no_head.bed",
             mybed = "../results/peaks/{CELL_TYPE}.hg19.noMHC.bed"
    output:  "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}.annot.gz"
    params:  out = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}.annot"
    message: "Creating LDSR annotation bed files for {wildcards.CHR}"
    shell:
             "module load bedtools; " 
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.out}; "
             "bedtools intersect -a {input.annot} -b {input.mybed} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.out}; "
             "gzip {params.out}"

rule ldsr:
    input:   annot = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase3 {wildcards.CHR} noMHC" 
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.noMHC.chr{CHR}_ldsr.log"
    shell:
             "python ../resources/ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule partitioned_heritability:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR_FC = expand("../results/LDSR/annotation_files/snATACseq.FC.UNION.noMHC.{CHR}.l2.ldscore.gz", CHR = range(1,23)),
             LDSR_GE = expand("../results/LDSR/annotation_files/snATACseq.GE.UNION.noMHC.{CHR}.l2.ldscore.gz", CHR = range(1,23)),
             LDSR_CELLS =  expand("../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.{CHR}.l2.ldscore.gz", CELL_TYPE = config['ATAC_CELL_TYPES'], CHR = range(1,23))
    output:  "../results/LDSR/part_herit/baseline_v1.2/{CELL_TYPE}_vs_UNION_noMHC_{GWAS}.done",
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE}.noMHC.",
             FC_anns = "../results/LDSR/annotation_files/snATACseq.FC.UNION.noMHC.",
             GE_anns = "../results/LDSR/annotation_files/snATACseq.GE.UNION.noMHC.",
             FC_out = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_FC.UNION_noMHC_{GWAS}_baseline.v1.2",
             GE_out = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE}_vs_GE.UNION_noMHC_{GWAS}_baseline.v1.2",
             out_dir = "../results/LDSR/part_herit/baseline_v1.2/"    
    message: "Running Prt Hrt conditional with {wildcards.CELL_TYPE}, the union peak sets and {wildcards.GWAS} GWAS, no MHC"
    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE}.vs.UNION.noMHC.{GWAS}_baseline.v1.2_partHerit.log"
    shell:
             """
             input_reads=({wildcards.CELL_TYPE})
             if [[ "$input_reads" == *"FC"* ]]; then
             
             echo $input_reads
             echo 'FC' 

             python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} --ref-ld-chr {params.baseline},{params.FC_anns},{params.LD_anns} \
             --overlap-annot --frqfile-chr {params.frqfile} --out {params.FC_out} --print-coefficients 2> {log} 
             
             touch {output}
 
             elif [[ "$input_reads" == *"GE"* ]]; then

       	     echo $input_reads
       	     echo 'GE'

             python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} --ref-ld-chr {params.baseline},{params.GE_anns},{params.LD_anns} \
             --overlap-annot --frqfile-chr {params.frqfile} --out {params.GE_out} --print-coefficients 2> {log}
             
             touch {output} 
             else
          
             echo error
             
             fi
             """

#    run:
#             if wildcards.CELL_TYPE in ['FC.ExN', 'FC.InN', 'FC.RG', 'FC.MG', 'FC.undef']: 

#                 shell("python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
#                 "--ref-ld-chr {params.baseline},{input.LDSR_FC},{params.LD_anns} --overlap-annot "
#                 "--frqfile-chr {params.frqfile} --out {params.FC_out} --print-coefficients 2> {log}")

#             else:

#                 shell("python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
#                 "--ref-ld-chr {params.baseline},{input.LDSR_GE},{params.LD_anns} --overlap-annot "
#                 "--frqfile-chr {params.frqfile} --out {params.GE_out} --print-coefficients 2> {log}")


rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv
    input:   expand("../results/LDSR/part_herit/baseline_v1.2/{CELL_TYPE}_vs_UNION_noMHC_{GWAS}.done", CELL_TYPE = config["ATAC_CELL_TYPES"], GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR/baseline_v1.2/snATACseq_LDSR_baseline.v1.2_summary_conditional_noMHC_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR/baseline_v1.2/",
             results_dir = "../results/LDSR/",
             cell_types = "../resources/sheets/atac_celltypes.tsv"
    log:     "../results/logs/LDSR/snATACseq.{GWAS}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_CGE.InN_vs_GE.UNION_noMHC_ADHD_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)

             for Line in $Lines; do

             grep L2_2 {params.ph_dir}snATACseq_LDSR_"$Line"_noMHC_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_2/$Line/g" >> {output}

             done

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
