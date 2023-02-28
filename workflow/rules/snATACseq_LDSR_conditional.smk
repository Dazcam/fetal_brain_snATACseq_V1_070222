# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC conditional analyses on snATAC-seq data
#
#
# -------------------------------------------------------------------------------------

configfile: "../config/config.yaml"
from pandas import read_table

rule partitioned_heritability_conditional:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{key}_vs_{value}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.{key}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{key}_vs_{value}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.{value}.",
    message: "Running Prt Hrt conditional with {wildcards.key}_vs_{wildcards.value} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.{key}_vs_{value}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_conditional_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand("../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{sample[1][key]}_vs_{sample[1][value]}_{GWAS}_baseline.v1.2.results", sample=read_table('../resources/sheets/LDSR_conditions.tsv').iterrows(), GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_baseline.v1.2_summary_conditional_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR/part_herit/baseline_v1.2/conditional/",
             results_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/LDSR_conditions_for_summary.tsv"
    log:     "../results/logs/LDSR/snATACseq.{GWAS}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_FC.ExN_vs_FC.UNION_SCZ_baseline.v1.2.results > {output}
             
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 {params.ph_dir}snATACseq_LDSR_"$Line"_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_2/$Line/g" >> {output}
             done
             
             cp {output} ../results/LDSR/
             """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


