# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC conditional analyses on snATAC-seq data
#
#
# -------------------------------------------------------------------------------------

LDSR_CONDITIONAL = {'FC.ExN': ["FC.InN", "FC.RG"],
                    'FC.InN': ["FC.ExN", "FC.RG"],
                    'FC.RG': ["FC.ExN", "FC.InN"]}

rule partitioned_heritability_conditional:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/SCZ_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{key}_vs_{value}_ext250bp_SCZ_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.{key}.ext250bp.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{key}_vs_{value}_ext250bp_SCZ_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.{value}.ext250bp.",
    message: "Running Prt Hrt conditional with {wildcards.key}_vs_{wildcards.value}, ext250bp and SCZ GWAS"
    log:     "../results/logs/LDSR/snATACseq.{key}_vs_{value}.ext250bp.SCZ_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_conditional_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   [f"../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_{key}_vs_{value}_ext250bp_SCZ_baseline.v1.2.results" for key, values in LDSR_CONDITIONAL.items() for value in values]
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_baseline.v1.2_ext250bp_summary_SCZ.tsv"
    message: "Creating summary file for ext250bp SCZ GWAS"
    params:  ph_dir = "../results/LDSR/part_herit/baseline_v1.2/conditional/",
             results_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/LDSR_conditions.tsv"
    log:     "../results/logs/LDSR/snATACseq.SCZ_ext250bp_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_FC.RG_vs_FC.InN_ext250bp_SCZ_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 {params.ph_dir}snATACseq_LDSR_"$Line"_ext250bp_SCZ_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


