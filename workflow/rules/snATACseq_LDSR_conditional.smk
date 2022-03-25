# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC conditional analyses  on snATAC-seq data
#
#
# -------------------------------------------------------------------------------------

rule partitioned_heritability_conditional_FC_ExN_vs_FC_InN:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.ExN_vs_FC.InN_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.ExN.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.ExN_vs_FC.InN_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.InN.{EXT}.",
    message: "Running Prt Hrt conditional with FC.ExN_vs_FC.InN, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.ExN_vs_FC.InN.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule partitioned_heritability_conditional_FC_ExN_vs_FC_RG:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.ExN_vs_FC.RG_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.ExN.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.ExN_vs_FC.RG_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.RG.{EXT}.",
    message: "Running Prt Hrt conditional with FC.ExN_vs_FC.RG, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.ExN_vs_FC.RG.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule partitioned_heritability_conditional_FC_InN_vs_FC_RG:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.InN_vs_FC.RG_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.InN.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.InN_vs_FC.RG_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.RG.{EXT}.",
    message: "Running Prt Hrt conditional with FC.InN_vs_FC.RG, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.InN_vs_FC.RG.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule partitioned_heritability_conditional_FC_InN_vs_FC_ExN:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.InN_vs_FC.ExN_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.InN.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.InN_vs_FC.ExN_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.ExN.{EXT}.",
    message: "Running Prt Hrt conditional with FC.InN_vs_FC.ExN, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.InN_vs_FC.ExN.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule partitioned_heritability_conditional_FC_RG_vs_FC_ExN:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.RG_vs_FC.ExN_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.RG.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.RG_vs_FC.ExN_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.ExN.{EXT}.",
    message: "Running Prt Hrt conditional with FC.RG_vs_FC.ExN, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.RG_vs_FC.ExN.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule partitioned_heritability_conditional_FC_RG_vs_FC_InN:
    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.RG_vs_FC.InN_{EXT}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.RG.{EXT}.",
             out_file = "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.RG_vs_FC.InN_{EXT}_{GWAS}_baseline.v1.2",
             condition = "../results/LDSR/annotation_files/snATACseq.FC.InN.{EXT}.",
    message: "Running Prt Hrt conditional with FC.RG_vs_FC.InN, {wildcards.EXT} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snATACseq.FC.RG_vs_FC.InN.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_conditional_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand("../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_FC.RG_vs_FC.InN_{EXT}_{GWAS}_baseline.v1.2.results", EXT = config['PEAK_EXTENSION'], GWAS = config['LDSC_GWAS'])
    output:  "../results/LDSR/part_herit/baseline_v1.2/conditional/snATACseq_LDSR_baseline.v1.2_{EXT}_summary_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.EXT} {wildcards.GWAS} GWAS"
    params:  ph_dir = "../results/LDSR/part_herit/baseline_v1.2/conditional/",
             results_dir = "../results/LDSR/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/LDSR_conditions.tsv"
    log:     "../results/logs/LDSR/snATACseq.{GWAS}_{EXT}_baseline_v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSR_FC.RG_vs_FC.InN_ext250bp_SCZ_baseline.v1.2.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 {params.ph_dir}snATACseq_LDSR_"$Line"_{wildcards.EXT}_{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """


#rule partitioned_heritability_conditional:
#    input:   SUMSTATS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
#             LDSR = expand("../results/LDSR/annotation_files/snATACseq.{CELL_TYPE_COND_A}.{EXT}.{CHR}.l2.ldscore.gz", CELL_TYPE_COND_A, = config['LDSR_CONDITIONAL'], EXT = config['PEAK_EXTENSION'], GWAS = config['LDSC_GWAS'])
#    output:  "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE_COND_A}.vs.{CELL_TYPE_COND_B}_{EXT}_{GWAS}_baseline.v1.2.results"
#    conda:   "../envs/ldsc.yml"
#    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
#             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
#             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
#             LD_anns = "../results/LDSR/annotation_files/snATACseq.FC.ExN.{CELL_TYPE_COND_A}.",
#             out_file = "../results/LDSR/part_herit/baseline_v1.2/snATACseq_LDSR_{CELL_TYPE_COND_A}.vs.{CELL_TYPE_COND_B}_{EXT}_{GWAS}_baseline.v1.2",
#             condition = "../results/LDSR/annotation_files/snATACseq.{CELL_TYPE_COND_B}.{EXT}.",
#    message: "Running Prt Hrt conditional with ExN.vs.FC.InN, {wildcards.EXT} and {wildcards.GWAS} GWAS"
#    log:     "../results/logs/LDSR/snATACseq.{CELL_TYPE_COND_A}.vs.{CELL_TYPE_COND_B}.{EXT}.{GWAS}_baseline.v1.2_partHerit_conditional.log"
#    shell:
#             "python ../resources/ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
#             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
#             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


