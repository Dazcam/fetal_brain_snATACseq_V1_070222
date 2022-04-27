# -------------------------------------------------------------------------------------
#
#
#    Enrichment of fetal eQTL and mQTL in fetal snATACseq peaks using Garfield
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule prepare_QTL_data:
    input:   "../resources/garfield/All_Imputed_BonfSignificant_mQTLs.csv"
    output:  "../results/garfield/QTL_for_Garfield/Hannon_mQTLs.tsv"
    message: "Preparing mQTL data"
    shell:
             "cut -d, -f2,3,10 {input} | sed 's/\"//g' | sed 's/,/\t /g' > {output}"

rule format_QTL_data:
    input:   "../results/garfield/QTL_for_Garfield/Hannon_mQTLs.tsv"
    output:  "../results/garfield/QTL_for_Garfield/Hannon_mQTL/chr22"
    params:  "../results/garfield/QTL_for_Garfield/Hannon_mQTL/"
    message: "Transform mQTL data into Garfield friendy format"
    log:    "../results/logs/garfield/format_QTL_data.log"
    shell:
             "scripts/snATACseq_create_input_QTL.sh {input} {params} 2> {log}" 

rule format_uk10k_SNP_data:
    input:   "../resources/garfield/garfield-data/maftssd/chr22"
    output:  "../results/garfield/uk10_SNPs/chr22"
    params:  "../results/garfield/uk10_SNPs/"
    message: "Transform uk10k SNP data into bedtools friendy format"
    log:     "../results/logs/garfield/format_uk10k_SNP_data.log"
    shell:
             "scripts/snATACseq_garfield_munge_UK10K_SNP_list.sh {params} 2> {log}"

rule create_annotation_file:
    input:   bed_file = "../results/peaks/{CELL_TYPE}.hg19.ext250bp.bed",
             uk10k_snp_file = "../results/garfield/uk10_SNPs/chr22" 
    output:  "../results/garfield/anns/{CELL_TYPE}/chr22"
    params:  "../results/garfield/anns/{CELL_TYPE}/"
    message: "Create annotation file and munge into Garfield friendy format"
    log:     "../results/logs/garfield/create_annotation_file_{CELL_TYPE}.log"
    shell:
             """
             
             module load bedtools 
             scripts/snATACseq_garfield_create_annotation_file.sh  {input.bed_file} {params} 2> {log}

             """

rule run_garfield:
    input:   ann_file = "../results/garfield/anns/{CELL_TYPE}/chr22",
             uk10k_snp_file = "../results/garfield/uk10_SNPs/chr22"
    output:  "../results/garfield/output/{CELL_TYPE}/garfield.test.Hannon_mQTL.out"
    params:  annot_dir = "../results/garfield/anns/{CELL_TYPE}/",
             out_dir = "../results/garfield/output/{CELL_TYPE}/"
    message: "Run Garfield on to test for QTL enrichment in snATACseq peaks"
    log:     "../results/logs/garfield/run_garfield_{CELL_TYPE}.log"
    shell:
             "scripts/snATACseq_garfield.sh {params.annot_dir} {params.out_dir} 2> {log}"
