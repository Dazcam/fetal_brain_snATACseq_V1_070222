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

rule format_annotation_data:
    input:   "../results/peaks/{CELL_TYPE}.hg19.ext250bp.bed"
    output:  "../results/garfield/anns/{CELL_TYPE}/chr22"
    params:  "../results/garfield/anns/{CELL_TYPE}/"
    message: "Transform annotation data into Garfield friendy format"
    log:     "../results/logs/garfield/format_annotation_data_{CELL_TYPE}.log"
    shell:
             "scripts/snATACseq_garfield_annotate_uk10k.sh  {input} {params} 2> {log}"

rule format_uk10k_SNP_data:
    input:   "../resources/garfield/garfield-data/maftssd/chr22"
    output:  "../results/garfield/uk10_SNPs/chr22"
    params:  "../results/garfield/uk10_SNPs/"
    message: "Transform uk10k SNP data into bedtools friendy format"
    log:     "../results/logs/garfield/format_uk10k_SNP_data.log"
    shell:
             "scripts/snATACseq_garfield_munge_UK10K_SNP_list.sh {params} 2> {log}"

