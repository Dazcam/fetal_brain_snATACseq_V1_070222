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
    output:  "../results/QTL_for_Garfield/Hannon_mQTLs.tsv"
    message: "Preparing mQTL data"
    shell:
             "cut -d, -f2,3,10 {input} | sed 's/\"//g' | sed 's/,/\t /g' > {output}"

rule format_QTL_data:
    input:   "../results/QTL_for_Garfield/Hannon_mQTLs.tsv"
    output:  "../results/QTL_for_Garfield/Hannon_mQTL/chrX"
    params:  "../results/QTL_for_Garfield/Hannon_mQTL/"
    message: "Transform mQTL data into Garfield friendy format"
    log:    "../results/logs/garfield/format_QTL_data.log"
    shell:
             "scripts/snATACseq_create_input_QTL.sh {input} {params} 2>{log}" 
