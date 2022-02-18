# -------------------------------------------------------------------------------------
#
#
#    Script for processing snATAC-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"


# -------------  RULES  ---------------

rule CR_cnt_ATAC:
    # Diminishing returns > 128Gs
    output: "../results/snATACseq_CR-atac_1.2.0/{SAMPLE}.stamp" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: FASTQ_DIR=config["FASTQ_DIR"],
            REFERENCE=config["REFERENCE_ATAC"]
    log:    "../results/logs/{SAMPLE}.log"
    shell:
            """
            cellranger-atac count --id={wildcards.SAMPLE} \
            --fastqs={params.FASTQ_DIR} \
            --sample={wildcards.SAMPLE} \
            --reference={params.REFERENCE} \
            --localcores=32 \
            --localmem=128 2> {log}
            
            touch {output}
            """

rule CR_aggr_ATAC:
    # Diminishing returns > 128Gs
#    input:  expand("../results/snATACseq_CR-atac_1.2.0/{SAMPLE}.stamp", SAMPLE=config["SAMPLES_ATAC"])
    output: "../results/snATACseq_CR-atac_1.2.0/{SAMPLE_AGGR}.aggr" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: FASTQ_DIR=config["FASTQ_DIR"],
            REFERENCE=config["REFERENCE_ATAC"]
    log:    "../results/logs/{SAMPLE_AGGR}.log"
    shell:
            """
            cellranger-atac aggr --id={wildcards.SAMPLE_AGGR} \
            --csv=../resources/sheets/snATACseq_pfc_aggr.csv \
            --normalize=none \
            --reference={params.REFERENCE} \
            --localcores=32 \
            --localmem=128 2> {log}

            touch {output}
            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
