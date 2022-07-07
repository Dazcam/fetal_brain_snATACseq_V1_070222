# -------------------------------------------------------------------------------------
#
#
#    Script for identifying overlaps between human accelerated regions (HARs) and 
#    snATACseq peaks
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule munge_HARs_file:
    input:   "../resources/sheets/HARs_capra_2019_supp_table_5.csv",
    output:  "../results/HARs/HARs_hg19.bed"
    message: "Munge {input}: add tab delims and remove header"
    log:     "../results/logs/HARs/munge_HARs_file.log"
    shell:
             """

             cat {input} | sed 's/,/\t/g' | tail -n +2 > {output} 2> {log}

             """

rule lift_over_HARs:
    input:   mybed = "../results/HARs/HARs_hg19.bed",
             chain_file = "../resources/liftover/hg19ToHg38.over.chain.gz"
    output:  "../results/HARs/HARs_hg38.bed"
    message: "Lifting {input.mybed} to hg38"
    log:     "../results/logs/HARs/munge_HARs_file.log"
    params:  "../results/HARs/HARs_hg38_unlifted.bed"
    shell:
             """

             ../resources/liftover/liftOver {input.mybed} {input.chain_file} {output} {params} 2> {log}

             """

rule sort_HARs:
    input:   "../results/HARs/HARs_hg38.bed"
    output:  "../results/HARs/HARs_hg38_srtd.bed"
    message: "Sorting {input}"
    log:     "../results/logs/HARs/sort_HARs.log"
    shell:
             """

             sort -k 1,1 -k2,2n {input} > {output} 2> {log}

             """

rule sort_beds:
    input:   "../results/peaks/{CELL_TYPE}.hg38.bed"
    output:  "../results/peaks/{CELL_TYPE}_hg38_srtd.bed"
    message: "Sorting {input}"
    log:     "../results/logs/HARs/sort_beds_{CELL_TYPE}.log"
    shell:
             """

             sort -k 1,1 -k2,2n {input} > {output}

             """

rule bedtools_fishers:
    input:   mybed = "../results/peaks/{CELL_TYPE}_hg38_srtd.bed",
             myHARs = "../results/HARs/HARs_hg38_srtd.bed"
    output:  "../results/HARs/{CELL_TYPE}_fishers.txt"
    message: "Running fishers exact test for overlap between {input.mybed} and {input.myHARs}"
    params:  "../resources/sheets/hg38_chrom_sizes.tsv"
    log:     "../results/logs/HARs/bedtools_fishers_{CELL_TYPE}.log"
    shell:
             """
             
             module load bedtools
             bedtools fisher -a {input.myHARs} -b {input.mybed} -g {params} -n 0.5 -r > {output} 2> {log} 

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
