# -------------------------------------------------------------------------------------
#
#
#    Script for identifying overlaps between human accelerated regions (HARs) and 
#    snATACseq peaks
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"
localrules: bedtools_fishers, bedtools_intersect, munge_overlap_df

## FIRST TWO RULES ONLY REQUIRED IF USING CAPRA HG19 FILE

# -------------  RULES  ---------------
#rule munge_HARs_file:
#    input:   "../resources/sheets/HARs_capra_2019_supp_table_5.csv",
#    output:  "../results/HARs/HARs_hg19.bed"
#    message: "Munge {input}: add tab delims and remove header"
#    log:     "../results/logs/HARs/munge_HARs_file.log"
#    shell:
#             """

#             cat {input} | sed 's/,/\t/g' | tail -n +2 > {output} 2> {log}

#             """

#rule lift_over_HARs:
#    input:   mybed = "../results/HARs/HARs_hg19.bed",
#             chain_file = "../resources/liftover/hg19ToHg38.over.chain.gz"
#    output:  "../results/HARs/HARs_hg38.bed"
#    message: "Lifting {input.mybed} to hg38"
#    log:     "../results/logs/HARs/munge_HARs_file.log"
#    params:  "../results/HARs/HARs_hg38_unlifted.bed"
#    shell:
#             """

#             ../resources/liftover/liftOver {input.mybed} {input.chain_file} {output} {params} 2> {log}

#             """

rule sort_HARs:
    input:   "../resources/public_datasets/girskis_2021/GSE180714_HARs.bed"
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
             bedtools fisher -a {input.myHARs} -b {input.mybed} -g {params} > {output} 2> {log} 

             """

rule bedtools_intersect:
    input:   mybed = "../results/peaks/{CELL_TYPE}_hg38_srtd.bed",
             myHARs = "../results/HARs/HARs_hg38_srtd.bed"
    output:  "../results/HARs/{CELL_TYPE}_overlaps.txt"
    message: "Running fishers intersect to calculate bp overlap between {input.mybed} and {input.myHARs}"
    log:     "../results/logs/HARs/bedtools_intersect_{CELL_TYPE}.log"
    shell:
             """

             module load bedtools
             bedtools intersect -a {input.myHARs} -b {input.mybed} -wo > {output} 2> {log}

             """

rule munge_overlap_df:
    input:   "../results/HARs/{CELL_TYPE}_overlaps.txt"
    output:  "../results/HARs/{CELL_TYPE}_overlaps_munged.txt"
    message: "Removing cols from overlap df for {wildcards.CELL_TYPE}"
    log:     "../results/logs/HARs/munge_overlap_df_{CELL_TYPE}.log"
    shell:
             """

             echo $'HARSID\tchr\tstart\tend\tnearest_gene\tchrom\tchrom_start\tchrom_end\toverlap_bp' > temp.{wildcards.CELL_TYPE}
             cat {input} | cut -f 4,1,2,3,15-18,22 |\
             awk '{{print $4"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' >> temp.{wildcards.CELL_TYPE}
             mv temp.{wildcards.CELL_TYPE} {output} 

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
