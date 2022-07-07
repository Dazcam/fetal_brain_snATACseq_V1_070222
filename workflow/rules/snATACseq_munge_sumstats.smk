# -------------------------------------------------------------------------------------
#
#    Script to standarise and munge GWAS sumstats for LDSC analysis
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  --------------
rule standardise_phenotype_sumstats:
    input:   lambda wildcards: config['LDSC_GWAS'][wildcards.GWAS]
    output:  "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_sumstats_fixHeader.tsv.gz"
    message: "Fixing header for {input} (If also AD adding rsIDs)"
    params:  out_dir = "../results/gwas_sumstats/GWAS_sumstats_standardised/",
             snps = "../resources/gwas_sumstats/00-All.vcf.gz"     
    log:     "../results/logs/munge_sumstats/{GWAS}_standardise_sumstats.log"
    benchmark: "../results/logs/benchmarks/{GWAS}_standardise_sumstats.txt"
    run:
             if "AD" in wildcards.GWAS:
             
                 shell("""

                 module load htslib/1.9                  # For tabix
                 module load snpeff                      # For SnpSift    
                 
                 if [ -f {params.snps}.tbi ]; then
                   echo "{params.snps}.tbi exists."
                 else 
                   echo "Creating 00-All.vcf.gz.tbi"
                   tabix -p vcf {params.snps}
                 fi
                 
                 zcat {input} |\
                 cut -f1-4 |\
                 tail -n+2 |\
                 awk '{{print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\t."}}' |\
                 sort -k1,1V -k2,2g |\
                 java -jar $SNPSIFT annotate -id {params.snps} > {params.out_dir}AD_hg19_sumstats_fixHeader.tsv 2> {log}
                 
                 gzip {output} 2> {log}
                  
                 """)

             elif "NRTCSM" in wildcards.GWAS:

                 shell("""
 
                 zcat {input} | sed 's/POS/BP/g' | sed 's/REF/A1/g' | sed 's/ALT/A2/g' > {params.out_dir}NRTCSM_hg19_sumstats_fixHeader.tsv
                 gzip {output} 2> {log}

                 """)

             elif "BV" in wildcards.GWAS:

                 shell("""
 
                 zcat {input} | sed 's/SNP/UNIQUE_ID/g' | sed 's/RSID/SNP/g' > {params.out_dir}BV_hg19_sumstats_fixHeader.tsv
                 gzip {output} 2> {log}

                 """)

             elif "INTEL" in wildcards.GWAS:

                 shell("""
 
                 zcat {input} | sed 's/Zscore/Z/g' | sed 's/N_analyzed/N/g' > {params.out_dir}INTEL_hg19_sumstats_fixHeader.tsv
                 gzip {output} 2> {log}

                 """)

             else:

                 shell("cp {input} {output}")

rule standardise_sumstats2:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    input:   "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_sumstats_fixHeader.tsv.gz"
    output:  "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.tsv"
    message: "Formatting {input}"
    log:     "../results/logs/munge_sumstats/{GWAS}_standardise_sumstats2.log"
    shell:
             """ 
       
             python ../resources/python_convert/sumstats.py csv \
             --sumstats {input} \
       	     --out {output} --force --auto --head 5 \
             --log {log}

             """

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking (SCZ and BPD)
    input:   "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.tsv"
    output:  "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_withZ_sumstats.tsv"
    message: "Adding Z score to {input}"
    log:     "../results/logs/munge_sumstats/{GWAS}_addZscore.log"
    run:
             if wildcards.GWAS in ("INSMNA", "NRTCSM"):

                 shell("""

                 python ../resources/python_convert/sumstats.py zscore \
                 --sumstats {input} \
                 --out {output} --force \
                 --log {log} \
                 --a1-inc

                 """)

             else: 

                 shell("cp {input} {output}")

rule final_standardisation_edits:
    # SCZ - removing SNP with 13 cols, adding N, sort cols to match rest of sumstats files
    # BPD - Adding N
    # Format all GWAS to SNP, CHR, BP, PVAL, A1, A2, N, Z, OR, BETA/INFO SE
    input:   "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_withZ_sumstats.tsv"
    output:  "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
    message: "Adding Z score to {input}"
    log:     "../results/logs/munge_sumstats/{GWAS}_SCZ_edits.log"
    run:
             if "SCZ" in wildcards.GWAS:
             
                 shell("""
 
                 sed -i '/rs148878475/d' {input}
                 awk '{{s=(NR==1)?"N":"161405";$0=$0 OFS s}}1' {input} |\
                 awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$12"\t"$11"\t"$7"\t"$9"\t"$8}}' >\
                 {output} 2> {log} 
                 
                 """)

             elif "BPD" in wildcards.GWAS:

                 shell("""
 
                 awk '{{s=(NR==1)?"N":"413466";$0=$0 OFS s}}1' {input} |\
                 awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$12"\t"$11"\t"$7"\t"$9"\t"$8}}' >\
                 {output}
                 
                 """)
 
             elif "NRTCSM" in wildcards.GWAS:

                 shell("""
 
                 awk '{{s=(NR==1)?"N":"313467";$0=$0 OFS s}}1' {input} > {output}
                 
                 """)
             
             else:

                 shell("cp {input} {output}")


#rule sumstats_to_bed:
#     # Convert GWAS sumstats to bed format to run liftover to hg38
#    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
#    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.bed"
#    message: "Creating {input} bed file"
#    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_create_bed.log"
#    shell:
#             """
#             awk '{{print "chr"$2"\t"$3"\t"($3 + 1)"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' {input} |\
#             sed -e '1s/^/#/' > {output} 2> {log}
#             """


#rule lift_over:
#    # Lift GWAS data from hg19 to hg38
#    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.bed"
#    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.bed"
#    message: "Lifting {input} to hg38"
#    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_liftOver.log"
#    params:  chain = GWAS_DIR + "hg19ToHg38.over.chain.gz",
#             unlift = GWAS_DIR + "{wildcards.GWAS}_hg38_unlifted.bed"
#    shell:
#             """
#             ./liftover/liftOver {input} {params.chain} {output} {params.unlift} -bedPlus=4 2> {log}
#             """


#rule bed_to_sumstats:
#    # Convert hg38 bed file to sumstats format
#    input:   bed = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.bed",
#             header = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
#    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.tsv"
#    message: "Revert {input} to sumstats"
#    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_revert2sumstats.log"
#    shell:
#             """
             
#             cut --complement -d$'\\t' -f3 {input.bed} |\
#             sed -e 's/^chr//g' |\
#             awk '{{print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' |\
#             sed -e "1s/^/$(head -n1 {input.header})\\n/" > {output} 2> {log} 
#             """

#rule sumstats_for_magma:
#    # MAGMA needs PVAL to be labeled P
#    input:   hg19 = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv",
#             hg38 = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.tsv"
#    output:  hg19 = MAGMA_GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv",
#             hg38 = MAGMA_GWAS_DIR + "{GWAS}_hg38_magma_ready.sumstats.tsv"
#    message: "Munging sumstats for {input} for magma compatibility"
#    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_sumstats_for_magma.log"
#    shell:   
#             """
#             sed 's/PVAL/P/g' {input.hg19} > {output.hg19}
#             sed 's/PVAL/P/g' {input.hg38} > {output.hg38}
#             """

rule sumstats_for_ldsc_hg19:
    input:   snps = "../resources/ldsc/reference_files/w_hm3.snplist",
             gwas = "../results/gwas_sumstats/GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
    output:  "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz"
    conda:   "../envs/ldsc.yml"
    message: "Munging sumstats for {input.gwas} for ldsc compatibility"
    params:  out = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready"
    log:     "../results/logs/munge_sumstats/{GWAS}_hg19_sumstats_for_ldsc.log"
    shell:
        """

        python ../resources/ldsc/munge_sumstats.py --sumstats {input.gwas} \
        --merge-alleles {input.snps} \
        --out {params.out} \
        --a1-inc 2> {log}

        """

#rule sumstats_for_ldsc_hg38:
#    input:   snps = SCRATCH + "ldsc/reference_files/w_hm3.snplist",
#             gwas = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_withZv2_sumstats.tsv"
#    output:  LDSC_GWAS_DIR + "{GWAS}_hg38_ldsc_ready.sumstats.gz"
#    conda:   "envs/ldsc.yml"
#    message: "Munging sumstats for {input.gwas} for ldsc compatibility"
#    params:  out = LDSC_GWAS_DIR + "{GWAS}_hg38_ldsc_ready"
#    log:     "logs/munge_sumstats/{GWAS}_hg38_sumstats_for_ldsc.log"
#    shell:
#        """
#        python ldsc/munge_sumstats.py --sumstats {input.gwas} \
#        --merge-alleles {input.snps} \
#        --out {params.out} \
#        --a1-inc 2> {log}
#        """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
