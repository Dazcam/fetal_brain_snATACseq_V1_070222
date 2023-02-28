# -------------------------------------------------------------------------------------
#
#    Remove MHC file from peaks (hg19)
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule remove_MHC_from_peaks:
    input:   bed = "../results/peaks/{CELL_TYPE}.hg19.bed",
             MHC = "../results/peaks/MHC.hg19.bed"
    output:  "../results/peaks/{CELL_TYPE}.hg19.noMHC.bed"
    message: "Removing MHC regions from {wildcards.CELL_TYPE}"
    log:     "../results/logs/remove_MHC/{CELL_TYPE}.hg19.noMHC.log"
    shell:
        """
        
        module load bedtools
	bedtools subtract -a {input.bed} -b {input.MHC} -A > {output}

        """ 
