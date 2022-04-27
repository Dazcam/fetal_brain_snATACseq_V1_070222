# ------------------------------------------------------------------------------
#
#   Garfield - Create custom annotation file using UK10K genome sequence data
#
#   Intersect formatted UK10K SNPs with annotation file
#
# ------------------------------------------------------------------------------

# Set variables  ---------------------------------------------------------------
BED_FILE=$1
OUTDIR=$2

#------------- Intersect -------------------------
for peakFile in ${BED_FILE}; do

  PREFIX=$(basename "$peakFile" .bed)
  echo -e '\nIntersecting' ${peakFile}

  for file in $(ls ../results/garfield/uk10_SNPs/*); do
    
    # Intersect 
    echo ${file}
    NUM=$(basename $file)
    bedtools intersect -a ${file} -b ${peakFile} -c > ${OUTDIR}${NUM}
     
    # Retain columns 2 and 4 (Format required for Garfield)
    cut -f2,4 ${OUTDIR}${NUM} > ${OUTDIR}${NUM}.temp
    mv ${OUTDIR}${NUM}.temp ${OUTDIR}${NUM}
    
    # Create link file
    echo "Index Annotation Celltype Tissue Type Category" > ${OUTDIR}link_file.txt
    echo "0" ${OUTDIR} "NA NA NA NA" >> ${OUTDIR}link_file.txt

  done

done


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
