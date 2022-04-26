#!/bin/bash

BED_FILE=$1
OUTPUT_DIR=$2
VARIANTS_DIR=../resources/garfield-data/maftssd 

## loop over chromosomes
for chr in {1..22}
do
	
	INFO=$VARIANTS_DIR/chr${chr}
	echo "CHROM POS" > $OUTPUT_DIR/chr${chr}
	awk -v chr=$chr '{print chr,$1}' $INFO > $OUTPUT_DIR/chr${chr}

	## loop over annotation files
	for f in $BED_FILES
	do
	
		## number of columns in the .bed file
		NCOL=$(head -1 $f | awk '{print NF}')
		## make a local copy removing possible headers
		grep "chr" $f |awk -v chr=$chr '$1=="chr"chr {print}' | sort -k2n  > $OUTPUT_DIR/tmp.chr$chr.bed

	        echo "Processing $f file"
		## running annotation part
	        ../resources/garfield/garfield-v2/garfield_annotate_uk10k_helper --ncol $NCOL --o $OUTPUT_DIR/tmp_chr${chr} --peaks $OUTPUT_DIR/tmp.chr$chr.bed --norsid --chunk 1000 --info $INFO
		## merge
		cp $OUTPUT_DIR/chr${chr} $OUTPUT_DIR/tmp0_chr${chr}
        	paste $OUTPUT_DIR/tmp0_chr${chr} <(awk '{print $4}' $OUTPUT_DIR/tmp_chr${chr}) > $OUTPUT_DIR/chr${chr}

		## clean up
		rm $OUTPUT_DIR/tmp0_chr${chr}
		rm $OUTPUT_DIR/tmp_chr${chr}
		rm $OUTPUT_DIR/tmp.chr$chr.bed
	done
done

## create link_file.txt (needed for running GARFIELD)
LINK_FILE=$OUTPUT_DIR/link_file.txt
echo "Index Annotation Celltype Tissue Type Category" > $LINK_FILE
i=-1
for f in $BED_FILES;do i=$[$i+1]; echo $i $f "NA" "NA" "NA" "NA"  >> $LINK_FILE ; done
fi

## reformat data for GARFIELD usage
f=$OUTPUT_DIR/chr${chr}
cat $f > $f.tmp
paste -d" " <(awk '{print $2}' $f.tmp | sed 1d) <(awk '{$1=$2=""; print $0}' $f.tmp | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > $f
rm $f.tmp

