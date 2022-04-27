#!/bin/bash

# This script was modified from the garfield executable

INPUTNAME=Hannon_mQTL 
DATADIR=../resources/garfield/garfield-data

PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd
PVALDIR=../results/garfield/QTL_for_Garfield/Hannon_mQTL
ANNOTDIR=$1
OUTDIR=$2

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
PTHRESH=1e-5,1e-8
BINNING=m5,n5,t5
CONDITION=0
SUBSET="1-1005"

F1=$OUTDIR/garfield.prep.$INPUTNAME.out
F0=$OUTDIR/garfield.Meff.$INPUTNAME.out

echo 'Prune and Clump'
echo -n > $F1
for CHR in `seq 1 22` #X
do
	echo 'CHR'$CHR
	../resources/garfield/garfield-v2/garfield-prep-chr -ptags $PRUNETAGSDIR/chr$CHR -ctags $CLUMPTAGSDIR/chr$CHR -maftss $MAFTSSDDIR/chr$CHR -pval $PVALDIR/chr$CHR -ann $ANNOTDIR/chr$CHR -excl 895,975,976,977,978,979,980 -chr $CHR -o $F1 || { echo 'Failure!'; } 
done

echo 'Calculate effective number of annotations'
/apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript ../resources/garfield/garfield-v2/garfield-Meff-Padj.R -i $F1 -o $F0
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

echo 'Calculate Enrichment and Significance'
F2=$OUTDIR/garfield.test.$INPUTNAME.out
/apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript ../resources/garfield/garfield-v2/garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION
echo 'GARFIELD single annotation analysis complete'

echo 'Create Plots'
/apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript  ../resources/garfield/garfield-v2/garfield-plot.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

# echo 'Prioritize relevant annotations by conditional analysis'
# Additional prioritization of annotations
# CONDITION=1
# CONDITIONTHRESH=0.05
# Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION -ct $CONDITIONTHRESH -padj $Padj
#	echo 'GARFIELD model selection complete'


# Step to extract variants driving the enrichment signals
# echo 'Extracting variants driving enrichment analysis signals'
# PTHRESH=1e-8 ## GWAS threshold
# PENRICH=1e-9 ## Enrichment significance threshold - by default you might want to use $Padj
# GARFIELD_significant_annotations=$F2.significant.annotations.$PTHRESH.$PENRICH
# GARFIELD_VARS=${GARFIELD_significant_annotations}.variants
#./garfield_extract_variants_overlapping_enriched_annotations.sh $PRUNETAGSDIR $CLUMPTAGSDIR $ANNOTDIR $PVALDIR $F1 $F2 $GARFIELD_significant_annotations $GARFIELD_VARS $PTHRESH $PENRICH
# echo 'Extraction of Variants Complete'

echo 'GARFIELD Analysis Complete!'
