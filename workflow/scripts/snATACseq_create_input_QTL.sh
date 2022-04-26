## This script can be used to put GWAS summary statistics files into GARFIELD input format
## NOTE: example data is not supplied and script should be modified so that

## This script was modified for testing QTL enrichment - DC 260422


chrcol=1 ## the column in GWAS file containing chormosome information
poscol=2 ## the column in GWAS file containing genomic position information
pvalcol=3 ## the column in GWAS file containing GWAS p-value information
GWASFILENAME=$1 ## name of file containing GWAS summary statistics
OUTDIR=$2

for CHR in {1..22} #'X'
do
awk -v chr=$CHR -v chrcol=$chrcol -v poscol=$poscol -v pvalcol=$pvalcol '$chrcol==chr {print $poscol,$pvalcol}' $GWASFILENAME | sort -k1n > $OUTDIR/chr$CHR
echo $chr
done
