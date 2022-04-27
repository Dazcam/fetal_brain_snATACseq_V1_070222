#--------------------------------------------------------------------------------------
#
#   Munge Garfield UK10K SNPs in separte folders into Bedtools readable
#
#--------------------------------------------------------------------------------------

# Info  -------------------------------------------------------------------------------

# Garfiled SNP lists are contained in separate space delimted chr folders 
# Bedtools does not like this so need to create bed files for each SNP list for 
# Intersection with annotation file

# Set variables  ----------------------------------------------------------------------
OUTDIR=$1

echo -e '\nThe outdir is' $OUTDIR

# Format SNP list  --------------------------------------------------------------------
echo -e '\nFormatting SNP list ... \n'
for file in ../resources/garfield/garfield-data/maftssd/chr*
do

  VAR=$(basename $file)
  NUM=$(echo $file | cut -d'r' -f 2)

  awk '{print $1 "\t" $1+1}' ${file} > col1.2.${NUM} # Extract col 1 / +1 for col 2 
  awk -v var="$VAR" '{ print var"\t"$0}' col1.2.${NUM} > ${OUTDIR}${VAR} # Add chr number to col 1

  # Check if file lengths match  -------------------------------------------------------
  # If so delete intermediate files
    
  if [ "$(wc -l < ${file})" -eq "$(wc -l < ${OUTDIR}${VAR})" ]

      then
          echo ${VAR}' Match!'
          rm col1.2.${NUM}
      else echo 'Warning: '${VAR}' No Match!'
  fi

done

# Quick visual check of text formatting of bed files
head ${OUTDIR}* | sed -n 'l'

echo 'Done.'
# If reformatting is required
# Change delimiter 
# for file in *.txt; do cat $file | tr '[,]' '[\t]' > new_*.txt; done
# See - https://www.biostars.org/p/330003/

# When the error message "has non positional records, which are only valid for the groupBy too" appears,do
# sed -n 'l' <filename>
# to check invisible characters
# If line ending is r$, to remove it do
# sed 's/\r//' <filename>



#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
