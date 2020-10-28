#!/bin/sh

# The program obtains the proper names from the Swiss-Prot database for the proteins, which file names in ../pdbs/pdbnames are based on accession numbers.
# The result is saved in propernames file.

cd ../pdbs/

rm propernames

for I in `cat pdbnames | cut -f1 -d.`
do
    LINK="http://www.expasy.ch/cgi-bin/sprot-search-de?S=1&T=1&SEARCH=$I"
    NAME=`wget -O - --quiet $LINK | perl -we '$f=0; while(<>){$f=1 if (/Entry name/); if ($f==1){if(/<b>(.*?)<\/b>/){print "$1";$f=0;}}}'`
    echo $NAME >> propernames

done
