#! /bin/sh

# To use this script, malign3d directory should be present
# also TEMPLATE.pdb file should be present in this directory 
# The program aligns the proteins to the target sequence. 

cd ../pdbs/
mkdir ../malign3d/pdbs/ # used to store new pdb files
rm -f ../malign3d/whatif.out # remove old log file
rm -f ../malign3d/out.pdb # remove old src file.
for I in `ls *.pdb`
do

    cp $I ../malign3d/in.pdb
    cd ../malign3d/
    echo "Running whatif on file $I ..."
    echo "-------------------\n FILE: $I \n" >>whatif.out
    whatif <../source/malign3d.whatif.in >>whatif.out 2>&1
    if [ ! -f out.pdb ]; then
	echo "WARNING: Calculations for $I were unsuccessful."
    else
	mv out.pdb "pdbs/$I"
    fi
    rm -f in.pdb
    cd ../pdbs/
done

rm in.pdb

echo "Finished."