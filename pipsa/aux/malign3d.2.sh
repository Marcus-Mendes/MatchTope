#! /bin/sh

# This program aligns the sequences to the TEMPLATE.pdb sequence, using modeller and modeller.2.in script.


cd ../pdbs
ls -1 *.pdb > names

# mkdir ../malign3d
cd ../malign3d/
mkdir pdbs
rm -f u.pdb modeller.log 


for I in `cat ../pdbs/names`
do
    echo "Now aligning $I ..."
    ln -s ../pdbs/$I u.pdb
    mod ../source/modeller.2.in
    if [ -f tmp.pdb ]; then
	mv tmp.pdb pdbs/$I
    else
	echo "Aligning $I was unsuccessful... "
    fi
    echo "---------\n$I\n\n" >> modeller.log
    cat ../source/modeller.2.in.log >> modeller.log
    rm -f u.pdb
done
echo "Done."
