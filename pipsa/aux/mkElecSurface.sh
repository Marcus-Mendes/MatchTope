#!/bin/sh 

# The program will make surface with -2.5 & 2.5 from a protein .pdb file and electrostatic potential grid made by uhbd. 
# The potential grid should be present in ../uhdb folder and should have the name of the protein file with extension .grd. The pdb files should be located in ../pdbs folder. 
# The program uses insightII and grid2i65 from the uhbd package.

mkdir ../elec/

cd ../elec/

# Remove old stuff...
rm -f in.pdb in.ins outfile.wrl grid2i65.log insightII.log outfile.wrl

# Make dir for the converted insight files and convert them. 
#rm -f ins/* 

echo "Converting the uhbd file format into insightII file format..."
mkdir ins
for I in `cut -f1 -d. ../pdbs/names`
do
    echo "Now converting file $I.grd..."
    echo "B\n../uhbd/$I.grd\nins/$I.ins" | ../source/grid2i65 >>grid2i65.log # Convert the stuff
done


# Create the .wrl files that can be viewed with netscape.
rm -f wrl/*
mkdir wrl

echo "Running insightII..."
for I in `cut -f1 -d. ../pdbs/names`
do
    ln -s ../pdbs/$I.pdb in.pdb
    ln -s ins/$I.ins in.ins
    echo "Now running insightII on file $I..."
# A bug... insightII produces different output with and without -at (batch mode) switch on.
    insightII < ../source/insightII.electSurface.in >>insightII.log # Make the files
    mv outfile.wrl wrl/$I.wrl
    rm -f in.pdb in.ins
done

