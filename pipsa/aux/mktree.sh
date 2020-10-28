#! /bin/sh -xv

# The program should create an unrooted tree which represents the relationships between the proteins. 
# For the program to work you should have neighbour / drawtree / fontfile / ps2pdf in source directory or in the path (in case of ps2pdf)

# As an argument either ../uhbd or ../grid should be given

DIR=$@
if [ -z DIR ]; then
echo "You should give a name of the directory with the sims.mat file"
exit 1
fi

# remove the old files
rm -f $DIR/infile $DIR/fontfile $DIR/plotfile $DIR/plotfile.ps ../pdbs/treefile ../pdbs/outfile ../pdbs/plotfile.pdf

# Do the stuff
echo "Converting the distance matrix into the Phylip file format..."
mat2phyl.pl < $DIR/sims.mat > $DIR/infile # convert the distance matrix into something digestable by phylip
cd $DIR
# ln -s ../source/fontfile fontfile
echo "Preparing treefile..."
../source/neighbour < ../source/neighbour.in # Creates 
echo "Drawing the tree..."
# sleep 3
cp ../source/fontfile fontfile # required by drawtree
# ../source/drawtree <../source/drawtree.in # Creates plotfile
../source/drawgram <../source/drawgram.in # Creates plotfile
echo "Changing the links..."
../source/addlinks.pl plotfile 'http://www.expasy.ch/cgi-bin/sprot-search-de?SEARCH=' > plotfile.ps
echo "Converting the file to .pdf file..."
ps2pdf plotfile.ps plotfile.pdf

echo "Done."
echo "The files plotfile.ps and plotfile.pdf contain the tree."
