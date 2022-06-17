#!/bin/bash
for pdbfile in *.pdb ; do
	sed -ri 's/vina/    /g'  $pdbfile
done
x=1
while [ $x -le 60 ]
do
	for pdbfile in *.pdb ; do
		sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
x=$(( $x + 1 ))
	done
done

