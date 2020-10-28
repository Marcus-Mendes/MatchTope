#!/bin/bash
x=1
while [ $x -le 5 ]
do
	for pdbfile in *.pdb ; do
		sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
x=$(( $x + 1 ))
	done
done

