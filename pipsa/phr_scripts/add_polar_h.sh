#!/bin/bash
for pdbfile in *.pdb ; do 
              echo -e "13\r\n1\r" | pdb2gmx -f $pdbfile -o $pdbfile-withH.pdb
              rm *.itp *.top
              awk -v values=". HOH" 'BEGIN{split(values,array)};{flag=0; for(val in array) if (array[val] == $4) flag=1; if (flag==0) print}'  $pdbfile-withH.pdb > Prepared-$pdbfile
              rm *-withH.pdb
done
