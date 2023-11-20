#!/bin/bash

cd /pipsa/exa/PHR/pdbs
rm *.pdb
cd ../../../../
cp -r /PDBs/. pipsa/exa/PHR/pdbs
cd /pipsa/exa/PHR/pdbs
# ./add_polar_h.sh
python fitting.py
# ./remove_characters_pdb.sh
printf '%s\n' *.pdb > pdbnames
cd ../
./do_PHR_com
cd uhbd/
./../../../aux/pipsa2R.pl -s sims.log -t h -m matrix -g
mv matrix_2_h.pdf Results.pdf
mv Results.pdf ../../../../
