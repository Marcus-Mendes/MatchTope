#!/bin/bash

# Navigate to pipsa/exa/PHR/pdbs and remove all PDB files
cd pipsa/exa/PHR/pdbs
rm -f *.pdb

# Copy files from PDBs to the current directory
cp -r /MatchTope/PDBs/. .

# Run fitting.py in the current directory
python fitting.py

# Further operations
printf '%s\n' *.pdb > pdbnames
cd ../
./do_PHR_com
cd uhbd/
../../aux/pipsa2R.pl -s sims.log -t h -m matrix -g
mv matrix_2_h.pdf Results.pdf
mv Results.pdf /MatchTope/
