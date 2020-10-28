#!/bin/bash
echo -e "Where is the SW folder?"; read caminho
cd $caminho/sw/pipsa/exa/PHR/pdbs/
./add_polar_h.sh
python fitting.py
python fitting.py
python fitting.py
#./remove_characters_pdb.sh
printf '%s\n' *.pdb > pdbnames
cd ../
./do_PHR_com
cd uhbd/
./../../../aux/pipsa2R.pl -s sims.log -t h -m matrix -g

