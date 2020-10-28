#!/bin/sh

# I am using modeller's command orient to alignt the sequences together... just quick hack to make it work...
# This command translates the MODEL so that its gravity center is at the 
# origin of the coordinate system and that the three principle axes of the 
# model's inertia ellipsoid correspond to the x, y, and z axes of the 
# coordinate system. It may even be used for approximate superposition if 
# molecules have a similar non-spherical shape. Information about the 
# principal axes is written to the log file. 

mkdir ../pdbs_orient

cd ../pdbs

for I in `ls *.pdb`
do
    echo "READ_MODEL FILE = '$I', MODEL_FORMAT = 'PDB'
ORIENT_MODEL
WRITE_MODEL FILE = '../pdbs_orient/$I'" > modeller.top

mod ../pdbs/modeller.top




done
