# MatchTope
MatchTope is a tool developed for predicting peptide similarity, which can trigger cross-reactivity events by computing and analyzing the electrostatic potentials of pMHC complexes. MatchTope uses a modified version of the PIPSA tool. The PIPSA standalone can be obtained at http://pipsa.h-its.org.

# Requirements:
- Python 3.7 version or higher
- R version 3.6.3 or higher
- pvclust (R package) version 2.2-0 or higher (apt install r-cran-pvclust)
- Pymol version 2.5.2 or higher (It can be installed with conda: conda install -c conda-forge pymol-open-source)
- Modeled pMHC or crystal files

# How to run:
- Put your pdb files into the PDBs folder. It already has three files that can be used to test the tool.
- In your terminal, run 'bash run_pipsa.sh'.
- When it finishes, it will generate a pdf file called Results.pdf in the same folder of run_pipsa.sh.
- There is a file called 'Results_example.pdf'. It is an output from the three PDBs examples in the PDBs folder.


