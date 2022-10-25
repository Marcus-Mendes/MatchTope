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

# How to cite:
How to cite:

Matchtope: (https://doi.org/10.3389/fimmu.2022.930590).


PIPSA: https://projects.h-its.org/mcmsoft/pipsa/4.0.2/references.html

Blomberg N, Gabdoulline RR, Nilges M, and Wade RC.
Classification of protein sequences by homology modeling and quantitative analysis of electrostatic similarity.
Proteins: Str., Function and Genetics 1999, 37: 379-387

Wade RC, Gabdoulline RR and De Rienzo F.
Protein Interaction Property Similarity Analysis.
Intl. J. Quant. Chem. 2001, 83: 122-127.

UHBD:

Madura, Jeffry D., et al.

Electrostatics and diffusion of molecules in solution: simulations with the University of Houston Brownian Dynamics program. Computer Physics Communications 1995, 91 (1-3): 57-95.

