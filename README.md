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

## Running MatchTope with Docker

MatchTope can be easily run in a Docker container, which encapsulates all its dependencies and provides a consistent running environment. This section guides you through the process of using Docker to run MatchTope.

### Prerequisites:

- Docker Desktop (for Windows/Mac) or Docker Engine (for Linux)
- Git (for cloning the repository)

### Steps:

1. **Clone the Repository:**

First, clone the MatchTope repository to your local machine:

```bash
git clone https://github.com/Marcus-Mendes/MatchTope.git
cd MatchTope
```

2. **Build the Docker Image:**

In the root directory of the project, build the Docker image:

```bash
docker build -t matchtope .
```

This command reads the Dockerfile in the current directory and builds an image named "matchtope".

2.1 **Windons OS**

Install Docker Desktop on Windows: If you haven't already, download and install Docker Desktop for Windows from the official Docker website.

* Enable WSL 2 Integration in Docker Desktop:
  * Open Docker Desktop.
  * Go to Settings (the gear icon).
  * Click on Resources > WSL Integration.
  * Enable integration with your WSL 2 distro (e.g., Ubuntu 22.04).
  
  * Apply & Restart Docker Desktop.
        
    Restart WSL 2: Sometimes, a restart of the WSL 2 instance is required for changes to take effect. 

    You can restart WSL 2 by running the following command in your Windows Command Prompt or PowerShell:

    ```bash
    wsl --shutdown
    ```
    
    ![IMAGE DEMONSTRATION.IMG]()

  Then, reopen your WSL 2 Ubuntu terminal.
  


* Verify Docker Installation:

After the restart, in your WSL 2 Ubuntu terminal, check if Docker is accessible by running:

```bash
docker --version
```

You should see the Docker version if it's installed correctly.

2.2 **Build Your Docker Image om Windows WSL:**

Now, try building your Docker image again in your WSL 2 environment:

```bash
docker build -t matchtope .
```

By enabling WSL 2 integration in Docker Desktop and ensuring Docker commands are accessible within your WSL 2 Ubuntu distro, you should be able to build Docker images directly from WSL 2.

3. **Prepare Input Data:**

Place your PDB files into the PDBs folder. If this folder does not exist, create it in the root directory of the project.

4. **Run the Docker Container:**

Execute the following command to run MatchTope in a Docker container:

```bash
docker run -v "${PWD}/PDBs:/MatchTope/PDBs" -v "${PWD}/Results:/MatchTope/Results" matchtope
```

This command mounts the PDBs folder from your local machine to the container and sets up a Results folder for the output.

5. **Accessing Results:**

After the container finishes running, the output (e.g., Results.pdf) will be available in the Results folder on your local machine.

**Note for Linux Users:**

Replace `${PWD}` with `$PWD` in the Docker run command to reference the current directory.


# How to cite:

### [Matchtope article](https://doi.org/10.3389/fimmu.2022.930590)
> M. F. de A. Mendes et al., “MatchTope: A tool to predict the cross reactivity of peptides complexed with Major Histocompatibility Complex I,” Front. Immunol., vol. 13, 2022, Accessed: Jan. 26, 2023. [Online]. Available: https://www.frontiersin.org/articles/10.3389/fimmu.2022.930590


### [PIPSA article](https://projects.h-its.org/mcmsoft/pipsa/4.0.2/references.html)
> Blomberg N, Gabdoulline RR, Nilges M, and Wade RC. Classification of protein sequences by homology modeling and quantitative analysis of electrostatic similarity. Proteins: Str., Function and Genetics 1999, 37: 379-387

> Wade RC, Gabdoulline RR and De Rienzo F. Protein Interaction Property Similarity Analysis. Intl. J. Quant. Chem. 2001, 83: 122-127.

### UHBD:
> Madura, Jeffry D., et al. Electrostatics and diffusion of molecules in solution: simulations with the University of Houston Brownian Dynamics program. Computer Physics Communications 1995, 91 (1-3): 57-95.