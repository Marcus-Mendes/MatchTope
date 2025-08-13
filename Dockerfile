# MatchTope Dockerfile
## Author: BragatteMAS
## Edits: Norgent
## Last Update: 2025-08

FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Installing general dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    dos2unix \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget -qO Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/miniconda/bin:${PATH}"

# Accept channel TOS and add conda-forge channe;
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN /bin/bash -c "conda config --prepend channels conda-forge"

# Creating a Conda environment with Python 3.10 and installing PyMOL
RUN conda create -n matchtope_env python==3.10 pymol-open-source -c conda-forge \
    && echo "source activate matchtope_env" > ~/.bashrc

# Activate conda environment 
RUN /bin/bash -c "source activate matchtope_env"

# Install R and R packages
RUN apt-get update && apt-get install -y r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/* \
    && R -e "install.packages('pvclust', repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('gplots', dependencies=TRUE, repos='http://cran.rstudio.com/')"


# Copying the project files into the container and ensuring Unix-style line endings
COPY . /MatchTope

WORKDIR /MatchTope

CMD ["/bin/bash", "-c", "source activate matchtope_env && bash run_pipsa.sh"]
