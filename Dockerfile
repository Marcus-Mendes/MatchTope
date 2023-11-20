# MatchTope Dockerfile
## Author: BragatteMAS
## Last Update: 2023-11

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

# Creating a Conda environment with Python 3.7 and installing PyMOL
RUN conda create -n pymol_env python=3.7 -y \
    && echo "source activate pymol_env" > ~/.bashrc

# Activate conda environment and install PyMOL
RUN /bin/bash -c "source activate pymol_env && conda install -c schrodinger pymol -y"

# Install R and R packages
RUN apt-get update && apt-get install -y r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/* \
    && R -e "install.packages('pvclust', repos='http://cran.rstudio.com/')"

# Copying the project files into the container and ensuring Unix-style line endings
COPY . /MatchTope
RUN dos2unix /MatchTope/run_pipsa.sh \
    && dos2unix /MatchTope/do_PHR_com \
    && dos2unix /MatchTope/pipsa2R.pl \
    && chmod +x /MatchTope/run_pipsa.sh \
    && chmod +x /MatchTope/do_PHR_com \
    && chmod +x /MatchTope/pipsa2R.pl

WORKDIR /MatchTope

CMD ["/bin/bash", "-c", "source activate pymol_env && bash run_pipsa.sh"]