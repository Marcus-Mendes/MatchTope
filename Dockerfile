# MatchTope Dockerfile
## Author: BragatteMAS
## Last Update: 2023-11

FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Installing dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget -qO Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/miniconda/bin:${PATH}"

# Creating a Conda environment with Python 3.7 and installing PyMOL in one RUN command
RUN conda create -n pymol_env python=3.7 -y \
    && echo "source activate pymol_env" > ~/.bashrc \
    && /bin/bash -c "source activate pymol_env && conda install -c schrodinger pymol -y"

# Rest of the installations
RUN apt-get update && apt-get install -y r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('pvclust', repos='http://cran.rstudio.com/')"

COPY . /MatchTope
WORKDIR /MatchTope

CMD ["bash", "run_pipsa.sh"]
