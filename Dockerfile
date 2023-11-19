# MatchTope Dockerfile
## Author: BragatteMAS
## Last Update: 2023-11

# Base image with Python suitable for installing R and PyMOL
FROM ubuntu:22.04

## Avoid interactive prompts during build
ARG DEBIAN_FRONTEND=noninteractive

## Update and install necessary dependencies, including Python
RUN apt-get update && apt-get install -y \
    python3.7 \
    python3-pip \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

## Install specific version of R (replace '3.6.3' with the desired version)
RUN apt-get update && apt-get install -y software-properties-common \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 \
    && apt-get install -y r-base=3.6.3-1focal \
    && rm -rf /var/lib/apt/lists/*

## Install specific R packages
RUN R -e "install.packages('pvclust', repos='http://cran.rstudio.com/')"

## Install PyMOL open-source (replace '2.5.2' with the desired version)
RUN wget -qO- https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh | bash \
    && /root/miniconda3/bin/conda install -c schrodinger pymol=2.5.2

## Copy the repository files into the container
COPY . /MatchTope

## Set the working directory
WORKDIR /MatchTope

## Uncomment the line below to install Python dependencies from a requirements.txt file
# RUN python3.7 -m pip install -r requirements.txt

## Execute the main script when starting the container
CMD ["bash", "run_pipsa.sh"]
