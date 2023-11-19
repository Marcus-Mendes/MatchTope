# MatchTope Dockerfile
## Author: BragatteMAS
## Last Update: 2023-11

# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

## Avoid interactive prompts during build
ARG DEBIAN_FRONTEND=noninteractive

## Update and install necessary dependencies, including Python and R
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    r-base-dev \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

## Install specific R packages (if needed)
RUN R -e "install.packages('pvclust', repos='http://cran.rstudio.com/')"

## Install PyMOL open-source
RUN wget -qO- https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh | bash \
    && /root/miniconda3/bin/conda install -c schrodinger pymol

## Copy the repository files into the container
COPY . /MatchTope

## Set the working directory
WORKDIR /MatchTope

## Uncomment the line below to install Python dependencies from a requirements.txt file
# RUN python3 -m pip install -r requirements.txt

## Execute the main script when starting the container
CMD ["bash", "run_pipsa.sh"]
