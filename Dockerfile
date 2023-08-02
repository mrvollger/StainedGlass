FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fd0eeef546ed00fa95ff85a8baf143dfccd52a6183b57e9bb8e972bf196948d5"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/ee4bfa88159ef77575445596ec2221ce
#   name: R
#   channels:
#     - conda-forge
#     - r
#     - defaults
#   dependencies:
#     #- xorg-libx11
#     #- xorg-libxau
#     #- r::r-ggplot2
#     - r-base>=4.0
#     - r-essentials
#     - r-cairo
#     - r-data.table
#     - r-cowplot
#     - r-argparse>=2.1.2
#     - r-glue
#     - r-r.utils
#     - r-rcolorbrewer
#     - r-scales
#     - r-tidyverse>=1.3.0
RUN mkdir -p /conda-envs/ee4bfa88159ef77575445596ec2221ce
COPY workflow/envs/R.yaml /conda-envs/ee4bfa88159ef77575445596ec2221ce/environment.yaml

# Conda environment:
#   source: workflow/envs/env.yaml
#   prefix: /conda-envs/5ea0207c595d8a051ab13a13766b14f6
#   name: python_and_cli_env
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - numpy<1.20
#     - numba
#     - cooler<=0.8.11
#     - pandas<=1.5
#     - minimap2==2.18
#     - bioconda::bedtools
#     - bioconda::samtools>=1.14
#     - bioconda::htslib>=1.14
#     - bioconda::pysam>=0.15.0
#     - bioconda::bwa
#     - pigz
RUN mkdir -p /conda-envs/5ea0207c595d8a051ab13a13766b14f6
COPY workflow/envs/env.yaml /conda-envs/5ea0207c595d8a051ab13a13766b14f6/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/ee4bfa88159ef77575445596ec2221ce --file /conda-envs/ee4bfa88159ef77575445596ec2221ce/environment.yaml && \
    mamba env create --prefix /conda-envs/5ea0207c595d8a051ab13a13766b14f6 --file /conda-envs/5ea0207c595d8a051ab13a13766b14f6/environment.yaml && \
    mamba clean --all -y
