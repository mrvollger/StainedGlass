#!/usr/bin/env bash
set -euo pipefail
smk=/net/eichler/vol26/projects/sda_assemblies/nobackups/software/miniconda3/envs/higlass/bin/snakemake
which snakemake 
which $smk

if [ $1 == "local" ]; then
  $smk -s ./dot_aln.smk \
    -p -j 250 -k \
    "${@:2}"
else 
  snakemake -s ./dot_aln.smk \
    --drmaa-log-dir logs \
    --drmaa " -w n -P eichlerlab -q eichler-short.q -l h_rt=48:00:00  -l mfree={resources.mem}G -pe serial {threads} -V -cwd -S /bin/bash" \
    -p -j 500 -k \
    "${@:2}"

fi

    #--use-conda \
    #--use-conda --conda-frontend mamba \
