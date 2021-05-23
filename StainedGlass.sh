#!/usr/bin/env bash
set -euo pipefail

snakemake -p --use-conda --configfile StainedGlass.yaml $@

