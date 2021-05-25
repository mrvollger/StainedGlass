#!/usr/bin/env bash
set -euo pipefail

snakemake --use-conda --configfile StainedGlass.yaml -p $@

