#!/usr/bin/env bash
set -euo pipefail

snakemake --use-conda --configfile config/config.yaml -p $@

