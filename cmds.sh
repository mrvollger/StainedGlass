#!/usr/bin/env bash
set -euo pipefail
function useage {
    echo "useage: ./cmds.sh prefix fasta windowsize "
    exit 1
}

if [ $# -lt 3 ]; then
    useage
    exit 1
fi

ODIR="results"
PRE="$ODIR/$1"
FASTA=$2
W=$3
mkdir -p $ODIR


if [ ! -f $PRE.tbl ]; then
  bedtools makewindows -g $FASTA.fai -w $W > $PRE.bed
  bedtools getfasta -fi $FASTA -bed $PRE.bed > $PRE.fasta

  minimap2 -t 80 -f 0.0001 --eqx -ax ava-ont $PRE.fasta $PRE.fasta \
    | samtools view -b -F 4 - \
    | samtools sort -@ 6 -m 20G - \
    > $PRE.bam

  ./samIdentity.py --header $PRE.bam > $PRE.tbl
fi

./refmt.py --window $W --fai $FASTA.fai $PRE.tbl > $PRE.aln.bed

