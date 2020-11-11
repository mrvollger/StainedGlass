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

awk '{print $1 "\t" $2 }' $FASTA.fai \
            > $FASTA.chrom.sizes
./refmt.py --window $W --fai $FASTA.fai $PRE.tbl > $PRE.aln.bed
cat $PRE.aln.bed | tail -n +2 |  cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 \
            --field count=7:agg=mean,dtype=float \
            --chunksize 50000000000 \
            $FASTA.fai:5000 \
	        --zero-based \
            - $PRE.cool
cooler zoomify --field count:agg=mean,dtype=float $PRE.cool

