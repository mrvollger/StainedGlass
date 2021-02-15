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
F=10000 
#0.00001
PRE="$ODIR/$1"
FASTA=$2
W=$3
mkdir -p $ODIR

samtools faidx $FASTA

awk '{print $1 "\t" $2 }' $FASTA.fai \
            > $FASTA.chrom.sizes

alntype="cram"

if [ ! -f $PRE.$F.$alntype ]; then
  bedtools makewindows -g $FASTA.fai -w $W > $PRE.bed
  bedtools getfasta -fi $FASTA -bed $PRE.bed > $PRE.fasta

  minimap2 -t 80 -f $F --eqx -ax ava-ont $PRE.fasta $PRE.fasta \
    | samtools view -u -F 4 - \
    | samtools sort -m 100G --write-index \
    -O CRAM -T $FASTA \
    -o $PRE.$F.$alntype
fi

if [ ! -f $PRE.$F.tbl ]; then
  ./samIdentity.py --matches 1000 --header $PRE.$F.$alntype > $PRE.$F.tbl
fi

if [ ! -f $PRE.$F.aln.bed ]; then
  ./refmt.py --window $W --fai $FASTA.fai $PRE.$F.tbl > $PRE.$F.aln.bed
fi


if [ "x" == "y" ]; then
cat $PRE.$F.aln.bed | tail -n +2 |  cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 \
            --field count=7:agg=mean,dtype=float \
            --chunksize 50000000000 \
            $FASTA.fai:5000 \
	        --zero-based \
            - $PRE.$F.cool

cooler zoomify --field count:agg=mean,dtype=float $PRE.$F.cool
fi
