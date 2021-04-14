# Commands that replicate Logsdon 2021 cen 8 dotplot figure
```
./cmds.sh chr8 ./test/chr8_cen.fasta 5000
# change the inputs in aln.R and then run
```

# Running on the full T2T CHM13 V1.0 assembly
Make `./dot_aln.yaml` with the following text:
```
sample: chm13.draft_v1.0
nbatch: 1000
alnthreads: 24
fasta: ../assemblies/chm13.draft_v1.0.fasta
window: 5000
```

Run the snakemake pipeline with:
```
./dot_aln.sh
```

