# Commands that replicate Logsdon 2021 cen 8 dotplot figure
```
./cmds.sh chr8 ./test/chr8_cen.fasta 5000
# change the inputs in aln.R and then run
```

# Running on the full T2T CHM13 V1.0 assembly
Make `./dot_aln.yaml` with the following feilds:
```
sample: chm13.draft_v1.0  # a prefix for your output
nbatch: 1000              # number of batches to split your input alignments across
alnthreads: 24            # alignment threads per batch
window: 5000              # window size for alignments 
fasta: ../assemblies/chm13.draft_v1.0.fasta
```

Run the snakemake pipeline with:
```
./dot_aln.sh
```

Unless you need to use resgen.io ignore the `cooler` rules that fail.

Once the pipeline is done open `./aln_plot.R` and look for the section 
where it tells you to edit the inputs. Then run the script in Rstudio.

