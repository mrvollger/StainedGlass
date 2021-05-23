# StainedGlass

This is a repository for making colorful dot-plots of genomic sequence.
![](test/chr13_p_arm.png "chr13 P arm")


## Installation 

You will need a current version of `snakemake` to run the code. To get `snakemake` please follow the install [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on their website.

Afterwards you can download the repository and all additional dependencies will be handled by `snakemake`.
```
git clone https://github.com/mrvollger/StainedGlass.git
```

## Running
Choose a sample identifier for your run e.g. `chr8` and a fasta file on which you want to show the colorful alignments and the modify the config file `StainedGlass.yaml` accordingly.
```
sample: small
fasta: test/small.fasta
```

Once this is done you can run the pipeline like so:
```
./StainedGlass.sh --cores 24 
```
Or do a dry run of the pipeline:
```
./StainedGlass.sh --cores 24 -n
```
In fact, all parameters passed to `StainedGlass.sh` are passed to `snakemake` so you can use all the `snakemake` command line options.

Please try the test case with the default configuration file before submitting issues.

### Making figures for a small number of regions
To make pdfs and pngs for a particular set of regions just add `make_figures` to your command.
```
./StainedGlass.sh --cores 24 make_figures
```

### output
The file `results/{sample}.{\d+}.{\d+}.bed` will contain all the alignments identified by the pipeline, and is the main input for figure generation. Under the same prefix there will also be a bam file that contains the unprocessed alignments. Note the bam will contain additional alignments not present in the bed file because redundant alignments with lower scores are removed in the figures.

If you ran the pipeline with `make_figures` then there will be an an output directory under `results/{SM}.{\d+}.{\d+}_figures` with a variety of dot plots in `pdf` and `png` format. 

If you see `tri.TRUE` in the output pdf/png it means that the dot plot is rotated and cropped into a triangle. If you see `onecolorscale.FALSE` it means that between different facets in the same plot different color scales are being used. 


### Making a visualization for the whole genome
This requires the use of the program [HiGlass](https://higlass.io/) and a web browser. But to make the necessary input files for that you can run:
```
./StainedGlass.sh --cores 24 cooler
```
