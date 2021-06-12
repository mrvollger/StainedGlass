# StainedGlass
[![Actions Status](https://github.com/mrvollger/StainedGlass/workflows/CI/badge.svg)](https://github.com/mrvollger/StainedGlass/actions) 

This is a repository for making colorful dot-plots of genomic sequence.
![](images/chr8.png "chr8 cen")


## Installation 

You will need a current version of `snakemake` to run this workflow. To get `snakemake` please follow the install [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on their website, but in brief once `conda` and `mamba` are installed you can install `snakemake` with:
```
mamba create -n snakemake -c conda-forge -c bioconda snakemake
```
Afterwards you can activate the `conda` environment and download the repository. And all additional dependencies will be handled by `snakemake`.
```
conda activate snakemake
git clone https://github.com/mrvollger/StainedGlass.git
```

## Running
Choose a sample identifier for your run e.g. `chr8` and a fasta file on which you want to show the colorful alignments and the modify the config file `config/config.yaml` accordingly.

Once this is done this and activated your `conda` env with `snakemake` you can run the pipeline like so:
```
./StainedGlass.sh --cores 24 
```
Or do a dry run of the pipeline:
```
./StainedGlass.sh --cores 24 -n
```
In fact, all parameters passed to `StainedGlass.sh` are passed to `snakemake` so you can use all the `snakemake` command line options.

Please try the test case with the default configuration file before submitting issues.
If you are familiar with `snakemake` and want to trouble shoot yourself you can find the `Snakefile` in the directory `workflow`.

Descriptions of additional configuration parameters can be found in `config/README.md`.

### Making figures for a small number of regions
To make pdfs and pngs for a particular set of regions just add `make_figures` to your command.
```
./StainedGlass.sh --cores 24 make_figures
```

### Output
The file `results/{sample}.{\d+}.{\d+}.bed` will contain all the alignments identified by the pipeline, and is the main input for figure generation. Under the same prefix there will also be a bam file that contains the unprocessed alignments. Note the bam will contain additional alignments not present in the bed file because redundant alignments with lower scores are removed in the figures.

If you ran the pipeline with `make_figures` then there will be an an output directory under `results/{SM}.{\d+}.{\d+}_figures` with a variety of dot plots in `pdf` and `png` format. 

If you see `tri.TRUE` in the output pdf/png it means that the dot plot is rotated and cropped into a triangle. If you see `onecolorscale.FALSE` it means that between different facets in the same plot different color scales are being used. 


### Making a visualization for the whole genome
Making an interactive whole genome visualization requires the use of the program [HiGlass](https://higlass.io/) and a web browser. However, this pipeline will make the necessary input files with the following command:
```
./StainedGlass.sh --cores 24 cooler
```

To view locally, use `higlass-manage`:

```
pip install higlass-manage
higlass-manage view results/small.5000.10000.strand.mcool
```

See the [T2T CHM13 v1.0 StainedGlass](https://resgen.io/paper-data/T2T/views/MtjcVgrlQmymnHIvdck5-g) for an example. 
