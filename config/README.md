## Configuration options

Below are the options that can be changed in `config.yaml`.

A sample/prefix identifier to append to your result files

```
sample: test
```

Path to a fasta file to deploy the workflow on

```
fasta: .test/small.fasta
```

Size of the window in which to breakup the input fasta before all by all alignment.

```
window: 5000
```

Setting for the minimap2 `-f` parameter. A smaller number will increase sensitivity at the cost of runtime.
See the `minimap2` man page for more details.

```
mm_f: 10000
```

The number of alignment jobs to distribute the workflow across. Does not change final output.

```
nbatch: 1
```

The number of alignment threads per job.

```
alnthreads: 4
```

Path for a temp dir to be used by pipeline.

```
tempdir: temp
```

### Configuration options for high-res dot plots (cooler_density)

This defines the smallest bin size (highest resolution) used in the cooler file.

```
cooler_window: 100
```

Since `window` defines the read length, it should be made to be smaller than `cooler_window`. I like to use `window: 32` and `cooler_window: 100`.

This is the maximum number of alignments to output for each read. From the `bwa` help page: Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written.

```
num_dups: 100
```
