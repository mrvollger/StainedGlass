## Configuration options
Below are the options that can be changed in `config.yaml`.

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