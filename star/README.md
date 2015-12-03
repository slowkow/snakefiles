# star

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site.
Collate outputs from multiple samples.

[STAR]: https://github.com/alexdobin/STAR

## View the job graph

```bash
snakemake --forceall --dag | dot -Tpng > dag.png
```

![Snakemake directed acyclic graph (DAG).][dag]

[dag]: https://github.com/slowkow/snakefiles/blob/master/star/dag.png

## Run the Snakefile

You can run snakemake like this. Notice the `--cluster` option, used to launch
jobs on an [LSF] cluster. The Python script [bsub.py][bsub] receives job
scripts from Snakemake and automatically submits them to an appropriate LSF
queue based on job requirements.

[LSF]: https://en.wikipedia.org/wiki/Platform_LSF
[bsub]: https://github.com/slowkow/snakefiles/tree/master/bsub.py

```bash
snakemake --jobs 999 --cluster '../bsub.py -o stdout'
```

## Output

The collated output looks like this:

```bash
zcat counts.tsv.gz | head | column -t
```

```
sample   gene_id          counts_unstranded  counts_strand1  counts_strand2
Sample1  N_unmapped       571                571             571
Sample1  N_multimapping   1384               1384            1384
Sample1  N_noFeature      1178               2059            2087
Sample1  N_ambiguous      141                33              36
Sample1  ENSG00000225630  1                  0               1
Sample1  ENSG00000237973  1                  1               0
Sample1  ENSG00000248527  1                  0               1
Sample1  ENSG00000069424  1                  0               1
Sample1  ENSG00000116251  1                  1               0
```

```bash
zcat junctions.tsv.gz | head | column -t
```

```
sample   chrom  intron_first  intron_last  strand  intron_motif  annotated  uniquely_mapped_reads  multimapped_reads  max_spliced_overhang
Sample1  1      6095444       6095529      1       1             1          1                      0                  19
Sample1  1      6095625       6096635      1       1             1          1                      0                  30
Sample1  1      10336743      10337073     1       1             1          1                      0                  33
Sample1  1      16206084      16206179     2       2             1          1                      0                  24
Sample1  1      16206341      16206947     2       2             1          1                      0                  8
Sample1  1      23691830      23692608     1       1             1          1                      0                  22
Sample1  1      23692760      23693806     1       1             1          1                      0                  31
Sample1  1      23804526      23808134     2       2             1          1                      0                  12
Sample1  1      27350133      27350388     1       1             1          1                      0                  28
```
