# star.snakefile

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site.
Collate outputs from multiple samples.

[STAR]: https://github.com/alexdobin/STAR

## View the job graph

```bash
snakemake \
  --snakefile star.snakefile \
  --configfile config.yml \
  --forceall \
  --dag \
  | dot -Tpng > dag.png
```

![Snakemake directed acyclic graph (DAG).][dag]

[dag]: https://github.com/slowkow/snakefiles/blob/master/star/dag.png

## Run the Snakefile

You can run snakemake like this. Notice the `--cluster 'bsub ...` option, used
to launch jobs on an [LSF] cluster.

[LSF]: https://en.wikipedia.org/wiki/Platform_LSF

```bash
snakemake \
  --snakefile star.snakefile \
  --configfile config.yml \
  --jobs 999 \
  --cluster 'bsub -q big-multi -n 16 -R "rusage[mem=35000]"'
```

## Output

The collated output looks like this:

```bash
zcat counts.tsv.gz | grep N_ | column -t
```

```
Sample2  N_unmapped      480   480   480
Sample2  N_multimapping  1653  1653  1653
Sample2  N_noFeature     973   1853  1889
Sample2  N_ambiguous     155   37    31
Sample1  N_unmapped      571   571   571
Sample1  N_multimapping  1383  1383  1383
Sample1  N_noFeature     1178  2059  2088
Sample1  N_ambiguous     142   34    36
```

```bash
zcat junctions.tsv.gz | head | column -t
```

```
sample   chrom  intron_first  intron_last  strand  intron_motif  annotated  uniquely_mapped_reads  multimapped_reads  max_spliced_overhang
Sample2  1      2056574       2059540      1       1             1          1                      0                  29
Sample2  1      6462175       6462628      2       2             1          1                      0                  24
Sample2  1      6462667       6462862      2       2             1          1                      0                  30
Sample2  1      8861430       8862886      2       2             1          1                      0                  10
Sample2  1      9570359       9573345      1       1             1          1                      0                  20
Sample2  1      9573413       9579953      1       1             1          1                      0                  4
Sample2  1      11080550      11080763     2       2             1          1                      0                  22
Sample2  1      25242706      25243549     2       2             1          1                      0                  21
Sample2  1      26318072      26320170     1       1             1          1                      0                  30
```
