# star_express

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site. Also
produce a BAM file with coordinates relative to transcripts. Quantify
transcripts in TPM with [eXpress]. Collate outputs from multiple samples.

[STAR]: https://github.com/alexdobin/STAR
[eXpress]: http://bio.math.berkeley.edu/eXpress/overview.html

## View the job graph

```bash
snakemake --forceall --dag | dot -Tpng > dag.png
```

![Snakemake directed acyclic graph (DAG).][dag]

[dag]: https://github.com/slowkow/snakefiles/blob/master/star_express/dag.png

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

```bash
zcat express.tsv.gz | head | column -t
```

```
sample   bundle_id  target_id        length  eff_length    tot_counts  uniq_counts  est_counts  eff_counts  ambig_distr_alpha  ambig_distr_beta  fpkm          fpkm_conf_low  fpkm_conf_high  solvable  tpm
Sample1  381        ENST00000216281  3379    3204.675807   1           0            0.500000    0.527198    1.000000e+00       1.000000e+00      7.436703e+01  0.000000e+00   2.231011e+02    F         8.208808e+01
Sample1  381        ENST00000334701  3510    3335.674784   1           0            0.500000    0.526130    1.000000e+00       1.000000e+00      7.144648e+01  0.000000e+00   2.143394e+02    F         7.886431e+01
Sample1  1594       ENST00000215832  11022   10847.616110  2           2            2.000000    2.032152    0.000000e+00       0.000000e+00      8.788003e+01  8.788003e+01   8.788003e+01    T         9.700405e+01
Sample1  1750       ENST00000310144  1579    1404.689867   1           0            0.500000    0.562046    1.000000e+00       1.000000e+00      1.696618e+02  0.000000e+00   5.089854e+02    F         1.872767e+02
Sample1  1750       ENST00000582130  828     653.695733    1           0            0.500000    0.633322    1.000000e+00       1.000000e+00      3.645767e+02  0.000000e+00   1.093730e+03    F         4.024283e+02
Sample1  1762       ENST00000616577  4814    4639.664599   4           0            0.280499    0.291039    1.000000e+00       1.000000e+00      2.881640e+01  0.000000e+00   2.386844e+02    F         3.180822e+01
Sample1  1762       ENST00000530705  4735    4560.665216   4           0            0.310656    0.322531    1.000000e+00       1.000000e+00      3.246727e+01  0.000000e+00   2.562418e+02    F         3.583814e+01
Sample1  1762       ENST00000379056  1101    926.693600    4           0            2.230352    2.649870    1.000000e+00       1.000000e+00      1.147181e+03  0.000000e+00   3.190887e+03    F         1.266285e+03
Sample1  1762       ENST00000528619  731     556.858785    1           0            0.000796    0.001044    1.000000e+00       1.000000e+00      6.809137e-01  0.000000e+00   2.007712e+01    F         7.516087e-01
```
