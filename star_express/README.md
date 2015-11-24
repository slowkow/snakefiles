# star_express.snakefile

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site. Also
produce a BAM file with coordinates relative to transcripts. Quantify
transcripts in TPM with [eXpress]. Collate outputs from multiple samples.

[STAR]: https://github.com/alexdobin/STAR
[eXpress]: http://bio.math.berkeley.edu/eXpress/overview.html

## View the job graph

```bash
snakemake \
  --snakefile star_express.snakefile \
  --configfile config.yml \
  --forceall \
  --dag \
  | dot -Tpng > dag.png
```

![Snakemake directed acyclic graph (DAG).][dag]

[dag]: https://github.com/slowkow/snakefiles/blob/master/star_express/dag.png

## Run the Snakefile

You can run snakemake like this. Notice the `--cluster 'bsub ...` option, used
to launch jobs on an [LSF] cluster.

[LSF]: https://en.wikipedia.org/wiki/Platform_LSF

```bash
snakemake \
  --snakefile star_express.snakefile \
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

```bash
zcat express.tsv.gz | awk '$6 > 0' | head | column -t
```

```
sample   bundle_id  target_id        length  eff_length   tot_counts  uniq_counts  est_counts  eff_counts  ambig_distr_alpha  ambig_distr_beta  fpkm          fpkm_conf_low  fpkm_conf_high  solvable  tpm
Sample2  78         ENST00000348428  9368    9170.331276  1           0            0.500000    0.510778    1.000000e+00       1.000000e+00      2.467134e+01  0.000000e+00   7.401401e+01    F         2.726641e+01
Sample2  78         ENST00000587011  625     438.998418   1           0            0.500000    0.711848    1.000000e+00       1.000000e+00      5.153648e+02  0.000000e+00   1.546094e+03    F         5.695738e+02
Sample2  320        ENST00000318006  4844    4652.205229  1           0            0.333333    0.347076    1.000000e+00       1.000000e+00      3.242109e+01  0.000000e+00   1.241218e+02    F         3.583132e+01
Sample2  320        ENST00000614932  2835    2645.813710  1           0            0.333333    0.357168    1.000000e+00       1.000000e+00      5.700687e+01  0.000000e+00   2.182467e+02    F         6.300318e+01
Sample2  320        ENST00000562258  4160    3969.093333  1           0            0.333333    0.349366    1.000000e+00       1.000000e+00      3.800101e+01  0.000000e+00   1.454841e+02    F         4.199817e+01
Sample2  357        ENST00000356978  4242    4050.986864  1           1            1.000000    1.047152    0.000000e+00       0.000000e+00      1.116984e+02  1.116984e+02   1.116984e+02    T         1.234475e+02
Sample2  437        ENST00000339818  2233    2044.595346  1           0            0.166667    0.182025    1.000000e+00       1.000000e+00      3.688494e+01  0.000000e+00   2.018394e+02    F         4.076471e+01
Sample2  437        ENST00000470196  654     467.906615   1           0            0.166667    0.232952    1.000000e+00       1.000000e+00      1.611749e+02  0.000000e+00   8.819707e+02    F         1.781281e+02
Sample2  437        ENST00000496321  1191    1003.948277  1           0            0.166667    0.197719    1.000000e+00       1.000000e+00      7.511819e+01  0.000000e+00   4.110570e+02    F         8.301955e+01
```

