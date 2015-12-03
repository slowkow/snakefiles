# kallisto

Quantify gene isoform expression in transcripts per million (TPM) with
[kallisto] and collate outputs from multiple samples into one file.

[kallisto]: https://github.com/pachterlab/kallisto

## View the job graph

```bash
snakemake --forceall --dag | dot -Tpng > dag.png
```

![Snakemake directed acyclic graph (DAG).][dag]

[dag]: https://github.com/slowkow/snakefiles/blob/master/kallisto/dag.png

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
zcat abundance.tsv.gz | head | column -t
```

```
sample   target_id        length  eff_length  est_counts  tpm
Sample1  ENST00000390469  520     367.098     1           980.612
Sample1  ENST00000453496  2469    2313.85     1           155.576
Sample1  ENST00000620987  945     789.852     1.68005     765.696
Sample1  ENST00000614992  948     792.852     1.67574     760.841
Sample1  ENST00000633705  760     604.852     4.64421     2764.03
Sample1  ENST00000436911  1013    857.852     1           419.631
Sample1  ENST00000390290  400     252.302     1           1426.78
Sample1  ENST00000633188  51      8.05        0.5         22359.1
Sample1  ENST00000390413  51      8.05        0.5         22359.1
```

```bash
zcat n_processed.tsv.gz
```

```
sample  n_processed
Sample1 5000
Sample2 4500
```
