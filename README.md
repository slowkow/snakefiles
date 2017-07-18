# snakefiles

This repository has Snakefiles for common RNA-seq data analysis workflows.
Please feel free to copy them and modify them to suit your needs.


## Getting started

If you are new to [Snakemake], you might like to start by walking through my
[tutorial for beginners][beginners]. Next, have a look at Johannes Koster's
[introductory slides][slides], [tutorial], [documentation], and [FAQ].

Quick start:

```bash
# Copy the files
git clone https://github.com/slowkow/snakefiles.git

# Go to the kallisto directory
cd snakefiles/kallisto

# Run snakemake
snakemake
```

## Data

This repository includes 6 FASTQ files in [data/fastq/][fastq] to illustrate
the usage of each of the RNA-seq workflows.

- Sample1
    - `Sample1.R1.fastq.gz` has the first mates of sequenced fragments.
    - `Sample1.R2.fastq.gz` has the second mates of sequenced fragments.
- Sample2
    - `Sample2.L1.R1.fastq.gz`
    - `Sample2.L2.R1.fastq.gz`
        - The first mate reads (R1), split across two files (L1 and L2). Some
          software such as STAR requires these reads to be merged into one file.
    - `Sample2.L1.R2.fastq.gz`
    - `Sample2.L2.R2.fastq.gz`
        - Likewise, the second mate reads (R2) are also split across two files
          (L1 and L2). To make matters worse, `Sample2.L2.R2.fastq.gz` has only
          2000 reads, whereas `Sample2.L2.R1.fastq.gz` has 2500 reads. The
          Snakefiles in this repository can handle this without any problems.

[fastq]: https://github.com/slowkow/snakefiles/tree/master/data/fastq


## Scripts

- [make_samples.py][make_samples] creates the [samples.json][samples] file.
- [bsub.py][bsub] receives job scripts from Snakemake and automatically
  submits them to an appropriate LSF queue based on job requirements.

[make_samples]: https://github.com/slowkow/snakefiles/tree/master/make_samples.py
[samples]: https://github.com/slowkow/snakefiles/tree/master/samples.json
[bsub]: https://github.com/slowkow/snakefiles/tree/master/bsub.py


## RNA-seq workflows


### [kallisto/][1]

[1]: https://github.com/slowkow/snakefiles/tree/master/kallisto

Quantify gene isoform expression in transcripts per million (TPM) with
[kallisto] and collate outputs from multiple samples into one file.


### [star_express/][3]

[3]: https://github.com/slowkow/snakefiles/tree/master/star_express

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site. Also
produce a BAM file with coordinates relative to transcripts. Quantify
transcripts in TPM with [eXpress]. Collate outputs from multiple samples.


## Contributing

Please [submit an issue][issues] to report bugs or ask questions.

Please contribute bug fixes or new features with a [pull request][pull] to
this repository.

[issues]: https://github.com/slowkow/snakefiles/issues
[pull]: https://help.github.com/articles/using-pull-requests/

[kallisto]: https://github.com/pachterlab/kallisto
[STAR]: https://github.com/alexdobin/STAR
[eXpress]: http://bio.math.berkeley.edu/eXpress/overview.html
[cc0]: https://creativecommons.org/publicdomain/zero/1.0/

[beginners]: http://slowkow.com/notes/snakemake-tutorial/

[Snakemake]: https://bitbucket.org/snakemake/snakemake/wiki/Home
[slides]: http://slides.com/johanneskoester/deck-1
[tutorial]: http://htmlpreview.github.io/?https://bitbucket.org/snakemake/snakemake/raw/master/snakemake-tutorial.html
[documentation]: https://bitbucket.org/snakemake/snakemake/wiki/Documentation
[FAQ]: https://bitbucket.org/snakemake/snakemake/wiki/FAQ
