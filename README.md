# snakefiles

This repository has Snakefiles for common RNA-seq data analysis workflows.
Please feel free to copy them and modify them to suit your needs.

Unless stated otherwise, everything in this repository is dedicated to the
public domain through the [Creative Commons Zero license][cc0].

If you are new to [Snakemake], you might like to start by walking through my
[tutorial for beginners][beginners]. Next, have a look at Johannes Koster's
[introductory slides][slides], [tutorial], [documentation], and [FAQ].


### [kallisto/][1]

[1]: https://github.com/slowkow/snakefiles/tree/master/kallisto

Quantify gene isoform expression in transcripts per million (TPM) with
[kallisto] and collate outputs from multiple samples into one file.


### [star/][2]

[2]: https://github.com/slowkow/snakefiles/tree/master/star

Execute a multi-sample 2-pass [STAR] alignment, sharing the splice junctions
across samples. Count fragments per gene and fragments per splice site.
Collate outputs from multiple samples.


[kallisto]: https://github.com/pachterlab/kallisto
[STAR]: https://github.com/alexdobin/STAR
[cc0]: https://creativecommons.org/publicdomain/zero/1.0/

[beginners]: http://slowkow.com/notes/snakemake-tutorial/

[Snakemake]: https://bitbucket.org/snakemake/snakemake/wiki/Home
[slides]: http://slides.com/johanneskoester/deck-1
[tutorial]: http://htmlpreview.github.io/?https://bitbucket.org/snakemake/snakemake/raw/master/snakemake-tutorial.html
[documentation]: https://bitbucket.org/snakemake/snakemake/wiki/Documentation
[FAQ]: https://bitbucket.org/snakemake/snakemake/wiki/FAQ
