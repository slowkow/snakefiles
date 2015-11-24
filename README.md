# snakefiles

This repository has Snakefiles for common RNA-seq data analysis workflows.
Please feel free to copy them and modify them to suit your needs.

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
