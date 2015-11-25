'''
star_express.snakefile
Kamil Slowikowski

Map paired-end RNA-seq reads with STAR and quantify transcripts with eXpress
----------------------------------------------------------------------------

Requirements:

  samtools
      http://www.htslib.org/download/

  STAR
      https://github.com/alexdobin/STAR/releases

  gffread
      https://cole-trapnell-lab.github.io/cufflinks/install/

  express
      http://bio.math.berkeley.edu/eXpress/

Usage: 

  snakemake \
  	--snakefile star_express.snakefile \
  	--configfile config.yml \
  	--jobs 999 \
  	--cluster 'bsub -q big-multi -n 16 -R "rusage[mem=35000]"'
'''

from os.path import join, basename, dirname
from subprocess import check_output

# Globals ---------------------------------------------------------------------

# The config dictionary contains values defined in the JSON or YAML file
# specified with 'snakemake --configfile'.

# Run multi-threaded programs with this many threads.
THREADS = config['THREADS']

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']

# Full path to an uncompressed GTF file with all gene annotations.
GTF = config['GTF']

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config['FASTQ_DIR']

# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']

# A snakemake regular expression matching the forward mate FASTQ files.
# e.g. '{sample,[^/]+}.R1.fastq.gz'
SAMPLES, = glob_wildcards(join(FASTQ_DIR, config['PATTERN_SAMPLES']))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
# e.g. '{sample}.R1.fastq.gz'
PATTERN_R1 = config['PATTERN_R1']
PATTERN_R2 = config['PATTERN_R2']

STAR_VERSION = check_output('STAR | grep versionSTAR', shell=True).split()[1]

GFFREAD_VERSION = check_output('md5sum $(which gffread)', shell=True)

EXPRESS_VERSION = check_output('echo $(express 2>&1 | grep "express v")', shell=True).split()[1]

# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

rule all:
    input:
        'counts.tsv.gz',
        'junctions.tsv.gz',
        'express.tsv.gz'

# Make a CDNA file from the given GTF and DNA files.
rule make_cdna:
    input:
        dna = DNA,
        gtf = GTF
    output:
        cdna = join(dirname(DNA), 'gffread',
                    rstrip(DNA, '.fa') + '.gffread_transcripts.fa'),
        faidx = DNA + '.fai'
    log:
        join(dirname(DNA), 'gffread', 'gffread.log')
    benchmark:
        join(dirname(DNA), 'gffread', 'gffread.benchmark.tsv')
    version:
        GFFREAD_VERSION
    run:
        # Create a FASTA index file for gffread.
        shell('samtools faidx {input.dna}')
        # Extract a sequence for each transcript in the GTF file.
        cmd = 'gffread -F -w {output.cdna} -g {input.dna} {input.gtf}'
        # Log the command and its output.
        shell('echo ' + cmd + ' > {log}')
        shell(cmd + ' >> {log} 2>&1')

# Make an index of the genome for STAR.
rule star_index:
    input:
        dna = DNA
    output:
        index = join(dirname(DNA), 'star', 'Genome')
    log:
        join(dirname(DNA), 'star', 'star.index.log')
    benchmark:
        join(dirname(DNA), 'star', 'star.index.benchmark.tsv')
    version:
        STAR_VERSION
    threads:
        THREADS
    run:
        # Write stderr and stdout to the log file.
        shell('STAR'
              ' --runThreadN {threads}'
              ' --runMode genomeGenerate'
              ' --genomeDir ' + join(dirname(DNA), 'star') +
              ' --genomeFastaFiles {input.dna}'
              ' > {log} 2>&1')

# 1. Map paired-end RNA-seq reads to the genome.
# 2. Count the number of reads supporting each splice junction.
# 3. Delete the output SAM file.
rule star_pass1:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        genomeDir = dirname(rules.star_index.output.index),
        gtf = GTF
    output:
        sam = temp(join(OUT_DIR, '{sample}', 'pass1', 'Aligned.out.sam')),
        sj = join(OUT_DIR, '{sample}', 'pass1', 'SJ.out.tab')
    log:
        join(OUT_DIR, '{sample}', 'pass1', 'star.map.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'pass1', 'star.map.benchmark.tsv')
    version:
        STAR_VERSION
    threads:
        THREADS
    run:
        # Map reads with STAR.
        shell('cd ' + join(OUT_DIR, '{wildcards.sample}', 'pass1') +
              '&& STAR'
              ' --runThreadN {threads}'
              ' --genomeDir ' + join(dirname(DNA), 'star') +
              ' --sjdbGTFfile {input.gtf}'
              ' --readFilesCommand zcat'
              ' --readFilesIn {input.r1} {input.r2}'
              # By default, this prefix is "./".
              ' --outFileNamePrefix ' + join(OUT_DIR, '{wildcards.sample}', 'pass1') + '/'
              #
              # This option causes STAR to throw an error.
              # Output sorted by coordinate.
              # ' --outSAMtype BAM SortedByCoordinate'
              #
              # If exceeded, the read is considered unmapped.
              ' --outFilterMultimapNmax 20'
              # Minimum overhang for unannotated junctions.
              ' --alignSJoverhangMin 8'
              # Minimum overhang for annotated junctions.
              ' --alignSJDBoverhangMin 1'
              # Maximum number of mismatches per pair.
              ' --outFilterMismatchNmax 999'
              # Minimum intron length.
              ' --alignIntronMin 1'
              # Maximum intron length.
              ' --alignIntronMax 1000000'
              # Maximum genomic distance between mates.
              ' --alignMatesGapMax 1000000'
              ' > {log} 2>&1')

# 1. Map paired-end RNA-seq reads to the genome.
# 2. Make a coordinate sorted BAM with genomic coordinates.
# 3. Make a second BAM with coordinates relative to transcripts.
# 4. Count the number of reads mapped to each gene.
# 5. Count the number of reads supporting each splice junction.
rule star_pass2:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        genomeDir = dirname(rules.star_index.output.index),
        gtf = GTF,
        sjs = expand(join(OUT_DIR, '{sample}', 'pass1', 'SJ.out.tab'), sample = SAMPLES)
    output:
        sam = temp(join(OUT_DIR, '{sample}', 'pass2', 'Aligned.out.sam')),
        bam = join(OUT_DIR, '{sample}', 'pass2', 'Aligned.out.bam'),
        t_bam = temp(join(OUT_DIR, '{sample}', 'pass2', 'Aligned.toTranscriptome.out.bam')),
        t_bam_sorted = join(OUT_DIR, '{sample}', 'pass2', 'Aligned.toTranscriptome.out.sorted.bam'),
        counts = join(OUT_DIR, '{sample}', 'pass2', 'ReadsPerGene.out.tab'),
        sj = join(OUT_DIR, '{sample}', 'pass2', 'SJ.out.tab')
    log:
        join(OUT_DIR, '{sample}', 'pass2', 'star.map.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'pass2', 'star.map.benchmark.tsv')
    version:
        STAR_VERSION
    threads:
        THREADS
    run:
        # Map reads with STAR.
        shell('cd ' + join(OUT_DIR, '{wildcards.sample}', 'pass2') +
              '&& STAR'
              ' --runThreadN {threads}'
              ' --genomeDir ' + join(dirname(DNA), 'star') +
              ' --sjdbGTFfile {input.gtf}'
              ' --readFilesCommand zcat'
              ' --readFilesIn {input.r1} {input.r2}'
              # Multi-sample 2-pass alignment, sharing splice junctions across
              # samples.
              ' --sjdbFileChrStartEnd {input.sjs}'
              # BAM file in transcript coords, in addition to genomic BAM file.
              ' --quantMode TranscriptomeSAM GeneCounts'
              # Allow insertions, deletions and soft-clips in the transcriptomic
              # alignments, which can be used by eXpress.
              ' --quantTranscriptomeBan Singleend'
              # By default, this prefix is "./".
              ' --outFileNamePrefix ' + join(OUT_DIR, '{wildcards.sample}', 'pass2') + '/'
              #
              # This option causes STAR to throw an error.
              # Output sorted by coordinate.
              # ' --outSAMtype BAM SortedByCoordinate'
              #
              # If exceeded, the read is considered unmapped.
              ' --outFilterMultimapNmax 20'
              # Minimum overhang for unannotated junctions.
              ' --alignSJoverhangMin 8'
              # Minimum overhang for annotated junctions.
              ' --alignSJDBoverhangMin 1'
              # Maximum number of mismatches per pair.
              ' --outFilterMismatchNmax 999'
              # Minimum intron length.
              ' --alignIntronMin 1'
              # Maximum intron length.
              ' --alignIntronMax 1000000'
              # Maximum genomic distance between mates.
              ' --alignMatesGapMax 1000000'
              ' > {log} 2>&1')
        # Convert to BAM and sort by coordinate.
        shell('samtools view -b {output.sam}'
              ' | samtools sort -@{threads} -l9 - ' + rstrip(output.bam, '.bam'))
        # Sort by read name for eXpress.
        shell('samtools sort -@{threads} -l9 -n'
              ' -o {output.t_bam_sorted}'
              ' -O bam'
              ' -T ' + rstrip(output.t_bam, '.bam') +
              ' {output.t_bam}')

# Compute transcripts per million (TPM) with eXpress.
rule express:
    input:
        bam = rules.star_pass2.output.t_bam_sorted,
        cdna = rules.make_cdna.output.cdna
    output:
        results = join(OUT_DIR, '{sample}', 'pass2', 'results.xprs')
    log:
        join(OUT_DIR, '{sample}', 'pass2', 'express.log')
    benchmark:
        join(OUT_DIR, '{sample}', 'pass2', 'express.benchmark.tsv')
    version:
        EXPRESS_VERSION 
    run:
        shell('express'
              ' --no-update-check'
              ' --no-bias-correct'
              ' --logtostderr'
              ' --output-dir ' + join(OUT_DIR, '{wildcards.sample}', 'pass2') +
              ' {input.cdna} {input.bam}'
              ' > {log} 2>&1')

rule collate_express:
    input:
        expand(join(OUT_DIR, '{sample}', 'pass2', 'results.xprs'), sample = SAMPLES)
    output:
        'express.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(dirname(i)))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                    out.write(b(sample + '\t' + line))

rule collate_counts:
    input:
        expand(join(OUT_DIR, '{sample}', 'pass2', 'ReadsPerGene.out.tab'), sample = SAMPLES)
    output:
        'counts.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = [b'sample', b'gene_id',
                b'counts_unstranded', b'counts_strand1', b'counts_strand2']
            out.write(b'\t'.join(header) + b'\n')

            for i in input:
                sample = basename(dirname(dirname(i)))
                for line in open(i):
                    out.write(b(sample + '\t' + line))

rule collate_junctions:
    input:
        expand(join(OUT_DIR, '{sample}', 'pass2', 'SJ.out.tab'), sample = SAMPLES)
    output:
        'junctions.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = [b'sample', b'chrom', b'intron_first', b'intron_last',
                b'strand', b'intron_motif', b'annotated', b'uniquely_mapped_reads',
                b'multimapped_reads', b'max_spliced_overhang']
            out.write(b'\t'.join(header) + b'\n')

            for i in input:
                sample = basename(dirname(dirname(i)))
                for line in open(i):
                    out.write(b(sample + '\t' + line))

