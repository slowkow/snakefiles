'''
star.snakefile
Kamil Slowikowski

Map paired-end RNA-seq reads and count genes with STAR
------------------------------------------------------

Requirements:

  samtools
      http://www.htslib.org/download/

  STAR
      https://github.com/alexdobin/STAR/releases

Usage: 

  snakemake \
  	--snakefile star.snakefile \
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
        'junctions.tsv.gz'

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
# 3. Count the number of reads mapped to each gene.
# 4. Count the number of reads supporting each splice junction.
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
              ' --sjdbFileChrStartEnd {input.sjs}'
              # Count fragments per gene, similar to HTseq.
              ' --quantMode GeneCounts'
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

