'''
kallisto.snakefile
Kamil Slowikowski

Quantify transcript expression in paired-end RNA-seq data with kallisto
-----------------------------------------------------------------------

Requirements:

  kallisto
      https://pachterlab.github.io/kallisto/download.html

Usage: 

  snakemake \
  	--snakefile kallisto.snakefile \
  	--configfile config.yml \
  	--jobs 999 \
  	--cluster 'bsub -q normal -R "rusage[mem=4000]"'
'''

from os.path import join, basename, dirname
from subprocess import check_output

# Globals ---------------------------------------------------------------------

# The config dictionary contains values defined in the JSON or YAML file
# specified with 'snakemake --configfile'.

# Run multi-threaded programs with this many threads.
THREADS = config['THREADS']

# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['CDNA']

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

KALLISTO_VERSION = check_output("echo $(kallisto)", shell=True).split()[1]

# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

rule all:
    input:
        'abundance.tsv.gz'

rule kallisto_index:
    input:
        cdna = CDNA
    output:
        index = join(dirname(CDNA), 'kallisto', rstrip(basename(CDNA), '.fa'))
    log:
        join(dirname(CDNA), 'kallisto/kallisto.index.log')
    benchmark:
        join(dirname(CDNA), 'kallisto/kallisto.index.benchmark.tsv')
    version:
        KALLISTO_VERSION
    run:
        # Record the kallisto version number in the log file.
        shell('echo $(kallisto index) &> {log}')
        # Write stderr and stdout to the log file.
        shell('kallisto index'
              ' --index={output.index}'
              ' --kmer-size=21'
              ' --make-unique'
              ' {input.cdna}'
              ' >> {log} 2>&1')

rule kallisto_quant:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        index = rules.kallisto_index.output.index
    output:
        join(OUT_DIR, '{sample}/abundance.tsv')
    version:
        KALLISTO_VERSION
    threads:
        THREADS
    shell:
        'kallisto quant'
        ' --threads={threads}'
        ' --index={input.index}'
        ' --output-dir=' + join(OUT_DIR, '{wildcards.sample}') +
        ' {input.r1} {input.r2}'

rule collate_kallisto:
    input:
        expand(join(OUT_DIR, '{sample}/abundance.tsv'), sample=SAMPLES)
    output:
        'abundance.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(i))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                    out.write(b(sample + '\t' + line))

