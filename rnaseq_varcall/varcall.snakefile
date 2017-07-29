# @author: Jean Fan
# @date: November 23, 2015
# @desc: Snakemake pipeline for https://www.broadinstitute.org/gatk/guide/article?id=3891
#        Based off of code by https://github.com/slowkow/snakefiles/blob/master/star/star.snakefile
#
# Usage: snakemake --snakefile varcall.snakefile --jobs 999 --cluster 'bsub -q short -W 12:00 -R "rusage[mem=4000]"'


from os.path import join, dirname
from subprocess import check_output


# Globals ---------------------------------------------------------------------

# Full path to genome fasta.
GENOME = "/n/data1/hms/dbmi/park/sl325/forJean/resources/human_g1k_v37_decoy.fasta"

# Full path to gene model annotations for splice aware alignment.
GTF = "/groups/shared_databases/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = "fastq/"

# Full path to output folder.
OUTPUT_DIR = "varcall/"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# A snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,Samp[^/]+}.R1.fastq.gz'))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = '{sample}.R1.fastq.gz'
PATTERN_R2 = '{sample}.R2.fastq.gz'


# Rules -----------------------------------------------------------------------

rule all:
    input:
        join(dirname(GENOME), 'star', 'Genome'),
        expand(join(OUTPUT_DIR, '{sample}', 'pass1', 'Aligned.out.sam'), sample = SAMPLES),
        expand(join(OUTPUT_DIR, '{sample}', 'pass2', 'Aligned.out.sam'), sample = SAMPLES),
        expand(join(OUTPUT_DIR, '{sample}', 'Aligned.out.dedupped.bam'), sample = SAMPLES),
        expand(join(OUTPUT_DIR, '{sample}', 'Aligned.out.split.bam'), sample = SAMPLES),
        expand(join(OUTPUT_DIR, '{sample}', 'call.filtered.vcf'), sample = SAMPLES)

# Make an index of the genome for STAR.
rule star_index:
    input:
        genome = GENOME
    output:
        index = join(dirname(GENOME), 'star', 'Genome')
    log:
        join(dirname(GENOME), 'star', 'star.index.log')
    threads:
        THREADS
    run:
        # Write stderr and stdout to the log file.
        shell('STAR'
              ' --runThreadN {threads}'
              ' --runMode genomeGenerate'
              ' --genomeDir ' + join(dirname(GENOME), 'star') +
              ' --genomeFastaFiles {input.genome}'
              ' > {log} 2>&1')

# 1. Map paired-end RNA-seq reads to the genome.
# 2. Count the number of reads supporting each splice junction.
rule star_pass1:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        genomeDir = dirname(rules.star_index.output.index),
        gtf = GTF
    output:
        sam = join(OUTPUT_DIR, '{sample}', 'pass1', 'Aligned.out.sam'),
        sj = join(OUTPUT_DIR, '{sample}', 'pass1', 'SJ.out.tab')
    log:
        join(OUTPUT_DIR, '{sample}', 'pass1', 'star.map.log')
    threads:
        THREADS
    run:
        # Map reads with STAR.
        shell('STAR'
              ' --runThreadN {threads}'
              ' --genomeDir ' + join(dirname(GENOME), 'star') +
              ' --sjdbGTFfile {input.gtf}'
              ' --readFilesCommand zcat'
              ' --readFilesIn {input.r1} {input.r2}'
              # By default, this prefix is "./".
              ' --outFileNamePrefix ' + join(OUTPUT_DIR, '{wildcards.sample}', 'pass1') + '/'
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
        sjs = expand(join(OUTPUT_DIR, '{sample}', 'pass1', 'SJ.out.tab'), sample = SAMPLES)
    output:
        sam = join(OUTPUT_DIR, '{sample}', 'pass2', 'Aligned.out.sam'),
        counts = join(OUTPUT_DIR, '{sample}', 'pass2', 'ReadsPerGene.out.tab'),
        sj = join(OUTPUT_DIR, '{sample}', 'pass2', 'SJ.out.tab')
    log:
        join(OUTPUT_DIR, '{sample}', 'pass2', 'star.map.log')
    threads:
        THREADS
    run:
        # Map reads with STAR.
        shell('STAR'
              ' --runThreadN {threads}'
              ' --genomeDir ' + join(dirname(GENOME), 'star') +
              ' --sjdbGTFfile {input.gtf}'
              ' --readFilesCommand zcat'
              ' --readFilesIn {input.r1} {input.r2}'
              ' --sjdbFileChrStartEnd {input.sjs}'
              # BAM file in transcript coords, in addition to genomic BAM file.
              ' --quantMode GeneCounts'
              # Basic 2-pass mapping, with all 1st pass junctions inserted
              # into the genome indices on the fly.
              ' --twopassMode Basic'
              # By default, this prefix is "./".
              ' --outFileNamePrefix ' + join(OUTPUT_DIR, '{wildcards.sample}', 'pass2') + '/'
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

# add read groups, sort
rule picard_cleanstar:
    input:
        sam = rules.star_pass2.output.sam
    output:
        bam = join(OUTPUT_DIR, '{sample}', 'Aligned.out.bam'),
    log:
        join(OUTPUT_DIR, '{sample}', 'picard.clean.log')
    run:
        shell('java -jar /opt/picard-1.138/bin/picard.jar'
              ' AddOrReplaceReadGroups'
              ' I={input.sam}'
              ' O={output}'
              ' SORT_ORDER=coordinate CREATE_INDEX=false RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample'
              ' > {log} 2>&1')

# mark duplicates
rule picard_markdup:
    input:
        bam = rules.picard_cleanstar.output.bam
    output:
        dupmarked = join(OUTPUT_DIR, '{sample}', 'Aligned.out.dedupped.bam'),
        metrics = join(OUTPUT_DIR, '{sample}', 'Aligned.out.dedupped.metrics'),
    log:
        join(OUTPUT_DIR, '{sample}', 'picard.markdup.log')
    run:
        shell('java -jar /opt/picard-1.138/bin/picard.jar'
              ' MarkDuplicates'
              ' I={input.bam}'
              ' O={output.dupmarked}'
              ' CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.metrics}'
              ' > {log} 2>&1')  


rule gatk_SplitNCigarReads:
    input:
        dupmarked = rules.picard_markdup.output.dupmarked,
        genome = GENOME
    output:
        split = join(OUTPUT_DIR, '{sample}', 'Aligned.out.split.bam')
    log:
        join(OUTPUT_DIR, '{sample}', 'gatk.splitcigar.log')
    run:
        shell('/opt/java/jdk1.7.0_71/bin/java -Xmx16g -jar /opt/gatk-3.4-46/gatk/GenomeAnalysisTK.jar'
              ' -T SplitNCigarReads'
              ' -R {input.genome}'
              ' -I {input.dupmarked}'
              ' -o {output.split}'
              ' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
              ' > {log} 2>&1')  


# # indel realignment

# rule gatk_realign_info:
#     input:
#         "mapping/{reference}/{prefix}.bam.bai",
#         ref=_get_ref,
#         bam="mapping/{reference}/{prefix}.bam"
#     output:
#         temp("mapping/{reference,[^/]+}/{prefix}.realign.intervals")
#     params:
#         custom=config.get("params_gatk", "")
#     log:
#         "mapping/log/{reference}/{prefix}.realign_info.log"
#     threads: 8
#     shell:
#         "gatk -T RealignerTargetCreator -R {input.ref} {params.custom} "
#         "-nt {threads} "
#         "-I {input.bam} -known {config[known_variants][dbsnp]} "
#         "-o {output} >& {log}"


# rule gatk_realign_bam:
#     input:
#         ref=_get_ref,
#         bam="mapping/{reference}/{prefix}.bam",
#         intervals="mapping/{reference}/{prefix}.realign.intervals"
#     output:
#         "mapping/{reference,[^/]+}/{prefix}.realigned.bam"
#     params:
#         custom=config.get("params_gatk", "")
#     log:
#         "mapping/log/{reference}/{prefix}.realign.log"
#     shell:
#         "gatk -T IndelRealigner -R {input.ref} {params.custom} "
#         "--disable_bam_indexing "
#         "-I {input.bam} -targetIntervals {input.intervals} "
#         "-o {output} >& {log}"


# # base recalibration
# rule gatk_recalibrate_info:
#     input:
#         "mapping/{reference}/{prefix}.bam.bai",
#         ref=_get_ref,
#         bam="mapping/{reference}/{prefix}.bam"
#     output:
#         temp("mapping/{reference,[^/]+}/{prefix}.recalibrate.grp")
#     params:
#         custom=config.get("params_gatk", "")
#     log:
#         "mapping/log/{reference}/{prefix}.recalibrate_info.log"
#     threads: 8
#     shell:
#         "gatk -T BaseRecalibrator -R {input.ref} {params.custom} "
#         "-nct {threads} "
#         "-I {input.bam} -knownSites {config[known_variants][dbsnp]} "
#         "-o {output} >& {log}"


# rule gatk_recalibrate_bam:
#     input:
#         ref=_get_ref,
#         bam="mapping/{reference}/{prefix}.bam",
#         grp="mapping/{reference}/{prefix}.recalibrate.grp"
#     output:
#         "mapping/{reference,[^/]+}/{prefix}.recalibrated.bam"
#     params:
#         custom=config.get("params_gatk", "")
#     log:
#         "mapping/log/{reference}/{prefix}.recalibrate.log"
#     threads: 8
#     shell:
#         "gatk -T PrintReads -R {input.ref} {params.custom} "
#         "-nct {threads} "
#         "--disable_bam_indexing "
#         "-I {input.bam} -BQSR {input.grp} "
#         "-o {output} >& {log}"


# variant calling
rule gatk_HaplotypeCaller:
    input:
        bam = rules.gatk_SplitNCigarReads.output.split,
        genome = GENOME
    output:
        vcf = join(OUTPUT_DIR, '{sample}', 'call.vcf')
    log:
        join(OUTPUT_DIR, '{sample}', 'gatk.haplotypecaller.log')
    run:
        shell('/opt/java/jdk1.7.0_71/bin/java -Xmx16g -jar /opt/gatk-3.4-46/gatk/GenomeAnalysisTK.jar'
              ' -T HaplotypeCaller'
              ' -R {input.genome}'
              ' -I {input.bam}'
              ' -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0'
              ' -o {output.vcf}'
              ' > {log} 2>&1')  

# Joint calling followed by VQSR
# untested; unclear if will cause issues for large number of cells...may want to filter first?
# https://www.broadinstitute.org/gatk/guide/article?id=3893
#rule gatk_CombineVariants:
#    input:
#        vcfs = expand(join(OUTPUT_DIR, '{sample}', 'call.vcf'), sample=SAMPLES),
#        genome = GENOME
#    output:
#        vcf = join(OUTPUT_DIR, 'comb.vcf')
#    log:
#        join(OUTPUT_DIR, 'gatk.combvar.log')
#    run:
#        shell('/opt/java/jdk1.7.0_71/bin/java -Xmx16g -jar /opt/gatk-3.4-46/gatk/GenomeAnalysisTK.jar' 
#              ' -T CombineVariants'
#              ' -R reference.fasta'
#   --variant input1.vcf \
#   --variant input2.vcf \
#   -o output.vcf \
#   -genotypeMergeOptions UNIQUIFY

# variant filtering 
rule gatk_VariantFiltration:
    input:
        vcf = rules.gatk_HaplotypeCaller.output.vcf,
        genome = GENOME
    output:
        vcf = join(OUTPUT_DIR, '{sample}', 'call.filtered.vcf')
    log:
        join(OUTPUT_DIR, '{sample}', 'gatk.varfilter.log')
    run:
        shell('/opt/java/jdk1.7.0_71/bin/java -Xmx16g -jar /opt/gatk-3.4-46/gatk/GenomeAnalysisTK.jar'
              ' -T VariantFiltration'
              ' -R {input.genome}'
              ' -V {input.vcf}'
              ' -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0"'
              ' -o {output.vcf}'
              ' > {log} 2>&1')  

