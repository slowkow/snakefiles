#!/usr/bin/env bash
# download_files.sh
#
# Download files from Ensembl.

BASE=ftp://ftp.ensembl.org/pub/release-82

URL_DNA=$BASE/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
URL_GTF=$BASE/gtf/homo_sapiens/Homo_sapiens.GRCh38.82.gtf.gz
URL_CDNA=$BASE/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

curl -k -L $URL_CDNA | gzip -dc > Homo_sapiens.GRCh38.cdna.all.fa

curl -k -L $URL_DNA | gzip -dc > Homo_sapiens.GRCh38.dna.primary_assembly.fa

curl -k -L $URL_GTF | gzip -dc > Homo_sapiens.GRCh38.82.gtf

