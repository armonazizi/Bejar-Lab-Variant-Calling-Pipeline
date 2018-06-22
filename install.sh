#!/usr/bin/env bash

#install snakemake
conda install -c bioconda -c conda-forge snakemake

# install multiqc
conda install multiqc

#install trim_galore
conda install -c bioconda cutadapt

wget -P ~/bejar_variant_dependencies/trim_galore https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip

unzip -d ~/bejar_variant_dependencies/trim_galore/ ~/bejar_variant_dependencies/trim_galore/0.4.5.zip

#install Picard
wget -P ~/bejar_variant_dependencies/Picard https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar

#install Freebayes
conda install -c bioconda freebayes

#get bowtie index
wget -P ~/bejar_variant_dependencies/reference_genomes/bowtie ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip

unzip -d ~/bejar_variant_dependencies/reference_genomes/bowtie/ ~/bejar_variant_dependencies/reference_genomes/bowtie/hg19.zip

#get raw reference
wget -P ~/bejar_variant_dependencies/reference_genomes/raw http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

chmod 775 ~/bejar_variant_dependencies/reference_genomes/raw/twoBitToFa

wget -P ~/bejar_variant_dependencies/reference_genomes/raw http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

~/bejar_variant_dependencies/reference_genomes/raw/twoBitToFa ~/bejar_variant_dependencies/reference_genomes/raw/hg19.2bit ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta

#index reference genome
module load bcftools

bgzip ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta
tabix -p vcf ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta.gz



