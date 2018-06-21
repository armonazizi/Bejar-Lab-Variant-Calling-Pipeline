# Bejar-Lab-Variant-Calling-Pipeline
Pipeline designed to be run on the UC San Diego Triton Shared Computing Cluster to call somatic variants on unmatched tumor targeted exome sequencing samples.


### Dependencies

The pipeline requires the following dependencies to run:
* Python 3
* Snakemake (Installed using python 3)
* multiqc
* trim_galore
* Picard
* freebayes

The pipeline utilizes the following tools that are already installed on the TSCC
* fastqc
* bowtie2
* samtools


### Manual Installation

The following are instructions for installing all dependencies for the pipeline.

Follow the instructions below to install each tool separately or to install using custom parameters. Alternatively, skip to script installation for instructions on how to run the install script.

#### Snakemake

**The easiest way to install snakemake is to use conda.**

If conda is not already installed in the environment, follow instructions at: https://conda.io/docs/user-guide/install/linux.html to install miniconda.

Make sure that a python3 version of conda is running by running the following command.

```shell
python --version
```

If conda is not running a python 3 version, run the following command to update it.

```shell
conda install python=3.6
```

Then, install snakemake:
```shell
conda install -c bioconda -c conda-forge snakemake
```

#### multiqc
```shell
conda install multiqc
```

#### trim_galore
```shell
conda install -c bioconda cutadapt

wget -P ~/bejar_variant_dependencies/trim_galore https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip

unzip ~/bejar_variant_dependencies/trim_galore/0.4.5.zip
```


#### Picard

```shell
wget -P ~/bejar_variant_dependencies/Picard https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar
```

#### Freebayes
```shell
conda install -c bioconda freebayes
```

#### Bowtie index
```shell
wget -P ~/bejar_variant_dependencies/reference_genomes/bowtie ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip

unzip ~/bejar_variant_dependencies/reference_genomes/bowtie/hg19.ebwt.zip
```

#### UCSC Index
To get the raw hg19 fasta file, the compressed version must be downloaded, unzipped, and indexed.

```shell
wget -P ~/bejar_variant_dependencies/reference_genomes/raw http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

chmod 775 ~/bejar_variant_dependencies/reference_genomes/raw/twoBitToFa

wget -P ~/bejar_variant_dependencies/reference_genomes/raw http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

./~/bejar_variant_dependencies/reference_genomes/raw/twoBitToFa ~/bejar_variant_dependencies/reference_genomes/raw/hg19.2bit ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta

module load bcftools

bgzip ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta
tabix -p vcf ~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta.gz
```

#### vcf Filter Sites

In order to filter out unwanted variants, the pipeline takes one vcf file containing all unwanted sites and removes these variants from the final output, or a 'blacklist' of sites. In order to obtain this file, the user must manually generate it to the wanted specifications. Alternitively, there is a filter file containing common ExAC sites with the exception of common mutation sites found in cosmic that can be found in the bejar lab files.

After obtaining the vcf file, it's path should be the following:
```shell
~/bejar_variant_dependencies/vcfs/filter_sites.vcf
```

Note: if downloading the vcf file from ExAC or cosmic, the chromosomes must be renamed to match the freebayes vcf output. For example, in the chromosome column, '1' should be renamed to 'chr1' etc.

To do this, an awk one-liner can be used. The following script takes as input a compressed vcf file and appends 'chr' to the relevant lines, ignoring headers. It then indexes the file with tabix. (replace input and output with relevant filenames)

```shell
module load bcftools

zcat <Input.vcf.gz> | awk -v OFS='\t' '/^#/{print $0;} !/^#/{$1="chr"$1; print}' | bgzip > <output.vcf.gz>

tabix -p vcf <output from previous command>
```

Note: it's a good idea to do this step in the oasis filesystem because the vcf files can get very large. Oasis is optimized to handle large files, but your home directory can run out of space easily.

### Script Installation

To install using a script, conda with python 3 must already be installed on the system. Then, run install.sh. This script will install the dependencies into the users home directory in the bejar_variant_dependencies directory.

```shell
./install.sh
```

Note: Script installation will not download the blacklist vcf file. It must be downloaded separately before running the pipeline.


### Running The Pipeline
Running the pipeline is relatively simple. The script command is a one-liner and will submit all jobs to the TSCC using PBS scripts.

First, files must be appropriately named. The current naming scheme is for each paired end sample to have one name followed by paired identifiers.
* {NAME}.read1.gz
* {NAME}.read2.gz

Alternitively, The main script can be modified to accept differently named files.

All fastq files should be in a directory titled fastq, and the snakemake script should be run in the 'fastq' parent directory.

The '-n' flag can be added to the following command to do a "dry run" of the pipeline. This will output all steps that will be taken without actually executing them.

Command to run pipeline:
```shell
snakemake -s variant.py -j <Maximum Number of jobs to submit concurrently> --cluster "qsub -q {params.queue} -l {params.time} -l {params.nodes}"
```

If the pipeline crashes or times out midway through the analysis, it can be started again without redoing the analysis that was already completed. In order to do this, delete all files that are incomplete, or that were being generated when the pipeline crashed. Then, run the command that was used to start the pipeline.

