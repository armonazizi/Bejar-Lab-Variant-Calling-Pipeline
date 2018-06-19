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

Follow the instructions below to install each tool separately or to install using custom parameters. Alternatively, skip to scipt installation for instructions on how to run the install script.

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

wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip

unzip 0.4.5.zip
```

Add path to trim_galore executable to .bashrc


#### Picard

```shell
wget https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar
```

#### Freebayes
```shell
conda install -c bioconda freebayes
```

### Script Installation

To install using a script, conda with python 3 must already be installed on the system. Then, run install.sh. This script will install the dependencies into the users home directory in the bejar_variant_dependencies directory.

```shell
./install.sh
```


### Running The Pipeline
Running the pipeline is relatively simple. The script command is a one-liner and will submit all jobs to the TSCC using PBS scripts.

First, files must be appropriately named. The current naming scheme is for each paired end sample to have one name followed by paired identifiers.
* {NAME}.read1.gz
* {NAME}.read2.gz

Alternitively, The main script can be modified to accept differently named files.

All fastq files should be in a directory titled fastq, and the snakemake script should be run in the 'fastq' parent directory.

Command to run pipeline:
```shell
snakemake -s variant.py -j <Maximum Number of jobs to submit concurrently> --cluster "qsub -q {params.queue} -l {params.time} -l {params.nodes}"
```

