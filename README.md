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


### Installation

The following are instructions for installing all dependencies for the pipeline.

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
'''shell
conda install multiqc
'''

#### trim_galore
'''shell
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip

unzip 0.4.5.zip
'''

Add path to trim_galore executable to .bashrc


#### Picard

```shell
wget https://github.com/broadinstitute/picard/releases/download/2.18.7/picard.jar
```

#### Freebayes
```shell
conda install -c bioconda freebayes
'''

