'''
Armon Azizi

Snakemake pipeline for variant calling analysis

Place all paired fastq files in directory named 'fastq'
Run command while inside of parent directory of fastq

To run, execute command:
    snakemake -s variant.py -j 10 --cluster "qsub -q {params.queue} -l {params.time} -l {params.nodes}"

'''
import sys
from os.path import join



'''
Function to write messages to command line

Syntax is same as python
'''
def printcmd(cmd):
    sys.stderr.write(cmd)
    

'''
Path to directory containing reference genomes for
bowtie2 and mutect (or other variant caller)
'''
HG19_BT = "~/bejar_variant_dependencies/reference_genomes/bowtie/hg19"
HG19_RAW = "~/bejar_variant_dependencies/reference_genomes/raw/hg19.ucsc.fasta"

#Path to picard
PICARD = "~/bejar_variant_dependencies/Picard/picard.jar"

#Path to trim_galore
TRIM_GALORE = "~/bejar_variant_dependencies/trim_galore/TrimGalore-0.4.5/trim_galore"

#Path to vcf file containing variants to filter out
FILTER = "~/bejar_variant_dependencies/vcfs/filter_sites.vcf.gz"

'''
Obtain working set of fastq files 
-use glob command to find files in directory called 'fastq'
-sort filenames for consistency
-print filenames to command line
'''
FILES, = glob_wildcards("fastq/{id}.read1.gz")
sorted(FILES)

printcmd("Will process files with ids:" + '\n')

for f in FILES:
    printcmd('\t' + f + '\n')

printcmd('\n')




'''
Define files composing final output

The final output defines the dependencies for the rest of the rules
'''
rule final_output:
    input:
        ['analysis/multiqc_output/multiqc_report.html']+
        expand('analysis/filtered_variants/{sample}.freebayes.filtered.vcf', sample=FILES)

'''
QC analysis

First, fastqc all files, then compile output with multiqc

MultiQC rule is dependent on fastqc rule output

NOTE: Each rule has a section called params that defines the 
parameters used for submitting that rule as a job to the TSCC
params also includes any variables for the shell command
'''
rule multiqc:
    input:
        expand('analysis/fastqc_reports/{sample}.read{num}_fastqc.zip', sample=FILES, num=[1,2])
    output:
        'analysis/multiqc_output/multiqc_report.html'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Compiling fastqc reports with multiqc..."
    shell:
        'multiqc -o analysis/multiqc_output analysis/fastqc_reports'
        
rule fastqc:
    input:
        'fastq/{sample}.read1.gz',
        'fastq/{sample}.read2.gz'
    output:
        'analysis/fastqc_reports/{sample}.read1_fastqc.zip',
        'analysis/fastqc_reports/{sample}.read2_fastqc.zip'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Checking {input} quality with fastqc..."
    shell:
        "module load fastqc; "
        'fastqc -o analysis/fastqc_reports/ {input}'
        


'''
read trimming
'''
rule trim_reads:
    input:
        'fastq/{sample}.read1.gz',
        'fastq/{sample}.read2.gz'
    output:
        'analysis/trimmed_reads/{sample}.read1.gz_val_1.fq.gz',
        'analysis/trimmed_reads/{sample}.read2.gz_val_2.fq.gz'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Trimming {input} adapters"
    shell:
        '{TRIM_GALORE} --paired -o analysis/trimmed_reads {input}'


'''
read mapping
'''
rule map_reads:
    input:
        mate1='analysis/trimmed_reads/{sample}.read1.gz_val_1.fq.gz',
        mate2='analysis/trimmed_reads/{sample}.read2.gz_val_2.fq.gz'
    output:
        'analysis/sam_files/{sample}.sam'
    params:
        ref=HG19_BT,
        queue = "hotel",
        nodes = "nodes=1:ppn=4",
        time = "walltime=01:00:00",
        threads = "4"
    message:
        "Mapping {input} using bowtie2"
    shell:
        "module load bowtie2; "
        "bowtie2 -x {params.ref} -p {params.threads}  -1 {input.mate1} -2 {input.mate2} -S {output}"

'''
Convert .sam to .bam
'''
rule sam_to_bam:
    input:
        'analysis/sam_files/{sample}.sam'
    output:
        'analysis/bam_files/{sample}.bam'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Converting .sam to .bam"
    shell:
        "module load samtools; "
        "samtools view -Sb {input} > {output}"
        


'''
Sort and index .bam
'''
rule sort_bam:
    input:
        'analysis/bam_files/{sample}.bam'
    output:
        'analysis/bam_files/{sample}_sorted.bam'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Sorting .bam"
    shell:
        "module load samtools; "
        "samtools sort {input} > {output}; "
        "samtools index {output}"
        
'''
Remove Duplicates with Picard MarkDuplicates
'''
rule markdups:
    input:
        'analysis/bam_files/{sample}_sorted.bam'
    output:
        bam='analysis/bam_files/{sample}_nodup.bam',
        metrics='analysis/duplicate_removal_reports/{sample}_rmdup_metrics.txt'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Finding and Removing PCR Duplicates"
    shell:
        "java -jar {PICARD} "
        "MarkDuplicates I={input} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true"


'''
Call Variants with varscan

Minimum of 5 reads and >2% allele frequency
'''
'''
rule call_variants:
    input:
        'analysis/bam_files/{sample}_nodup.bam'
    output:
        'analysis/variants/{sample}.vcf'
    params:
        ref=HG19_RAW,
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=01:00:00"
    message:
        "Calling variants on {input}"
    shell:
        "module load samtools; "
        "samtools mpileup -B -f {params.ref} {input} | "
        "java -jar {VARSCAN} "
        "mpileup2cns --min-coverage 5 --min-var-freq 0.02 "
        "--output-vcf --variants > {output}" 
'''
        
'''
Call Variants with freebayes

Does local indel realignment automatically

Minimum of 5 reads and >2% allele frequency
'''
rule call_variants:
    input:
        'analysis/bam_files/{sample}_nodup.bam'
    output:
        'analysis/raw_variants/{sample}.freebayes.vcf'
    params:
        ref=HG19_RAW,
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=02:00:00"
    message:
        "Calling variants on {input}"
    shell:
        "freebayes -f {params.ref} -F 0.02 -C 5 -v {output} {input}"


'''
Filter variants by removing variants found in FILTER file.
'''
rule filter:
    input:
        'analysis/raw_variants/{sample}.freebayes.vcf'
    output:
        A = 'analysis/filtered_variants/{sample}.freebayes.filtered.vcf.gz',
        B = 'analysis/temp/{sample}'
    params:
        queue = "hotel",
        nodes = "nodes=1:ppn=1",
        time = "walltime=02:00:00",
    message:
        "Filtering variants out of {input}"
    shell:
        "module load bcftools; "
        "bgzip {input}; "
        "tabix {input}.gz; "
        "bcftools isec  -O z -p {output.B} -c all -w1 {input}.gz {FILTER}; "
        "mv {output.B}/0000.vcf.gz {output.A}; "
        "rm -r {output.B}/*"
        





