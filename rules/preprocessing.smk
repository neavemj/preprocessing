"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina or MinION sequencing reads
The outputs are:
    - cleaned Illumina or MinION sequencing reads
    
The technologies are differentiated by 
putting the program name in the file name
"""

import os, sys
from os.path import join

configfile: "config.yaml"

# 2023-09-26 updated to get sample names from directly from raw file directory
RAW_DIR = config["raw_dir"]

# need to differentiate between Illumina or MinION data
# will use the IRMA module config parameter for now
if config["irma_module"] == "FLU-avian":
    # if FLU-avain is used is must be Illumina data
    # assuming the following naming pattern
    # note: need the comma after SAMPLES, or a list is not created
    SAMPLES, = glob_wildcards(join(RAW_DIR, "{sample}_L001_R1_001.fastq.gz"))
    # if we print like below, the draw_dag rule throws an error
    # somehow this print statement gets put into the dag call..
    #print("Samples in this job are: ", SAMPLES)
elif config["irma_module"] == "FLU-minion":
        # assuming the fastq.gz file produced by the MinION will be cat'd together
        # won't work if they are uncompressed fastq files
    SAMPLES, = glob_wildcards(join(RAW_DIR, "{sample}.fastq.gz"))
    #print("Samples in this job are: ", SAMPLES)


# need a function here to differentiate MiSeq from MinION dataset
def get_reads(wildcards):
    if config["irma_module"] == "FLU-avian":
        # if FLU-avain is used is must be Illumina data
        # assuming the following naming pattern
        return([
            join(RAW_DIR, "{sample}_L001_R1_001.fastq.gz"),
            join(RAW_DIR, "{sample}_L001_R2_001.fastq.gz")
            ])
    elif config["irma_module"] == "FLU-minion":
        # assuming the fastq.gz file produced by the MinION will be cat'd together
        # won't work if they are uncompressed fastq files
        return([      
            join(RAW_DIR, "{sample}.fastq.gz")
            ])  


rule trim_all:
    input: 
        expand("01_preprocessing/{sample}_1P.fastq.gz", sample=SAMPLES)
        #expand("01_preprocessing/{sample}.porechop.nanofilt.fastq", sample=SAMPLES)


rule trimmomatic_PE:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        get_reads
    output:
        R1_P = "01_preprocessing/{sample}_1P.fastq.gz",
        R1_U = "01_preprocessing/{sample}_1U.fastq.gz",
        R2_P = "01_preprocessing/{sample}_2P.fastq.gz",
        R2_U = "01_preprocessing/{sample}_2U.fastq.gz"
    params:
        qual = config["trimmomatic_quality"],
        adapters = config["program_dir"] + "preprocessing/adapters/" + config["trimmomatic_adapters"],
        minlen = config["trimmomatic_minlen"]
    threads: 8
    log:
        "logs/trimmomatic_PE/{sample}.log"
    benchmark:
        "benchmarks/trimmomatic_PE/{sample}.txt"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """

rule summarise_trimmomatic_log:
    input:
        "logs/trimmomatic_PE/{sample}.log"
    output:
        "logs/trimmomatic_PE/trim_logs.summary"
    shell:
        """
        {config[program_dir]}/preprocessing/scripts/summarise_trimmomatic.py \
        -i {input} -o {output}
        """
        
rule porechop:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} MinION reads for adapters using Porechop
        """
    input:
        get_reads
    output:
        "01_preprocessing/{sample}.porechop.fastq",
    params:
 
    threads: 8
    log:
    benchmark:
        "benchmarks/porechop/{sample}.txt"
    shell:
        """
        porechop \
            -t {threads} \
            -i {input} \
            -o {output}
        """
        
rule nanofilt:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} MinION reads for quality using NanoFilt
        """
    input:
        "01_preprocessing/{sample}.porechop.fastq",
    output:
        "01_preprocessing/{sample}.porechop.nanofilt.fastq",
    params:
        nanofilt_quality = config["nanofilt_quality"]
    threads: 1
    log:
    benchmark:
        "benchmarks/nanofilt/{sample}.txt"
    shell:
        """
        NanoFilt \
            --quality {params.nanofilt_quality} \
            {input} \
            > {output}
        """
