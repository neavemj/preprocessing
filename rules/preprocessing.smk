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

configfile: "config.yaml"

rule trim_all:
    input: expand("01_preprocessing/{sample}.porechop.nanofilt.fastq", sample=config["samples"])


# function to get fastq file locations from the config file
def getFastq(wildcards):
    return config['samples'][wildcards.sample]


# TODO: might need a trimmomatic SE mode
rule trimmomatic_PE:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        reads = getFastq
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
            {input.reads} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """

rule summarise_trimmomatic_log:
    input:
        expand("logs/trimmomatic_PE/{sample}.log", sample=config["samples"])
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
        reads = getFastq
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
