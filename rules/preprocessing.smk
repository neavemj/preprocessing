"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned Illumina sequencing reads
"""

import os, sys

configfile: "../introduction/config.yaml"


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
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R1_U = config["sub_dirs"]["trim_dir"] + "/{sample}_1U.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz",
        R2_U = config["sub_dirs"]["trim_dir"] + "/{sample}_2U.fastq.gz"
    params:
        qual = config["trimmomatic_quality"],      
        adapters = config["program_dir"] + "preprocessing/adapters/" + config["trimmomatic_adapters"],
        minlen = config["trimmomatic_minlen"]
    threads: 8
    log:
        "logs/trimmomatic_PE/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["trim_dir"] + "/trimmomatic_PE/{sample}.txt"
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
        {config[program_dir]}/scripts/summarise_trimmomatic.py \
        -i {input} -o {output}
        """


rule all:
    input:
        expand(config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz", sample=config["samples"])
