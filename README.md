# AIV Pipeline Preprocessing Module
---

Module to take raw reads directly from the MiSeq and remove Illumina adapters and trim low-quality sequence. Could also remove any non-influenza reads (i.e., host) if that is better for the IRMA assembly module.

Input:
* Raw files from the MiSeq in fastq.gz format.

Output:
* Cleaned fastq.gz files.
