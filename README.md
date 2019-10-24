# Demuxlink [WIP]

This Nextflow pipeline basecalls and demultiplexes linked-reads from 10x Genomics. To run this pipeline, the 8-base sample indexes are needed, corresponding to the 10x Genomics indexes (e.g. `SI-GA-A1`).

This pipeline makes some assumptions about the input data. For example, it makes the assumption that it is paired-end sequencing, and therefore uses `--use-bases-mask=Y*,I*,Y*` in `bcl2fastq`, and assumes that the read lengths (and index length) is found in `RunInfo.xml`.

# Running on tiny-bcl

Here's how to run this pipeline on the "tiny-bcl" example dataset from 10x Genomics. First of all, download the tiny-bcl tar file and the Illumina Experiment Manager sample sheet: tiny-bcl-samplesheet-2.1.0.csv.

> Download the tiny-bcl data:
> https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/mkfastq#example_data

Next, edit the `[Data]` part of the samplesheet from the following:

```
Lane,Sample_ID,index,Sample_Project
5,Sample1,SI-GA-C5,tiny_bcl
```

To the following:

```
Lane,Sample_ID,index
5,Sample1,CGACTTGA
5,Sample1,TACAGACT
5,Sample1,ATTGCGTG
5,Sample1,GCGTACAC
```

Using that the index `SI-GA-C5` corresponds to the four octamers `CGACTTGA,TACAGACT,ATTGCGTG,GCGTACAC`. Notice that we also removed the `Sample_Project` column.

Provided that you've set installed all the software in `environment.yml` (or maybe used this pipeline's Docker container at https://hub.docker.com/r/olavurmortensen/demuxlink), you should be able to run the pipeline like this:

```
nextflow demuxlink/main.nf --rundir tiny-bcl-2.2.0 --outdir results --samplesheet tiny-bcl-samplesheet-2.1.0.csv
```

And get the FASTQ files in `results/fastq_out/Sample1/outs` and the FastQC reports in `results/fastqc/Sample1`.

