#!/usr/bin/env nextflow

params.run_dir = null
params.out_dir = null
params.samplesheet = null

helpMessage = """
    
    """.stripIntent()
interop_dir=$out_dir/InterOp

bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y151,I8,Y151 -R $run_dir --output-dir=$out_dir --interop-dir=$interop_dir --sample-sheet=$samplesheet -p 6 -r 6 -w 6
