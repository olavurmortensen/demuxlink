#!/usr/bin/env nextflow

params.run_dir = null
params.out_dir = null
params.samplesheet = null
params.help = false

helpMessage = """
    
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.run_dir != null, 'Input parameter "run_dir" cannot be unassigned.'
assert params.out_dir != null, 'Input parameter "out_dir" cannot be unassigned.'
assert params.samplesheet != null, 'Input parameter "samplesheet" cannot be unassigned.'

run_dir = file(params.run_dir)
out_dir = file(params.out_dir)
samplesheet = file(params.samplesheet)
interop_dir = file(out_dir + "InterOp")

println "D E M U X    L I N K    "
println "================================="
println "run_dir             : ${run_dir}"
println "out_dir             : ${out_dir}"
println "samplesheet         : ${samplesheet}"
println "interop_dir         : ${interop_dir}"

// Call bcl2fastq, performing simultaneous basecalling and demultiplexing.
// TODO:
// Do I need any of the extra parameters in longranger's call to bcl2fastq?
// reading, writing and processing threads (p, r and w) parameter.
process bcl2fastq {
    output:
    file "out/*fastq.gz" into fastq_ch1, fastq_ch2

    script:
    """
    mkdir out
    bcl2fastq \
        -R $run_dir \
        -o out \
        --interop-dir $interop_dir \
        --sample-sheet $samplesheet \
        -p 6 -r 6 -w 6
    # longranger bcl2fastq call is:
    #bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y151,I8,Y151 -R $run_dir --output-dir=$out_dir --interop-dir=$interop_dir --sample-sheet=$samplesheet -p 6 -r 6 -w 6
    """
}

// Get sample names from FASTQ filenames.
samplenames_ch = fastq_ch1.map { it.toString() }
    .map { it.replaceAll(/\/.*\//, "") }  // Remove everything leading up to the sample name.
    .map { it.replaceAll(/_.*/, "") }  // Remove everything after the sample name.

// Get (sample name, FASTQ paths) tuples.
// This assumes samplenames_ch and fastq_ch2 are in the same order.
fastq_ch = samplenames_ch.merge(fastq_ch2)
    .groupTuple()

// Merge FASTQs by indexes and lanes.
process merge {

    publishDir "$outdir/fastq_out/$sample", mode: 'copy'

    input:
    set sample, file(fastqs) from fastq_ch

    output:
    file "$sample*fastq.gz" into merged_ch

    script:
    """
    zcat $sample\\_*R1*.fastq.gz | gzip -c > $sample\\_R1.fastq.gz
    zcat $sample\\_*R2*.fastq.gz | gzip -c > $sample\\_R2.fastq.gz
    """
}






