#!/usr/bin/env nextflow

params.rundir = null
params.outdir = null
params.samplesheet = null
params.help = false

helpMessage = """
    
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.rundir != null, 'Input parameter "rundir" cannot be unassigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unassigned.'
assert params.samplesheet != null, 'Input parameter "samplesheet" cannot be unassigned.'

rundir = file(params.rundir)
outdir = file(params.outdir)
samplesheet = file(params.samplesheet)
interop_dir = file(outdir + "InterOp")

println "D E M U X    L I N K    "
println "================================="
println "rundir              : ${rundir}"
println "outdir              : ${outdir}"
println "samplesheet         : ${samplesheet}"
println "interop_dir         : ${interop_dir}"

// Call bcl2fastq, performing simultaneous basecalling and demultiplexing.
// TODO:
// Do I need any of the extra parameters in longranger's call to bcl2fastq?
// reading, writing and processing threads (p, r and w) parameter.
process bcl2fastq {
    output:
    file "out/*fastq.gz" into fastq_samplenames_ch, fastq_filegroups_ch, fastq_fastqc_ch

    script:
    """
    mkdir out
    bcl2fastq \
        -R $rundir \
        -o out \
        --interop-dir $interop_dir \
        --sample-sheet $samplesheet \
        -p 6 -r 6 -w 6
    # longranger bcl2fastq call is:
    #bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y151,I8,Y151 -R $rundir --output-dir=$outdir --interop-dir=$interop_dir --sample-sheet=$samplesheet -p 6 -r 6 -w 6
    """
}

// Get sample names from FASTQ filenames.
samplenames_ch = fastq_samplenames_ch.map { it.toString() }
    .map { it.replaceAll(/\/.*\//, "") }  // Remove everything leading up to the sample name.
    .map { it.replaceAll(/_.*/, "") }  // Remove everything after the sample name.

// Get (sample name, FASTQ paths) tuples.
// This assumes samplenames_ch and fastq_ch2 are in the same order.
fastq_ch = samplenames_ch.merge(fastq_filegroups_ch)
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


// Run FastQC for QC metrics of raw data.
// Note that FastQC will allocate 250 MB of memory per thread used.
// TODO: how many threads are needed?
//process fastqc_analysis {
//    memory = "250MB"
//    cpus = 1
//
//    publishDir "$outdir/fastqc", mode: 'copy',
//        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
//
//    input:
//    set fastq from fastq_fastqc_ch
//
//    output:
//    set sample, file('*.{zip,html}') into fastqc_report_ch
//    set sample, file('.command.out') into fastqc_stdout_ch
//
//    script:
//    """
//    unset DISPLAY
//    mkdir tmp
//    fastqc -q --dir tmp --outdir . $fastq
//    """
//}
//
//process multiqc {
//    publishDir "$outdir/multiqc", mode: 'copy', overwrite: true
//
//    input:
//    val status from status_ch
//
//    output:
//    file "multiqc_report.html" into multiqc_report_ch
//    file "multiqc_data" into multiqc_data_ch
//
//    script:
//    """
//    multiqc -f $outdir --config ${params.multiqc_config}
//    """
//}



