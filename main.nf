#!/usr/bin/env nextflow

params.rundir = null
params.outdir = null
params.samplesheet = null
params.threads = null
params.mem = null
params.help = false

helpMessage = """
    Input:
    rundir:         Path to FASTQ run directory.
    outdir:         Where to store output data.
    samplesheet:    Path to sample sheet.
    threads:        Number of threads to use.
    mem:            Amount of memory to use (e.g. "10 GB").
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.rundir != null, 'Input parameter "rundir" cannot be unassigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unassigned.'
assert params.samplesheet != null, 'Input parameter "samplesheet" cannot be unassigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unassigned.'

if(params.threads == null) {
    params.threads = 1
}

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
// As --use-bases-mask is not specified, RunInfo.xml will be used to determine the base masking.
// Adapter sequences (read 1 and read2) should be contained in the sample sheet.
// TODO:
// reading, writing and processing threads (p, r and w) parameter.
// Memory and cores?
// Do I need any of the extra parameters in longranger's call to bcl2fastq?
process bcl2fastq {
    output:
    file "outs/*/*fastq.gz" into fastq_samplenames_ch

    script:
    if(params.threads > 20) {
        p_threads = params.threads - 8
        w_threads = 4
        r_threads = 4
    } else if(params.threads > 10) {
        p_threads = params.threads - 2
        w_threads = 1
        r_threads = 1
    } else if(params.threads > 3) {
        p_threads = params.threads - 2
        w_threads = 1
        r_threads = 1
    } else {
        p_threads = 1
        w_threads = 1
        r_threads = 1
    }
    """
    bcl2fastq \
        -R $rundir \
        -o outs \
        --interop-dir interop \
        --sample-sheet $samplesheet \
        -p $p_threads -r $r_threads -w $w_threads
    # longranger bcl2fastq call is:
    #bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y151,I8,Y151 -R $rundir --output-dir=$outdir --interop-dir=$interop_dir --sample-sheet=$samplesheet -p 6 -r 6 -w 6
    """
}

// Get (key, FASTQ files) tuples, where key is a (sample names, lane, read) tuple.
// First flatten the channel because each instance of process "bcl2fastq" outputs a tuple.
// Then map the channel to (key, FASTQ path) tuples. This channel has one record per file.
// Then group by sample name to get (key, FASTQ list) tuples.
fastq_ch = fastq_samplenames_ch.flatten()
    .map { file ->
        def sample = file.name.toString().split('_')[0]
        def lane = file.name.toString().split('_')[2]
        def read = file.name.toString().split('_')[3]
        def key = tuple(sample, lane, read)
        return tuple(key, file)}
    .groupTuple()

// Since 10x samples have multiple indexes per sample, we merge these.
// Since this process is only concatenating (cat) and zipping files, it doesn't need much memory or many cores.
process merge {
    memory = "250MB"
    cpus = 1

    publishDir "$outdir/fastq_out/$sample", mode: 'copy'

    input:
    set key, file(fastqs) from fastq_ch

    output:
    set sample, file("$sample*fastq.gz") into fastq_qc_ch

    script:
    sample = key[0]
    lane = key[1]
    read = key[2]
    """
    zcat $sample\\_*$lane\\_$read*.fastq.gz | gzip -c > $sample\\_$lane\\_$read\\_merged.fastq.gz
    """
}

// The merged FASTQ files are grouped by sample, lane and read (but the tuple only contains sample
// and file path). Group tuple by key so that we have all files in the sample in one channel element.
fastq_qc_ch = fastq_qc_ch.groupTuple()

// FastQC allocates 250 MB of memory per thread used. Therefore, we need to calculate
// how many threads we can "afford", given how much memory we have available.
// The maximum number of threads we can use.
fastqc_threads = Math.floor((params.mem as nextflow.util.MemoryUnit).toMega() / 250) as Integer

// If the calculated number of threads is greater than the available amount, we default to that.
fastqc_threads =  Math.min(fastqc_threads, params.threads)

// The number of bytes of memory that corresponds to.
mem_b = 250 * 1048576 * (fastqc_threads as Long)

// Memory in megabytes.
fastqc_mem = (mem_b as nextflow.util.MemoryUnit).toMega()

// Run FastQC for QC metrics of raw data.
process fastqc_analysis {
    // Use the memory and threads calculated above.
    memory = fastqc_mem
    cpus = fastqc_threads

    publishDir "$outdir/fastqc/$sample", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set sample, file(fastqs) from fastq_qc_ch

    output:
    set sample, file('*.{zip,html}') into fastqc_report_ch
    set sample, file('.command.out') into fastqc_stdout_ch

    script:
    fastq_list = (fastqs as List).join(' ')
    """
    # We unset the DISPLAY variable to avoid having FastQC try to open the GUI.
    unset DISPLAY
    mkdir tmp
    fastqc -q --dir tmp --threads $fastqc_threads --outdir . $fastq_list
    """
}


