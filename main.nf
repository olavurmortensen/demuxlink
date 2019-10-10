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
    //echo true 

    output:
    file "outs/*/*fastq.gz" into fastq_samplenames_ch

    script:
    """
    bcl2fastq \
        -R $rundir \
        -o outs \
        --interop-dir interop \
        --sample-sheet $samplesheet \
        -p 6 -r 6 -w 6
    # longranger bcl2fastq call is:
    #bcl2fastq --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --create-fastq-for-index-reads --ignore-missing-positions --ignore-missing-filter --ignore-missing-bcls --use-bases-mask=Y151,I8,Y151 -R $rundir --output-dir=$outdir --interop-dir=$interop_dir --sample-sheet=$samplesheet -p 6 -r 6 -w 6
    """
}

// Get (sample names, FASTQ files) tuples.
// First flatten the channel because each instance of process "bcl2fastq" outputs a tuple.
// Then map the channel to (sample name, FASTQ path) tuples. This channel has one record per file.
// Then group by sample name to get (sample name, FASTQ list) tuples.
fastq_ch = fastq_samplenames_ch.flatten()
    .map { file ->
        def key = file.name.toString().split('_')[0]
        return tuple(key, file)}
    .groupTuple()

// Merge FASTQs by indexes and lanes.
process merge {
    publishDir "$outdir/fastq_out/$sample", mode: 'copy'

    input:
    set sample, file(fastqs) from fastq_ch

    output:
    //file "outs/$sample*fastq.gz" into fastq_out_ch
    set sample, file("outs/$sample*fastq.gz") into fastq_qc_ch

    script:
    """
    #touch $sample\\_blabla.fastq.gz
    # One merging operation for each lane/read pair.
    # If there are not multple indexes for each sample, then the "merged" file will have the
    # same name as the original. Therefore, we do the merging operation in another folder.
    mkdir outs
    zcat $sample\\_*L005\\_R1*.fastq.gz | gzip -c > outs/$sample\\_L005\\_R1.fastq.gz
    zcat $sample\\_*L005\\_R2*.fastq.gz | gzip -c > outs/$sample\\_L005\\_R2.fastq.gz
    #zcat $sample\\_*L001\\_R1*.fastq.gz | gzip -c > $sample\\_L001\\_R1.fastq.gz
    #zcat $sample\\_*L002\\_R1*.fastq.gz | gzip -c > $sample\\_L002\\_R1.fastq.gz
    #zcat $sample\\_*L003\\_R1*.fastq.gz | gzip -c > $sample\\_L003\\_R1.fastq.gz
    #zcat $sample\\_*L004\\_R1*.fastq.gz | gzip -c > $sample\\_L004\\_R1.fastq.gz
    #zcat $sample\\_*L001\\_R2*.fastq.gz | gzip -c > $sample\\_L001\\_R2.fastq.gz
    #zcat $sample\\_*L002\\_R2*.fastq.gz | gzip -c > $sample\\_L002\\_R2.fastq.gz
    #zcat $sample\\_*L003\\_R2*.fastq.gz | gzip -c > $sample\\_L003\\_R2.fastq.gz
    #zcat $sample\\_*L004\\_R2*.fastq.gz | gzip -c > $sample\\_L004\\_R2.fastq.gz
    """
}

// Run FastQC for QC metrics of raw data.
// Note that FastQC will allocate 250 MB of memory per thread used.
// TODO: how many threads are needed?
process fastqc_analysis {
    memory = "250MB"
    cpus = 1

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
    unset DISPLAY
    mkdir tmp
    fastqc -q --dir tmp --outdir . $fastq_list
    """
}

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



