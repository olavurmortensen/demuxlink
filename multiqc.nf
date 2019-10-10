#!/usr/bin/env nextflow

params.outdir = null
params.help = false

helpMessage = """
    
    """.stripIndent()

if (params.help) {
    log.info helpMessage
    exit 0
}

assert params.outdir != null, 'Input parameter "outdir" cannot be unassigned.'

rundir = file(params.rundir)
outdir = file(params.outdir)
samplesheet = file(params.samplesheet)
interop_dir = file(outdir + "InterOp")

println "M u l t i Q C    D E M U X    L I N K    "
println "================================="
println "outdir              : ${outdir}"

process multiqc {
    publishDir "$outdir/multiqc", mode: 'copy', overwrite: true

    input:
    val status from status_ch

    output:
    file "multiqc_report.html" into multiqc_report_ch
    file "multiqc_data" into multiqc_data_ch

    script:
    """
    multiqc -f $outdir --config ${params.multiqc_config}
    """
}



