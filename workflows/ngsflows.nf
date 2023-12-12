
nextflow.enable.dsl = 2

// CONFIG FILES
ch_multiqc_config          = Channel.fromPath("assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// nf-core modules
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC } from '../modules/nf-core/fastqc'
include { MULTIQC } from '../modules/nf-core/multiqc/main'


def read_manifest(manifest, check_bams=true) {
  Channel
    .fromPath( manifest )
    .ifEmpty( "Manifest empty." )
    .splitCsv( header: true)
    .map { row -> if (row.sample_id =~ /^[^#]/) { row } }             // remove records starting with #
    .map { row ->
      [
        "meta": ["id": row.sample_id, "patient_id": row.patient_id, "sample_id": row.sample_id],
        "bam": file(row.bam, checkIfExists: true),
        "bai": file(row.bai, checkIfExists: true)
      ]
    }
}


workflow NGSFLOWS {
    ch = read_manifest(params.input)
    ch_versions = Channel.empty()

    // MODULE: Run FastQC
    ch_fastqc = ch.map { row ->
      [
        "meta": row.meta,
        "reads": file(row.bam, checkIfExists: true)
      ]
    }

    FASTQC ( ch_fastqc )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // MODULE: dump suftware versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // // MODULE: MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}
