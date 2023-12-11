process BAM_READCOUNT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bam-readcount=0.8"
    container "docker.io/mgibio/bam-readcount"

    input:
    tuple val(meta), path(bam), path(bai)
    path(sites)
    path(fasta)

    output:
    path "*bamreadcount.txt", emit: readcounts
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filename = meta.id + ".bamreadcount.txt"
    """
    bam-readcount -f $fasta \\
        -l $sites \\
        -d 500000000 -w 100 \\
        $bam \\
        > $filename

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam-readcount: \$(bam-readcount -v)
    END_VERSIONS
    """
}
