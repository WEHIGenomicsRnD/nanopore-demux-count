process trim_primer {
    label = 'trim_primer'

    publishDir params.outdir, mode: 'copy'

    conda "${projectDir}/envs/cutadapt.yaml"

    input:
    path fastq
    val primer
    val mismatches
    val primer_type

    output:
    path "*.fastq"

    script:
    def fastqName = fastq.getSimpleName()
    def option = primer_type == 'fwd' ? '--front' : '--adapter'
    def outFile = "${fastqName}_${primer_type}.fastq"
    """
    cutadapt \
        ${option} ${primer} \
        --errors ${mismatches} \
        --cores ${task.cpus} \
        --discard-untrimmed \
        -o ${outFile} \
        ${fastq}
    """
}
