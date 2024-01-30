process trim_primer {
    label = 'trim_primer'

    publishDir params.outdir, mode: 'copy'

    conda "${projectDir}/envs/biopython.yaml"

    input:
    path fastq
    val fwd_primer
    val rev_primer
    val mismatches
    val barcode_length
    val log

    output:
    path "*.fastq.gz"

    script:
    def fastqName = fastq.getSimpleName()
    def outFile = "${fastqName}_trimmed.fastq.gz"
    def logArg = log ? "" : "--log"
    """
    trim_primer.py \
        --fwd_primer ${fwd_primer} \
        --rev_primer ${rev_primer} \
        --mismatches ${mismatches} \
        --barcode_length ${barcode_length} \
        ${logArg} ${fastq} | gzip > ${outFile}
    """
}
