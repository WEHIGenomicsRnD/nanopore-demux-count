process TrimPrimer {
    label = 'TrimPrimer'

    publishDir params.outdir, mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    path fastq
    val fwd_primer
    val rev_primer
    val mismatches
    val barcode_length
    val output_untrimmed

    output:
    path "*_trimmed.fastq*", emit: trimmed_ch
    path "*_untrimmed.fastq*", emit: untrimmed_ch, optional: true
    path "*.txt"

    script:
    def fastqName = fastq.getSimpleName()
    def outFile = "${fastqName}_trimmed.fastq.gz"
    def untrimmedFastq = output_untrimmed ? "--untrimmed_fastq ${fastqName}_untrimmed.fastq.gz" : ""
    """
    trim_primer.py \
        --fwd_primer ${fwd_primer} \
        --rev_primer ${rev_primer} \
        --mismatches ${mismatches} \
        --barcode_length ${barcode_length} \
        ${untrimmedFastq} \
        ${fastq} | gzip > ${outFile}
    """
}
