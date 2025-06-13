process TrimPrimer {
    label = 'TrimPrimer'

    publishDir params.outdir, mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(sampleName), path(fastq)
    val fwd_primer
    val rev_primer
    val mismatches
    val barcode_length
    val output_untrimmed

    output:
    tuple val(sampleName), path("*_trimmed.fastq*"), emit: trimmed_ch
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


process MergeFastq{
    label = 'MergeFastq'

    publishDir "${params.outdir}/merge_fastq", mode: 'copy'

    input:
    tuple val(key),val(fastqs)

    output:
    path "*.fastq.gz", emit: merge_ch

    script:
    def merged = "${key}.fastq.gz"
    """
       cat ${fastqs.join(' ')} > ${merged}
    """
}


process BuildReference{
    label = 'BuildReference'

    publishDir "${params.outdir}/reference", mode: 'copy'
 
    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(sampleName), path(fastq)
    val fwd_primer
    val rev_primer
    val mismatches

    output:
    path("*reference.txt"), emit: ref_ch

    script:
    """
    fetch_reference.py \
        --reads ${fastq} \
        --fwd_primer ${fwd_primer} \
        --rev_primer ${rev_primer} \
        --mismatches ${mismatches} \
    """
}
