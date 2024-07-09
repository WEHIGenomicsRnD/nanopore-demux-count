include { Medaka } from '../modules/medaka'
include { AlignPairs } from '../modules/align_pairs.nf'

process PrepareForConsensus {
    label = "PrepareForConsensus"
    // rezip fastq file using bgzip for medaka

    publishDir "${params.outdir}/count/${sampleName}"

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    tuple val(sampleName), path(fastq)

    output:
    tuple val(sampleName), path("rezip_*.fastq.gz"), emit: rezipped

    script:
    def sample = fastq.getSimpleName()
    """
    zcat < ${fastq} | bgzip -c - > rezip_${sample}.fastq.gz
    """
}

workflow Consensus {
    take:
        consensus_input
        reference
        medaka_model

    main:
        PrepareForConsensus(consensus_input)

        Medaka(PrepareForConsensus.out.rezipped, reference, medaka_model)

        consensus_sequences = Medaka.out.assembly

        AlignPairs(consensus_sequences, reference)

        pair_alignments = AlignPairs.out.pair_alignments

    emit:
        consensus_sequences
        pair_alignments
}
