include { Racon } from '../modules/racon'
include { Medaka } from '../modules/medaka'

process PrepareForConsensus {
    label = "PrepareForConsensus"
    // convert bam file to sam for racon and
    // rezip fastq file using bgzip for medaka

    publishDir "${params.outdir}/count/${sampleName}"

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    tuple val(sampleName), path(bam), path(fastq)

    output:
    tuple val(sampleName), path("*.sam"), path("rezip_*.fastq.gz"), emit: racon_input

    script:
    def sample = bam.getSimpleName()
    """
    samtools view -h ${bam} > ${sample}.sam
    zcat ${fastq} | bgzip -c - > rezip_${sample}.fastq.gz
    """
}

workflow Consensus {
    take:
        consensus_input
        reference
        medaka_model

    main:
        PrepareForConsensus(consensus_input)

        Racon(PrepareForConsensus.out.racon_input, reference)
        
        Medaka(Racon.out.racon_consensus, medaka_model)

        consensus_sequences = Medaka.out.assembly

    emit:
        consensus_sequences
}
