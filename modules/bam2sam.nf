process PrepareForConsensus {
    label = "PrepareForConsensus"

    publishDir "${params.outdir}/count/${sampleName}"

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    tuple val(sampleName), path(bam), path(fastq)

    output:
    tuple val(sampleName), path("*.sam"), path("*.bz")

    script:
    def sample = bam.getSimpleName()
    """
    samtools view -h ${bam} > ${sample}.sam
    zcat ${fastq} | bgzip -c - > ${sample}.fastq.bz
    """
}
