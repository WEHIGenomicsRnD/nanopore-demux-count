process IndexGuides {
    label = "IndexGuides"

    publishDir params.outdir, mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    path guides_fasta

    output:
    path "*.mmi"

    script:
    def guidesName = guides_fasta.getSimpleName()
    def outFile = "${guidesName}.mmi"
    """
    minimap2 -d ${outFile} ${guides_fasta}
    """
}

process CountGuides {
    label = "CountGuides"

    publishDir "${params.outdir}/count/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    tuple val(sampleName), path(fastqs)
    path guides_index

    output:
    path "*.bam*"
    path "*.txt"

    script:
    def outCounts = "${sampleName}_guide_counts.txt"
    """
    for fastq in ${fastqs};
    do
        sample=\${fastq%.fastq*}

        minimap2 -ax map-ont -N 1 \
            ${guides_index} \$fastq | \
            samtools view -S -b | \
            samtools sort -o \${sample}.bam

        samtools index \${sample}.bam

        samtools idxstats \${sample}.bam > \${sample}_counts.txt
    done
    """
}