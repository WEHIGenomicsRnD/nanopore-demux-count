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
    val true, emit: done

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
    val ready
    tuple val(sampleName), path(fastqs)
    path guides_index

    output:
    path "*.bam*"
    path "*.txt", emit: counts

    script:
    def lenientFlag = params.lenient_counts ? "--lenient" : ""
    def extraThreads = task.cpus - 1
    """
    for fastq in ${fastqs};
    do
        sample=\${fastq%.fastq*}

        minimap2 -ax map-ont -N 1 -t ${task.cpus} -f ${params.minimap_f} \
            ${guides_index} \$fastq | \
            samtools view -S -b -@ ${extraThreads} | \
            samtools sort -@ ${extraThreads} -o \${sample}.bam

        count_guides.py \${sample}.bam ${params.guides_fasta} ${sampleName} ${lenientFlag} > ${sampleName}_\${sample}_counts.txt
    done
    """
}

process CollateCounts {
    label = "CollateCounts"

    publishDir "${params.outdir}/count", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    path counts

    output:
    path "collated_counts.txt"

    script:
    """
    collate_counts.py ${counts} > collated_counts.txt
    """
}
