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

    publishDir "${params.outdir}/${primerName}/count/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/minimap-samtools' :
              projectDir + '/envs/minimap-samtools.yaml' }"

    input:
    val ready
    tuple val(primerName), val(sampleName), path(fastqs)
    path guides_index

    output:
    tuple val(sampleName), path("*.bam"), path(fastqs) , val(primerName) , emit: alignments
    tuple val(primerName) , path("*.txt") , emit: counts

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

        count_guides.py \${sample}.bam ${guides_index} ${sampleName} ${lenientFlag} > ${sampleName}_\${sample}_counts.txt
    done
    """
}

process CollateCounts {
    label = "CollateCounts"

    publishDir "${params.outdir}/${primerName}/count", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(primerName), path(counts)

    output:
    path "collated_counts.txt"
    path "collated_overall.txt"

    script:
    """
    collate_counts.py ${counts} > collated_counts.txt
    """
}


process CountDict{
    label = "CountDict"

    publishDir "${params.outdir}/${primerName}/count_dict/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(primerName), val(fwd_primer) ,val(rev_primer), val(sampleName) , path(fastqs)
    val mismatches

    output:
    tuple val(primerName), path("*.txt") , emit: counts_dict

    script:
    """
    for fastq in ${fastqs};
    do
       if [ -n "\$(gunzip < \${fastq} | head -c 1 | tr '\0\n' __)" ]; then
          custom_count.py \
              --reads \$fastq \
              --fwd_primer ${fwd_primer} \
              --rev_primer ${rev_primer} \
              --mismatches ${mismatches} \
              --sample ${sampleName}

       else
          touch ${sampleName}_ref-free_count.txt
       fi
    done
    """

}


process MergeDictCounts {
    label = "MergeDictCounts"

    publishDir "${params.outdir}/${primerName}/count_dict", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(primerName), path(counts)

    output:
    path "merged_filtered_counts.txt"
    path "merged_all_counts.txt"

    script:
    """
    merge_custom_count.py ${counts}
    """
}
