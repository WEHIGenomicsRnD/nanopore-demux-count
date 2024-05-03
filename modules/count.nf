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
    path "*.txt"

    script:
    def outCounts = "${sampleName}_guide_counts.txt"
    /*
    count collation is a hacky bash script to get collated output,
    the script pastes the count files together, and then cuts out
    the relevant count fields (keep the first as a label column)
    we will want to make this into a nice python script later
    */
    """
    for fastq in ${fastqs};
    do
        sample=\${fastq%.fastq*}

        minimap2 -ax map-ont -N 1 \
            ${guides_index} \$fastq | \
            samtools view -S -b | \
            samtools sort -o \${sample}.bam

        count_guides.py \${sample}.bam ${params.guides_fasta} > \${sample}_counts.txt
    done

    paste *_counts.txt > tmpfile
    ncounts=\$(expr \$(ls *_counts.txt | wc -l | xargs) \\* 2)
    fields_to_cut=\$(echo 1 \$(seq 2 2 \$ncounts) | sed 's/ /,/g')
    cut -f \$fields_to_cut tmpfile > ${outCounts}
    """
}
