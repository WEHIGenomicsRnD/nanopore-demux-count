// based on https://github.com/nf-core/modules/blob/master/modules/nf-core/medaka/main.nf

process Medaka {
    tag "${reads.getSimpleName()}"
    label 'Medaka'

    publishDir "${params.outdir}/consensus/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/medaka' :
              projectDir + '/envs/medaka.yaml' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.1--py310h87e71ce_0' :
        'biocontainers/medaka:1.11.1--py310h87e71ce_0' }"

    input:
    tuple val(sampleName), path(reads) 
    path(reference)
    val(medaka_model)

    output:
    tuple val(sampleName), path("*.fa.gz"), emit: assembly

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: reads.getSimpleName().replaceFirst("rezip_", "")
    // will create a blank fasta file if fastq file input is blank
    // as for some sample we will not be able to generate a consensus
    """
    if [ -s ${reads} ]; then
        medaka_consensus \\
            -t $task.cpus \\
            $args \\
            -i $reads \\
            -d $reference \\
            -m $medaka_model \\
            -o ./
    else
        touch consensus.fasta
    fi

    mv consensus.fasta ${prefix}.fa

    gzip -n ${prefix}.fa
    """
}


process Medaka_split {
    tag "${bam.getSimpleName()}"
    label 'Medaka'

    publishDir "${params.outdir}/consensus/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/medaka' :
              projectDir + '/envs/medaka.yaml' }"


    input:
    tuple val(sampleName), path(bam)
    path(reference)
    val(medaka_model)

    output:
    tuple val(sampleName), path("*.fa.gz"), emit: assembly

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: bam.getSimpleName()
    // will create a blank fasta file if fastq file input is blank
    // as for some sample we will not be able to generate a consensus
    """
    if [ -s ${bam} ]; then
         samtools index ${bam}
         medaka consensus $bam consensus.hdf --model ${medaka_model} --batch_size 100 --threads 2
         medaka stitch consensus.hdf ${reference} consensus.fasta --threads 2
         
    else
        touch consensus.fasta
    fi

    mv consensus.fasta ${prefix}.fa

    gzip -n ${prefix}.fa
    """
}
