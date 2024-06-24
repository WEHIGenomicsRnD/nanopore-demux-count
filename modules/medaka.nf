// based on https://github.com/nf-core/modules/blob/master/modules/nf-core/medaka/main.nf

process Medaka {
    tag "${reads.getSimpleName()}"
    label 'Medaka'

    publishDir "${params.outdir}/consensus/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/racon-medaka' :
              projectDir + '/envs/racon-medaka.yaml' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.1--py310h87e71ce_0' :
        'biocontainers/medaka:1.11.1--py310h87e71ce_0' }"

    input:
    tuple val(sampleName), path(reads), path(assembly)
    val(medaka_model)

    output:
    tuple val(sampleName), path("*.fa.gz"), emit: assembly

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: reads.getSimpleName().replaceFirst("rezip_", "")
    // will create a blank fasta file if racon assembly is blank
    // as for some sample we will not be able to generate a consensus
    """
    if [ -s ${assembly} ]; then
        medaka_consensus \\
            -t $task.cpus \\
            $args \\
            -i $reads \\
            -d $assembly \\
            -m $medaka_model \\
            -o ./
    else
        touch consensus.fasta
    fi

    mv consensus.fasta ${prefix}.fa

    gzip -n ${prefix}.fa
    """
}
