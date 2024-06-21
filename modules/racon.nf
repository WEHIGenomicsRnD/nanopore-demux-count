// adapted from https://github.com/nf-core/modules/blob/master/modules/nf-core/racon/main.nf

process Racon {
    tag "${sam.getSimpleName()}"
    label 'Racon'

    publishDir "${params.outdir}/consensus/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/racon-medaka' :
              projectDir + '/envs/racon-medaka.yaml' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' :
        'biocontainers/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(sampleName), path(sam), path(reads)
    path(assembly)

    output:
    tuple val(sampleName), path(reads), path('*_racon_consensus.fasta') , emit: racon_consensus
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sam.getSimpleName()
    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        "${sam}" \\
        $args \\
        "${assembly}" > \\
        ${prefix}_racon_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        racon: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
