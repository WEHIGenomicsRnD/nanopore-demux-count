process AlignPairs {
    tag "${consensusFasta.getSimpleName()}"
    label = "AlignPairs"

    publishDir "${params.outdir}/${primerName}/consensus/${sampleName}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/biopython' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    tuple val(sampleName), path(consensusFasta), val(primerName)
    path(reference)

    output:
    tuple val(sampleName), path("*.txt"), val(primerName), emit: pair_alignments

    script:
    def prefix = task.ext.prefix ?: consensusFasta.getSimpleName()
    """
    if [ -n "\$(gunzip < ${consensusFasta} | head -c 1 | tr '\0\n' __)" ]; then
        align_pairs.py ${consensusFasta} ${reference} > ${prefix}_ref_align.txt
    else
        touch ${prefix}_ref_align.txt
    fi
    """
}
