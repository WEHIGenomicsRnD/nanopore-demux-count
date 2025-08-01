process ProcessExcel{
    label = 'ProcessExcel'

    publishDir "${params.outdir}", mode: 'copy'

    conda "${ params.conda_env_location != null && params.conda_env_location != '' ?
              params.conda_env_location + '/plasmid-env' :
              projectDir + '/envs/biopython.yaml' }"

    input:
    path(excel)

    output:
    path("index.csv"), emit: index_ch
    path("primer.csv"), emit: primer_ch
    path("guides.fa"), optional:true, emit: guides

    script:
    """
       convert_xlsx.py -i ${excel}
    """
}
