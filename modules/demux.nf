// This module was written by WEHI Research Computing
// see https://github.com/WEHIGenomicsRnD/demultiplex-paired-end-library

process GenerateSelectFile {
    label 'CreateSFile'
    queue 'regular'
    cpus  2
    memory 1.GB
    time '1h'
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path primer_template_file

    output:
    path "select.txt"
    val true, emit: done

    script:
    """
    #!/usr/bin/env bash

    primer_template=\$(cat ${primer_template_file})
    fwd_primers=\$(cat ${primer_template_file} |tr -d '\r'|  awk -F ',' '{print \$1}' | grep Fwd)
    rev_primers=\$(cat ${primer_template_file} |tr -d '\r'| awk -F ',' '{print \$1}' | grep Rev)

    > select.txt
    for fwd_primer in \${fwd_primers[@]}; do
        for rev_primer in \${rev_primers[@]}; do
            echo -e "\${fwd_primer},\${rev_primer}\t\${fwd_primer}-\${rev_primer}" >> select.txt
        done
    done
    """
}

process CreateConfigFile {

        label 'CreateCFile'
        queue 'regular'
        cpus  2
        memory 1.GB
        time '1h'

        tag "${sampleId}"
        publishDir "${params.outdir}/", mode: 'copy'
        
        input:
        path(configtxt)

        output:
        path('config.txt') 
        val true, emit: done

        script:
        // add a header if not using a custom config file
        def header = params.splitcode_config_file == "" || params.splitcode_config_file == null ?
            "groups\tids\ttags\tdistances\tnext\tminFindsG\tmaxFindsG\tlocations\n" : ""
        """
        echo -e "${header}\$(cat ${configtxt})" > config.txt
        """
        
}
process SplitCode{
    label 'SplitCode'

    tag "${reads.getSimpleName()}"
    publishDir "${params.outdir}/${primerName}/split/${reads.getSimpleName()}", mode: 'copy'
    container 'oras://ghcr.io/wehi-researchcomputing/splitcode_container:latest'
   
    input:
    tuple val(primerName), val(fwd_primer) , val(rev_primer), val(sampleName), path(reads)
    val configFileReady
    val selectFileReady
    path(config)
    path(select)
    
    output:
    tuple val(primerName), val(sampleName), path("*.fastq*") , emit: fastq
    path "*.txt"

    script:
    """
    splitcode -c ${config} --keep=${select} -t ${task.cpus} --nFastqs=1 \
                --assign --summary summary.txt -o out.fastq.gz --gzip \
                --no-outb --mapping mapping.txt --seq-names \
                --mod-names --com-names --unassigned=unmapped.fastq.gz \
                ${reads}

    for file in *_0.fastq.gz; do
        newname="\${file/_0.fastq.gz/.fastq.gz}"
        mv "\$file" "\$newname"
    done
    """
}
