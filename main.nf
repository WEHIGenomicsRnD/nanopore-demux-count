#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

println "*********************************************************"
println "*                                                       *"
println "*           Nanopore Overhang process                   *"
println "*      Written by Marek Cmero, WEHI Genomics R&D        *"
println "*            genomicsrnd@wehi.edu.au                    *"
println "*                                                       *"
println "*********************************************************"

include { TrimPrimer } from './modules/trim.nf'
include { BuildReference } from './modules/trim.nf'
include { MergeFastq } from './modules/trim.nf'
include { GenerateSelectFile } from './modules/demux.nf'
include { CreateConfigFile } from './modules/demux.nf'
include { SplitCode } from './modules/demux.nf'
include { IndexGuides } from './modules/count.nf'
include { CountGuides } from './modules/count.nf'
include { CollateCounts } from './modules/count.nf'
include { MergeDictCounts } from './modules/count.nf'
include { CountDict } from './modules/count.nf'
include { Consensus } from './subworkflows/consensus'
include { ProcessExcel } from './modules/preprocess.nf'

if (params.use_db) {
    include { fromQuery } from 'plugin/nf-sqldb'
}

workflow {

    index_ch = Channel.empty()
    Channel
        .fromPath( "${params.input_dir}/**fastq.gz" )
           .ifEmpty {
                     error("""
                     No sample folder could be found! Please check whether your input directory
                     contains barcode* folder for concatenating the files for each sample.
                     """)
           }
           .flatten()
           .map { file-> [file.parent.name, file] }
           .groupTuple()
           .filter{ Sid, file ->
              !Sid.startsWith("unclassified")
           }.set{rawfastq_ch}


    // Merge fastq files
    merge_ch = MergeFastq(rawfastq_ch)

    merge_ch.map{ it -> [it.getSimpleName(), it]}
      .groupTuple()
      .set{base_input_ch}

    // Checking if excel file is present
    if (params.excel_file != ""){
       ProcessExcel(params.excel_file).set{excel_ch}

       excel_ch.primer_ch.splitCsv(header:true)
          .map { row -> tuple(row.name,row.fwd_primer,row.rev_primer)}
          .set{primer_list_ch}
  
       excel_ch.index_ch.set{index_file_ch}


    }else if( params.excel_file == ""){
        if (params.index_template_file != ''){
              index_file_ch = Channel.fromPath(params.index_template_file)
              primer_list_ch = Channel.of(["primer1",params.fwd_primer,params.rev_primer ])
        }else{
           error "Index Template file is required {params.index_template_file)"
        }
    }
   

    // Combining the input fastq and primer channel    
    primer_list_ch.combine(base_input_ch).set{base_ch}

    base_ch
         .map { primerName, fwd_primer, rev_primer, sampleName, fastq_files ->
             tuple( primerName, sampleName, fastq_files )
         }.set{input_ch}

    if (params.count_only && !params.demultiplex) {
        input_ch.set{ demux_ch }
    } else {
        trim_ch = TrimPrimer(base_ch,
                             params.primer_mismatches,
                             params.barcode_length,
                             params.output_untrimmed)
    }
 
    if (params.demultiplex) {
        if (params.splitcode_config_file != null && params.splitcode_config_file != '') {
            Channel.fromPath(params.splitcode_config_file).set{config_ch}
        } else if (params.use_db) {
            def where_ch = []
            // Construct the where clause for the query
            new File(excel_ch.index_ch).readLines().each { line ->
                if (line.trim() != 'index_name') {
                    where_ch << "'${line.trim()}'"
                }
            }
            def where_clause = where_ch.join(",")
            def query = """SELECT index_name,
                                  index_sequence,
                                  index_sequence_rc,
                                  index_direction
                            FROM amplicon_index
                            WHERE index_name IN (${where_clause});"""

            Channel.fromQuery(query, db: 'my-db', batchSize:100)
                .map { index ->
                        def id = index[0]
                        def direction = index[3]
                        def group = direction == "F" ? "Fwd" : "Rev"
                        def sequence = direction == "F" ? index[1] : index[2] // index_sequence or index_sequence_rc
                        def distances = direction == 'F' ? "${params.idx_5p_mismatch}" : "${params.idx_3p_mismatch}"
                        def nextTag = direction == 'F' ? '{{Rev}}' : '-'
                        def locations = direction == 'F' ? "0:0:${params.barcode_length}" : "0:-${params.barcode_length}:0"

                        return "$group\t$id\t$sequence\t$distances\t$nextTag\t1\t1\t$locations"
                }
                .collectFile(name: 'config_tmp.txt', newLine: true).set{config_ch}
        } else {
            // build the config file from the index template
            def indexes = []
            index_file_ch.view() 
            index_file_ch.splitCsv(header:true)
                   .map { row ->
                    def id = row.index_name
                    def direction = row.direction
                    def group = direction == "F" ? "Fwd" : "Rev"
                    def sequence = row.sequence
                    def distances = direction == 'F' ? "${params.idx_5p_mismatch}" : "${params.idx_3p_mismatch}"
                    def nextTag = direction == 'F' ? '{{Rev}}' : '-'
                    def locations = direction == 'F' ? "0:0:${params.barcode_length}" : "0:-${params.barcode_length}:0"

                    return "$group\t$id\t$sequence\t$distances\t$nextTag\t1\t1\t$locations"
                
            }
            .collectFile(name: 'config_tmp.txt', newLine: true).set{config_ch}
        }

        CreateConfigFile(config_ch).set{configFile}
        GenerateSelectFile(index_file_ch).set{selectTxt}
        SplitCode(trim_ch.trimmed_ch,
                  configFile.done.first(),
                  selectTxt.done,
                  file("${params.outdir}/config.txt"),
                  file("${params.outdir}/select.txt")).fastq.set{demux_ch}
    } else if (! params.count_only) {
        trim_ch.trimmed_ch.set{demux_ch}
    }

    if (params.count_dict) {
       primer_list_ch.combine(demux_ch, by: 0).set{count_dict_ch}
       count_dict_ch.view()
       CountDict(count_dict_ch,params.primer_mismatches).set{dict_ch}
       MergeDictCounts(dict_ch.counts_dict)
       
    }


    if (params.count_only | params.consensus) {
        if (params.guides_fasta != ''){
              guides_ch = Channel.fromPath(params.guides_fasta)
        }else{
              excel_ch.guides.set{guides_ch}
        }

        IndexGuides(guides_ch).set{index_ch}
        CountGuides(index_ch.done, demux_ch, guides_ch).set{count_ch}
        CollateCounts(count_ch.counts)

        if (params.consensus) {
            if (params.count_only) {
                demux_ch.set{ fastq_ch }
            } else {
                // in this case fastqFiles are arrays
                // not singular elements per sample
                // reformat channel for consensus input
                // to tuple (sampleName, fastqFile)
                // filter out unmapped and out files from splitcode
                demux_ch.flatMap { sample ->
                    def (primerName, sampleName, fastqFiles) = sample
                    return fastqFiles.indices.collect { index ->
                        [primerName, sampleName, fastqFiles[index]]
                    }
                }.filter{ primerName, sampleName, fastqFile ->
                    !fastqFile.getName().startsWith("unmapped.fastq") &&
                    !fastqFile.getName().startsWith("out.fastq")
                }.set{ fastq_ch }
            }
            Consensus(fastq_ch, guides_ch, params.medaka_model)
        }
    }
}
