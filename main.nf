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
include { GenerateSelectFile } from './modules/demux.nf'
include { CreateConfigFile } from './modules/demux.nf'
include { SplitCode } from './modules/demux.nf'
include { IndexGuides } from './modules/count.nf'
include { CountGuides } from './modules/count.nf'
if (params.use_db) {
    include { fromQuery } from 'plugin/nf-sqldb'
}

workflow {

    Channel.fromPath("${params.input_dir}/*.{fq,fastq}{,.gz}")
             .ifEmpty {
                     error("""
                     No samples could be found! Please check whether your input directory
                     contains any fastq files (.fq or .fastq with optional .gz compression).
                     """)
         }
         .map{ it -> [it.getSimpleName(), it]}
         .groupTuple()
         .set{input_ch}

    trim_ch = TrimPrimer(input_ch,
                         params.fwd_primer,
                         params.rev_primer,
                         params.primer_mismatches,
                         params.barcode_length,
                         params.output_untrimmed)

    if (params.demultiplex) {
        if (params.splitcode_config_file != null && params.splitcode_config_file != '') {
            Channel.fromPath(params.splitcode_config_file).set{config_ch}
        } else if (params.use_db) {
            def where_ch = []
            // Construct the where clause for the query
            new File(params.index_template_file).readLines().each { line ->
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
                        def locations = direction == 'F' ? "0:0:${params.bases_num_r1}" : "0:${params.bases_num_r2}:0"

                        return "$group\t$id\t$sequence\t$distances\t$nextTag\t1\t1\t$locations"
                }
                .collectFile(name: 'config.txt', newLine: true).set{config_ch}
        } else {
            // build the config file from the index template
            def indexes = []
            new File(params.index_template_file).readLines().each { line ->
                if (!line.startsWith('index_name')) {
                    def index = line.trim().split(',').each { it.trim() }
                    def id = index[0]
                    def direction = index[1]
                    def group = direction == "F" ? "Fwd" : "Rev"
                    def sequence = index[2]
                    def distances = direction == 'F' ? "${params.idx_5p_mismatch}" : "${params.idx_3p_mismatch}"
                    def nextTag = direction == 'F' ? '{{Rev}}' : '-'
                    def locations = direction == 'F' ? "0:0:${params.bases_num_r1}" : "0:${params.bases_num_r2}:0"

                    indexes << "$group\t$id\t$sequence\t$distances\t$nextTag\t1\t1\t$locations"
                }
            }
            Channel.from( indexes )
                .collectFile(name: 'config.txt', newLine: true).set{config_ch}
        }
        CreateConfigFile(config_ch).set{configFile}
        GenerateSelectFile(file(params.index_template_file)).set{selectTxt}
        SplitCode(trim_ch.trimmed_ch,
                  configFile.done.first(),
                  selectTxt.done,
                  file("${params.outdir}/config.txt"),
                  file("${params.outdir}/select.txt")).fastq.set{demux_ch}        
    } else {
        trim_ch.trimmed_ch.set{demux_ch}
    }

    if (params.guides_fasta != null && params.guides_fasta != '') {
        def guidesIndex = file(params.guides_fasta).getSimpleName() + ".mmi"
        IndexGuides(params.guides_fasta).set{index_ch}
        CountGuides(index_ch.done, demux_ch, file("${params.outdir}/${guidesIndex}"))
    }
}
