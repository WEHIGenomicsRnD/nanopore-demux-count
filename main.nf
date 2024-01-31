#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

println "*********************************************************"
println "*                                                       *"
println "*           Nanopore Overhang preprocess                *"
println "*      Written by Marek Cmero, WEHI Genomics R&D        *"
println "*            genomicsrnd@wehi.edu.au                    *"
println "*                                                       *"
println "*********************************************************"

include { trim_primer } from './modules/trim.nf'

workflow {

    Channel.fromPath("${params.input_dir}/*.{fq,fastq}{,.gz}")
             .ifEmpty {
                     error("""
                     No samples could be found! Please check whether your input directory
                     contains any fastq files (.fq or .fastq with optional .gz compression).
                     """)
         }.set{input_ch}

    trim_primer(input_ch, params.fwd_primer, params.rev_primer, params.mismatches, params.barcode_length, params.log)

}