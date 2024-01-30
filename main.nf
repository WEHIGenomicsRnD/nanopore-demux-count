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

    input_ch = Channel.fromPath(params.input)

    trim_primer(input_ch, params.fwd_primer, params.rev_primer, params.mismatches, params.barcode_length, params.trim_log)

}