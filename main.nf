#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

println "*********************************************************"
println "*                                                       *"
println "*           Nanopore Overhang preprocess                *"
println "*      Written by Marek Cmero, WEHI Genomics R&D        *"
println "*            genomicsrnd@wehi.edu.au                    *"
println "*                                                       *"
println "*********************************************************"

include { trim_primer as trim_fwd_primer } from './modules/trim.nf'
include { trim_primer as trim_rev_primer } from './modules/trim.nf'

workflow {

    input_ch = Channel.fromPath(params.input)

    trim_fwd_primer(input_ch, params.fwd_primer, params.mismatches, 'fwd')
        .set { fwd_ch }
    
    trim_rev_primer(fwd_ch, params.rev_primer, params.mismatches, 'rev')

}