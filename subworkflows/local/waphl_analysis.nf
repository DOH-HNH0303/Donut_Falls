#!/usr/bin/env nextflow

include { MASH_TAXA } from '../../modules/local/mash'
include { FINAL_SUMMARY } from '../../modules/local/final_summary'

workflow WAPHL_ANALYSIS {
    take:
    ch_consensus_files    // channel: [ val(meta), file(fasta) ]
    ch_donut_summary      // channel: file(donut_falls_summary.tsv)

    main:
    ch_versions = Channel.empty()

    // Run MASH_TAXA on all consensus files
    MASH_TAXA(ch_consensus_files)
    ch_versions = ch_versions.mix(MASH_TAXA.out.versions.first())

    // Collect all mash taxa files and consensus files for final summary
    MASH_TAXA.out.taxa
        .map { meta, taxa_file -> taxa_file }
        .collect()
        .set { ch_mash_taxa_files }

    ch_consensus_files
        .map { meta, fasta_file -> fasta_file }
        .collect()
        .set { ch_consensus_fasta_files }

    // Run FINAL_SUMMARY to create enhanced summary
    FINAL_SUMMARY(
        ch_donut_summary,
        ch_mash_taxa_files,
        ch_consensus_fasta_files
    )
    ch_versions = ch_versions.mix(FINAL_SUMMARY.out.versions)

    emit:
    mash_taxa = MASH_TAXA.out.taxa
    final_summary = FINAL_SUMMARY.out.summary
    versions = ch_versions
}