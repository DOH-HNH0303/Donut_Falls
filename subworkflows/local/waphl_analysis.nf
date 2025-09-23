#!/usr/bin/env nextflow

include { MASH_TAXA } from '../../modules/local/mash'
include { FINAL_SUMMARY } from '../../modules/local/final_summary'
include { ONT_COVERAGE } from '../../modules/local/ont_coverage'

workflow WAPHL_ANALYSIS {
    take:
    ch_consensus_meta     // channel: [ val(meta), file(fasta) ] - for MASH_TAXA
    ch_consensus_files    // channel: file(fasta) - from copy process, includes sub_fasta files
    ch_donut_summary      // channel: file(donut_falls_summary.tsv)
    ch_nanopore_input     // channel: [ val(meta), file(nanopore_fastq) ] - for ONT_COVERAGE

    main:
    ch_versions = Channel.empty()

    // Run MASH_TAXA on consensus files with meta information
    MASH_TAXA(ch_consensus_meta)
    ch_versions = ch_versions.mix(MASH_TAXA.out.versions.first())

    // Run ONT_COVERAGE on consensus files with nanopore reads
    ch_consensus_meta
        .join(ch_nanopore_input, by: 0, remainder: false)
        .set { ch_ont_coverage_input }

    ONT_COVERAGE(ch_ont_coverage_input)
    ch_versions = ch_versions.mix(ONT_COVERAGE.out.versions)




    // Collect all mash taxa files for final summary
    MASH_TAXA.out.taxa
        .map { meta, taxa_file -> taxa_file }
        .collect()
        .ifEmpty([])
        .set { ch_mash_taxa_files }

    // Collect all ONT coverage files for final summary
    ONT_COVERAGE.out.summary
        .map { meta, coverage_file -> coverage_file }
        .collect()
        .ifEmpty([])
        .set { ch_ont_coverage_files }

    // Collect all consensus files (including sub_fasta files) for final summary
    ch_consensus_files
        .collect()
        .set { ch_all_consensus_files }


    // Run FINAL_SUMMARY to create enhanced summary
    FINAL_SUMMARY(
        ch_donut_summary,
        ch_mash_taxa_files,
        ch_all_consensus_files,
        ch_ont_coverage_files
    )
    ch_versions = ch_versions.mix(FINAL_SUMMARY.out.versions)

    emit:
    mash_taxa = MASH_TAXA.out.taxa
    ont_coverage = ONT_COVERAGE.out.summary
    final_summary = FINAL_SUMMARY.out.summary
    versions = ch_versions
}