#!/usr/bin/env nextflow

include { MASH_TAXA } from '../../modules/local/mash'
include { MASH_HUMAN_CONTAMINATION } from '../../modules/local/mash'
include { FINAL_SUMMARY } from '../../modules/local/final_summary'
include { COVERAGE_ANALYSIS } from '../../modules/local/coverage_analysis'

workflow WAPHL_ANALYSIS {
    take:
    ch_consensus_meta     // channel: [ val(meta), file(fasta) ] - for MASH_TAXA
    ch_consensus_files    // channel: file(fasta) - from copy process, includes sub_fasta files
    ch_donut_summary      // channel: file(donut_falls_summary.tsv)
    ch_nanopore_input     // channel: [ val(meta), file(nanopore_fastq) ] - for COVERAGE_ANALYSIS
    ch_illumina_input     // channel: [ val(meta), [file(R1), file(R2)] ] - for COVERAGE_ANALYSIS

    main:
    ch_versions = Channel.empty()

    // Filter consensus files to only process the final/best version for each sample
    // Priority: pypolca > polypolish > clair3 > reoriented > unicycler
    ch_consensus_meta
        .groupTuple(by: 0)
        .map { meta, fastas ->
            // Find the best fasta file for this sample
            def best_fasta = null
            def priority_order = ['pypolca', 'polypolish', 'clair3', 'reoriented', 'unicycler']
            
            for (priority in priority_order) {
                def matching_fasta = fastas.find { it.name.contains("_${priority}.fasta") }
                if (matching_fasta) {
                    best_fasta = matching_fasta
                    break
                }
            }
            
            // If no priority match found, take the first one
            if (!best_fasta) {
                best_fasta = fastas[0]
            }
            
            tuple(meta, best_fasta)
        }
        .set { ch_consensus_final }

    // Run MASH_TAXA on final consensus files only
    MASH_TAXA(ch_consensus_final)
    ch_versions = ch_versions.mix(MASH_TAXA.out.versions.first())

    // Run MASH_HUMAN_CONTAMINATION on final consensus files only
    MASH_HUMAN_CONTAMINATION(ch_consensus_final)
    ch_versions = ch_versions.mix(MASH_HUMAN_CONTAMINATION.out.versions.first())

    // Run COVERAGE_ANALYSIS on consensus files with available reads (ONT or Illumina)
    // Combine ONT and Illumina channels, prioritizing ONT if both are available
    ch_nanopore_input
        .map { meta, reads -> tuple(meta, reads, 'ont') }
        .mix(
            ch_illumina_input
                .map { meta, reads -> tuple(meta, reads, 'illumina') }
        )
        .groupTuple(by: 0)
        .map { meta, reads_list, type_list ->
            // If both ONT and Illumina are available, prioritize ONT
            def ont_idx = type_list.findIndexOf { it == 'ont' }
            def selected_reads = ont_idx >= 0 ? reads_list[ont_idx] : reads_list[0]
            tuple(meta, selected_reads)
        }
        .join(ch_consensus_final, by: 0, remainder: false)
        .map { meta, reads, fasta -> tuple(meta, fasta, reads) }
        .set { ch_coverage_input }

    COVERAGE_ANALYSIS(ch_coverage_input)
    ch_versions = ch_versions.mix(COVERAGE_ANALYSIS.out.versions)




    // Collect all mash taxa files for final summary
    MASH_TAXA.out.taxa
        .map { meta, taxa_file -> taxa_file }
        .collect()
        .ifEmpty([])
        .set { ch_mash_taxa_files }

    // Collect all coverage analysis files for final summary
    COVERAGE_ANALYSIS.out.summary
        .map { meta, coverage_file -> coverage_file }
        .collect()
        .ifEmpty([])
        .set { ch_coverage_files }

    // Collect all human contamination files for final summary
    MASH_HUMAN_CONTAMINATION.out.human_summary
        .map { meta, human_file -> human_file }
        .collect()
        .ifEmpty([])
        .set { ch_human_contamination_files }

    // Collect all consensus files (including sub_fasta files) for final summary
    ch_consensus_files
        .collect()
        .set { ch_all_consensus_files }


    // Run FINAL_SUMMARY to create enhanced summary
    FINAL_SUMMARY(
        ch_donut_summary,
        ch_mash_taxa_files,
        ch_all_consensus_files,
        ch_coverage_files,
        ch_human_contamination_files
    )
    ch_versions = ch_versions.mix(FINAL_SUMMARY.out.versions)

    emit:
    mash_taxa = MASH_TAXA.out.taxa
    coverage_analysis = COVERAGE_ANALYSIS.out.summary
    human_contamination = MASH_HUMAN_CONTAMINATION.out.human_contamination
    human_summary = MASH_HUMAN_CONTAMINATION.out.human_summary
    final_summary = FINAL_SUMMARY.out.summary
    versions = ch_versions
}