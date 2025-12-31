#!/usr/bin/env nextflow

include { SOURMASH_TAXA } from '../../modules/local/sourmash'
include { MASH_HUMAN_CONTAMINATION } from '../../modules/local/mash'
include { FINAL_SUMMARY } from '../../modules/local/final_summary'
include { COVERAGE_ANALYSIS } from '../../modules/local/coverage_analysis'

workflow WAPHL_ANALYSIS {
    take:
    ch_consensus_meta     // channel: [ val(meta), file(fasta) ] - for SOURMASH_TAXA
    ch_consensus_files    // channel: file(fasta) - from copy process, includes sub_fasta files
    ch_donut_summary      // channel: file(donut_falls_summary.tsv)
    ch_nanopore_input     // channel: [ val(meta), file(nanopore_fastq) ] - for COVERAGE_ANALYSIS
    ch_illumina_input     // channel: [ val(meta), [file(R1), file(R2)] ] - for COVERAGE_ANALYSIS

    main:
    ch_versions = Channel.empty()

    // Select only pypolca assemblies (the best/final polished version)
    // Input comes as grouped tuples: [meta, [fasta1, fasta2, ...]]
    // We need to flatten and filter for pypolca files only
    ch_consensus_final = ch_consensus_meta
        .flatMap { meta, fastas ->
            // Ensure fastas is a list
            def fasta_list = fastas instanceof List ? fastas : [fastas]
            // Find pypolca files and return as individual tuples
            fasta_list.findAll { fasta -> fasta.name.contains('_pypolca.fasta') }
                      .collect { fasta -> tuple(meta, fasta) }
        }

    // Prepare database for SOURMASH_TAXA - handle remote URLs (S3, HTTP, etc.) and local files
    // Use collect() to create a value channel that can be consumed multiple times
    def database_path = params.sourmash_db_taxa ?: params.sourmash_db ?: null
    ch_database = database_path ? Channel.fromPath(database_path, checkIfExists: false).collect() : Channel.value([])
    
    // Run SOURMASH_TAXA on final consensus files only
    SOURMASH_TAXA(ch_consensus_final, ch_database)
    ch_versions = ch_versions.mix(SOURMASH_TAXA.out.versions.first())

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




    // Collect all sourmash taxa files for final summary
    SOURMASH_TAXA.out.taxa
        .map { _meta, taxa_file -> taxa_file }
        .collect()
        .ifEmpty([])
        .set { ch_mash_taxa_files }

    // Collect all coverage analysis files for final summary
    COVERAGE_ANALYSIS.out.summary
        .map { _meta, coverage_file -> coverage_file }
        .collect()
        .ifEmpty([])
        .set { ch_coverage_files }

    // Collect all human contamination files for final summary
    MASH_HUMAN_CONTAMINATION.out.human_summary
        .map { _meta, human_file -> human_file }
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
    mash_taxa = SOURMASH_TAXA.out.taxa
    coverage_analysis = COVERAGE_ANALYSIS.out.summary
    human_contamination = MASH_HUMAN_CONTAMINATION.out.human_contamination
    human_summary = MASH_HUMAN_CONTAMINATION.out.human_summary
    final_summary = FINAL_SUMMARY.out.summary
    versions = ch_versions
}