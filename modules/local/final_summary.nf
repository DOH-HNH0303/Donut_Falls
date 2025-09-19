process FINAL_SUMMARY {
  tag           "Creating final summary with WAPHL analysis"
  label         "process_low"
  publishDir    "${params.outdir}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.30'
  time          '30m'
  
  input:
  file(donut_falls_summary)
  file(mash_taxa_files)
  file(consensus_files)

  output:
  path "waphl_final_summary.tsv", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import glob
  import os
  from pathlib import Path

  def determine_assembly_type(sample, summary_df):
      '''Determine if assembly is illumina/ont/hybrid based on available data'''
      if sample not in summary_df['sample'].values:
          return 'unknown'
      
      sample_row = summary_df[summary_df['sample'] == sample].iloc[0]
      
      # Check for illumina data
      has_illumina = False
      for col in sample_row.index:
          if 'illumina' in col.lower() and pd.notna(sample_row[col]) and str(sample_row[col]).strip() != '':
              has_illumina = True
              break
      
      # Check for nanopore data
      has_nanopore = False
      for col in sample_row.index:
          if 'nanopore' in col.lower() and pd.notna(sample_row[col]) and str(sample_row[col]).strip() != '':
              has_nanopore = True
              break
      
      # Check for unicycler (hybrid assembly)
      has_unicycler = False
      for col in sample_row.index:
          if 'unicycler' in col.lower() and pd.notna(sample_row[col]) and str(sample_row[col]).strip() != '':
              has_unicycler = True
              break
      
      if has_unicycler:
          return 'hybrid'
      elif has_illumina and has_nanopore:
          return 'hybrid'
      elif has_nanopore:
          return 'ont'
      elif has_illumina:
          return 'illumina'
      else:
          return 'unknown'

  def get_mash_taxa(sample, mash_files):
      '''Extract the top mash taxa hit for a sample'''
      for mash_file in mash_files:
          if sample in mash_file:
              try:
                  with open(mash_file, 'r') as f:
                      lines = f.readlines()
                      if len(lines) > 1:  # Skip header
                          # Get the first data line (best match)
                          data_line = lines[1].strip().split('\\t')
                          if len(data_line) >= 6:
                              return data_line[5]  # genus_species column
              except Exception as e:
                  print("Error reading {}: {}".format(mash_file, e))
                  continue
      return 'unknown'

  def get_consensus_filepath(sample, consensus_files):
      '''Get the consensus genome filepath for a sample'''
      for consensus_file in consensus_files:
          if sample in consensus_file and consensus_file.endswith('.fasta'):
              return consensus_file
      return ''

  def get_sub_fasta_filepath(sample, consensus_files):
      '''Get the sub_fasta filepath if it exists'''
      for consensus_file in consensus_files:
          if sample in consensus_file and 'sub_' in consensus_file and consensus_file.endswith('.fasta'):
              return consensus_file
      return ''

  def main():
      # Read the original donut falls summary
      summary_df = pd.read_csv('${donut_falls_summary}', sep='\\t')
      
      # Get all mash taxa files
      mash_files = glob.glob('*_taxa.txt')
      
      # Get all consensus files
      consensus_files = glob.glob('*.fasta')
      
      # Create new columns
      summary_df['mash_taxa'] = ''
      summary_df['assembly_type'] = ''
      summary_df['consensus_filepath'] = ''
      summary_df['sub_fasta_filepath'] = ''
      
      # Fill in the new columns for each sample
      for idx, row in summary_df.iterrows():
          sample = row['sample']
          
          # Get mash taxa
          summary_df.at[idx, 'mash_taxa'] = get_mash_taxa(sample, mash_files)
          
          # Determine assembly type
          summary_df.at[idx, 'assembly_type'] = determine_assembly_type(sample, summary_df)
          
          # Get consensus filepath
          summary_df.at[idx, 'consensus_filepath'] = get_consensus_filepath(sample, consensus_files)
          
          # Get sub_fasta filepath
          summary_df.at[idx, 'sub_fasta_filepath'] = get_sub_fasta_filepath(sample, consensus_files)
      
      # Save the enhanced summary
      summary_df.to_csv('waphl_final_summary.tsv', sep='\\t', index=False)
      
      print("Enhanced summary created with {} samples".format(len(summary_df)))
      print("Found {} mash taxa files".format(len(mash_files)))
      print("Found {} consensus files".format(len(consensus_files)))

  if __name__ == '__main__':
      main()

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(python --version | awk '{print \$2}')
  END_VERSIONS
  """
}