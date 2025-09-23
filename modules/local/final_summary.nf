process FINAL_SUMMARY {
  tag           "Creating final summary with WAPHL analysis"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.30'
  time          '30m'
  
  input:
  file(donut_falls_summary)
  file(mash_taxa_files)
  file(consensus_files)
  file(coverage_files)

  output:
  path "waphl_final_summary.tsv", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  """
  #!/usr/bin/env python3
  import csv
  import glob
  import os
  import sys
  from pathlib import Path

  def read_tsv_file(filename):
      \"\"\"Read TSV file and return list of dictionaries\"\"\"
      data = []
      with open(filename, 'r', newline='') as f:
          reader = csv.DictReader(f, delimiter='\\t')
          for row in reader:
              data.append(row)
      return data

  def write_tsv_file(filename, data, fieldnames):
      \"\"\"Write data to TSV file\"\"\"
      with open(filename, 'w', newline='') as f:
          writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\\t')
          writer.writeheader()
          writer.writerows(data)

  def determine_assembly_type(sample_row):
      # Determine if assembly is illumina/ont/hybrid based on available data
      
      # Check for illumina data
      has_illumina = False
      for col, value in sample_row.items():
          if 'illumina' in col.lower() and value and str(value).strip() != '':
              has_illumina = True
              break
      
      # Check for nanopore data
      has_nanopore = False
      for col, value in sample_row.items():
          if 'nanopore' in col.lower() and value and str(value).strip() != '':
              has_nanopore = True
              break
      
      # Check for unicycler (hybrid assembly)
      has_unicycler = False
      for col, value in sample_row.items():
          if 'unicycler' in col.lower() and value and str(value).strip() != '':
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

  def get_mash_taxa_and_distance(sample, mash_files):
      # Extract the top mash taxa hit and distance for a sample
      for mash_file in mash_files:
          if sample in mash_file:
              try:
                  with open(mash_file, 'r') as f:
                      lines = f.readlines()
                      if len(lines) > 1:  # Skip header
                          # Get the first data line (best match)
                          data_line = lines[1].strip().split('\\t')
                          if len(data_line) >= 6:
                              genus_species = data_line[5]  # genus_species column
                              distance = data_line[2]       # distance column
                              return genus_species, distance
              except Exception as e:
                  print("Error reading {}: {}".format(mash_file, e))
                  continue
      return 'unknown', 'unknown'


  def get_coverage_data(sample, coverage_files):
      # Extract coverage data for a sample (ONT or Illumina)
      coverage_data = {}
      
      for coverage_file in coverage_files:
          if sample in coverage_file and '_coverage_summary.tsv' in coverage_file:
              try:
                  with open(coverage_file, 'r') as f:
                      lines = f.readlines()
                      read_platform = None
                      
                      # Parse the main coverage statistics
                      if len(lines) > 1:  # Skip header
                          # Get the first data line
                          data_line = lines[1].strip().split('\\t')
                          if len(data_line) >= 7:
                              read_platform = data_line[1]
                              coverage_data.update({
                                  f'{read_platform.lower()}_total_bases': data_line[2],
                                  f'{read_platform.lower()}_covered_bases': data_line[3], 
                                  f'{read_platform.lower()}_coverage_breadth_percent': data_line[4],
                                  f'{read_platform.lower()}_mean_depth': data_line[5],
                                  f'{read_platform.lower()}_max_depth': data_line[6],
                                  'read_platform': read_platform
                              })
                      
                      # Parse theoretical coverage data
                      for line in lines:
                          if line.startswith('genome_size'):
                              continue
                          elif '\\t' in line and not line.startswith('#') and 'genome_size' in line:
                              # This is the theoretical coverage line
                              parts = line.strip().split('\\t')
                              if len(parts) >= 4:
                                  platform = parts[3].strip()
                                  coverage_data.update({
                                      f'{platform.lower()}_genome_size': parts[0],
                                      f'{platform.lower()}_total_read_bases': parts[1],
                                      f'{platform.lower()}_theoretical_coverage': parts[2]
                                  })
                                  break
                      
                      # If we found data, return it
                      if coverage_data:
                          return coverage_data
                          
              except Exception as e:
                  print("Error reading {}: {}".format(coverage_file, e))
                  continue
      return coverage_data


  def get_consensus_filepath(sample, consensus_files):
      # Get the consensus genome filepath for a sample
      # Prioritize the most processed version available
      
      # First priority: pypolca (final polishing step)
      for consensus_file in consensus_files:
          if sample in consensus_file and '_pypolca.fasta' in consensus_file:
              return consensus_file
      
      # Second priority: polypolish
      for consensus_file in consensus_files:
          if sample in consensus_file and '_polypolish.fasta' in consensus_file:
              return consensus_file
      
      # Third priority: clair3
      for consensus_file in consensus_files:
          if sample in consensus_file and '_clair3.fasta' in consensus_file:
              return consensus_file
      
      # Fourth priority: reoriented
      for consensus_file in consensus_files:
          if sample in consensus_file and '_reoriented.fasta' in consensus_file:
              return consensus_file
      
      # Fifth priority: unicycler (hybrid assembly)
      for consensus_file in consensus_files:
          if sample in consensus_file and '_unicycler.fasta' in consensus_file:
              return consensus_file
      
      # Fallback: any fasta file for the sample
      for consensus_file in consensus_files:
          if sample in consensus_file and consensus_file.endswith('.fasta'):
              return consensus_file
      
      return ''

  def get_sub_fasta_filepath(sample, consensus_files):
      # Get the sub_fasta filepath if it exists
      # First priority: Look for files with 'sub_' prefix (circular assemblies with special headers)
      for consensus_file in consensus_files:
          if sample in consensus_file and 'sub_' in consensus_file and consensus_file.endswith('.fasta'):
              return consensus_file
      
      # Second priority: Look for reoriented files as they are submission-ready
      for consensus_file in consensus_files:
          if sample in consensus_file and '_reoriented' in consensus_file and consensus_file.endswith('.fasta'):
              return consensus_file
      
      # If neither sub_ files nor reoriented files are found, return empty string
      return ''

  def main():
      try:
          # Read the original donut falls summary
          summary_data = read_tsv_file('${donut_falls_summary}')
          
          # Get all mash taxa files
          mash_files = glob.glob('*_taxa.txt')
          
          # Get all consensus files
          consensus_files = glob.glob('*.fasta')

          # Get all coverage analysis files
          coverage_files = glob.glob('*_coverage_summary.tsv')
          
          print("Found {} samples in summary".format(len(summary_data)))
          print("Found {} mash taxa files: {}".format(len(mash_files), mash_files))
          print("Found {} consensus files: {}".format(len(consensus_files), consensus_files))
          print("Found {} coverage analysis files: {}".format(len(coverage_files), coverage_files))
          
          # Process each sample
          for row in summary_data:
              sample = row['sample']
              print("Processing sample: {}".format(sample))
              
              # Get mash taxa and distance
              mash_taxa, mash_distance = get_mash_taxa_and_distance(sample, mash_files)
              row['mash_taxa'] = mash_taxa
              row['mash_distance'] = mash_distance

              # Get coverage data (ONT or Illumina)
              coverage_data = get_coverage_data(sample, coverage_files)
              for key, value in coverage_data.items():
                  row[key] = value
              
              # Determine assembly type
              row['assembly_type'] = determine_assembly_type(row)
              
              # Get consensus filepath
              row['consensus_filepath'] = get_consensus_filepath(sample, consensus_files)
              
              # Get sub_fasta filepath
              row['reoriented_assembly_filepath'] = get_sub_fasta_filepath(sample, consensus_files)
          
          # Get all fieldnames (original + new columns)
          if summary_data:
              fieldnames = list(summary_data[0].keys())
          else:
              fieldnames = ['sample']
          
          # Write the enhanced summary
          write_tsv_file('waphl_final_summary.tsv', summary_data, fieldnames)
          
          print("Enhanced summary created with {} samples".format(len(summary_data)))
          print("Found {} mash taxa files".format(len(mash_files)))
          print("Found {} consensus files".format(len(consensus_files)))
          print("Found {} coverage analysis files".format(len(coverage_files)))
          
          # Write versions file
          with open('versions.yml', 'w') as f:
              f.write('"${task.process}":\\n')
              f.write('  python: {}\\n'.format(sys.version.split()[0]))
          
      except Exception as e:
          print("Error in main: {}".format(e))
          import traceback
          traceback.print_exc()
          sys.exit(1)

  if __name__ == '__main__':
      main()
  """
}