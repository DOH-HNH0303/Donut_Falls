process FINAL_SUMMARY {
  tag           "Creating final summary with WAPHL analysis"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.30'
  time          '30m'
  
  input:
  path(donut_falls_summary)
  path('sourmash_taxa/*', stageAs: 'sourmash_taxa/*')
  path('consensus/*', stageAs: 'consensus/*')
  path('coverage/*', stageAs: 'coverage/*')
  path('human_contamination/*', stageAs: 'human_contamination/*')

  output:
  path "waphl_final_summary.tsv", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def ont_cov_min = params.ont_coverage_min ?: 40
  def ill_cov_min = params.illumina_coverage_min ?: 40
  def ont_q_min = params.ont_q_min ?: 15
  def busco_min = params.busco_completeness_of_genome ?: 95
  def ont_breadth_min = params.ont_coverage_breadth_min ?: 98
  """
  #!/usr/bin/env python3
  import csv
  import glob
  import os
  import sys
  import re
  from pathlib import Path

  # QC thresholds passed from Nextflow params
  ONT_COVERAGE_MIN = ${ont_cov_min}
  ILLUMINA_COVERAGE_MIN = ${ill_cov_min}
  ONT_Q_MIN = ${ont_q_min}
  BUSCO_COMPLETENESS_MIN = ${busco_min}
  ONT_COVERAGE_BREADTH_MIN = ${ont_breadth_min}

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

  def rename_pypolca_columns(row):
      \"\"\"Rename pypolca columns to WAPHL standard names\"\"\"
      # Find and rename pypolca columns for all assemblers
      renamed_row = {}
      for key, value in row.items():
          # Match pattern: <assembler>_pypolca_<metric>
          if '_pypolca_input_contigs' in key:
              new_key = key.replace('_pypolca_input_contigs', '_num_contigs_before_pypolca')
              renamed_row[new_key] = value
          elif '_pypolca_output_contigs' in key:
              new_key = key.replace('_pypolca_output_contigs', '_num_contigs_after_pypolca')
              renamed_row[new_key] = value
          elif '_pypolca_perc_warnings' in key:
              new_key = key.replace('_pypolca_perc_warnings', '_pypolca_perc_warnings')
              renamed_row[new_key] = value
          else:
              renamed_row[key] = value
      return renamed_row

  def extract_busco_completeness(busco_string):
      \"\"\"Extract completeness percentage from BUSCO string\"\"\"
      # BUSCO format: C:98.5%[S:98.3%,D:0.2%],F:0.8%,M:0.7%,n:2472
      if not busco_string or busco_string == '':
          return None
      match = re.search(r'C:(\\d+\\.?\\d*)%', busco_string)
      if match:
          return float(match.group(1))
      return None

  def safe_float(value, default=0.0):
      \"\"\"Safely convert value to float\"\"\"
      if value is None or value == '' or value == 'unknown':
          return default
      try:
          return float(value)
      except (ValueError, TypeError):
          return default

  def perform_qc_checks(row):
      \"\"\"Perform QC checks and return qc_check status and reason_for_failure\"\"\"
      failures = []
      flags = []
      
      # Determine which assembler was used (check for populated columns)
      assemblers = []
      for assembler in ['flye', 'raven', 'myloasm', 'unicycler']:
          # Check if this assembler has data (look for busco or circulocov data)
          if any(key.startswith(f'{assembler}_busco_') for key in row.keys() if row.get(key, '')):
              assemblers.append(assembler)
      
      # If no assembler detected, try to infer from any column
      if not assemblers:
          for key in row.keys():
              for assembler in ['flye', 'raven', 'myloasm', 'unicycler']:
                  if key.startswith(f'{assembler}_') and row.get(key, ''):
                      if assembler not in assemblers:
                          assemblers.append(assembler)
                      break
      
      # For each assembler, perform checks
      for assembler in assemblers:
          # FAIL checks - critical failures
          
          # 1. Check ONT coverage depth
          ont_coverage_key = f'{assembler}_circulocov_nanopore_meandepth'
          if ont_coverage_key in row:
              ont_coverage = safe_float(row.get(ont_coverage_key))
              if ont_coverage > 0 and ont_coverage < ONT_COVERAGE_MIN:
                  failures.append(f'{ont_coverage_key}<{ONT_COVERAGE_MIN}')
          
          # 2. Check Illumina coverage depth (if available)
          ill_coverage_key = f'{assembler}_circulocov_illumina_meandepth'
          if ill_coverage_key in row and row.get(ill_coverage_key, ''):
              ill_coverage = safe_float(row.get(ill_coverage_key))
              if ill_coverage > 0 and ill_coverage < ILLUMINA_COVERAGE_MIN:
                  failures.append(f'{ill_coverage_key}<{ILLUMINA_COVERAGE_MIN}')
          
          # 3. Check ONT average quality
          ont_q_key = 'seqkit_AvgQual'
          if ont_q_key in row:
              ont_q = safe_float(row.get(ont_q_key))
              if ont_q > 0 and ont_q < ONT_Q_MIN:
                  failures.append(f'{ont_q_key}<{ONT_Q_MIN}')
          
          # 4. Check BUSCO completeness - find the most processed version
          # Priority: pypolca > polypolish > clair3 > reoriented > unicycler
          busco_completeness = None
          busco_key_used = None
          for step in ['pypolca', 'polypolish', 'clair3', 'reoriented', 'unicycler']:
              busco_key = f'{assembler}_busco_{step}'
              if busco_key in row and row.get(busco_key, ''):
                  busco_completeness = extract_busco_completeness(row[busco_key])
                  if busco_completeness is not None:
                      busco_key_used = busco_key
                      break
          
          if busco_completeness is not None and busco_completeness < BUSCO_COMPLETENESS_MIN:
              failures.append(f'{busco_key_used}<{BUSCO_COMPLETENESS_MIN}%')
          
          # 5. Check ONT coverage breadth
          ont_breadth_key = f'{assembler}_circulocov_nanopore_coverage_breadth_percent'
          # Try alternate key names
          if ont_breadth_key not in row:
              ont_breadth_key = 'nanopore_coverage_breadth_percent'
          if ont_breadth_key in row:
              ont_breadth = safe_float(row.get(ont_breadth_key))
              if ont_breadth > 0 and ont_breadth < ONT_COVERAGE_BREADTH_MIN:
                  failures.append(f'{ont_breadth_key}<{ONT_COVERAGE_BREADTH_MIN}%')
      
      # 6. Check human contamination (applies to all assemblers)
      human_contam_key = 'human_contamination'
      if human_contam_key in row:
          human_contam = safe_float(row.get(human_contam_key), default=1.0)
          if human_contam > 0 and human_contam < 0.95:
              failures.append(f'{human_contam_key}<0.95')
      
      # If we have failures, return FAIL
      if failures:
          return 'FAIL', '; '.join(failures)
      
      # FLAG checks - warnings
      
      # 1. Check match_containment_ani
      ani_key = 'match_containment_ani'
      if ani_key in row:
          ani = safe_float(row.get(ani_key), default=1.0)
          if ani > 0 and ani < 0.95:
              flags.append(f'{ani_key}<0.95')
      
      # 2. Check pypolca percent warnings (for each assembler)
      for assembler in assemblers:
          pypolca_warn_key = f'{assembler}_pypolca_perc_warnings'
          if pypolca_warn_key in row and row.get(pypolca_warn_key, ''):
              pypolca_warn = safe_float(row.get(pypolca_warn_key))
              if pypolca_warn >= 0.1:  # 0.1 = 10%
                  flags.append(f'{pypolca_warn_key}>=10%')
          
          # 3. Check if contig counts changed after pypolca
          before_key = f'{assembler}_num_contigs_before_pypolca'
          after_key = f'{assembler}_num_contigs_after_pypolca'
          if before_key in row and after_key in row:
              if row.get(before_key, '') and row.get(after_key, ''):
                  before_count = safe_float(row.get(before_key), default=-1)
                  after_count = safe_float(row.get(after_key), default=-1)
                  if before_count > 0 and after_count > 0 and before_count != after_count:
                      flags.append(f'{assembler}_contig_count_changed')
      
      # If we have flags, return FLAG
      if flags:
          return 'FLAG', '; '.join(flags)
      
      # Otherwise, PASS
      return 'PASS', ''

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
      # Extract the top sourmash taxa hit and match_containment_ani for a sample
      # New format: sample, gather_rank, taxa, match_containment_ani, f_unique_to_query, average_abund
      for mash_file in mash_files:
          # Extract sample name from file path
          file_sample = os.path.basename(mash_file).replace('_taxa.txt', '')
          if sample == file_sample:
              try:
                  with open(mash_file, 'r') as f:
                      lines = f.readlines()
                      if len(lines) > 1:  # Skip header
                          # Get the first data line (best match - gather_rank 1 or lowest rank)
                          data_line = lines[1].strip().split('\\t')
                          if len(data_line) >= 4:
                              taxa = data_line[2]                  # taxa column (cleaned name)
                              match_ani = data_line[3]            # match_containment_ani column
                              # Handle "No matches" case
                              if taxa == "No matches >= 95% ANI":
                                  return 'No matches >= 95% ANI', '0'
                              return taxa, match_ani
              except Exception as e:
                  print("Error reading {}: {}".format(mash_file, e))
                  continue
      return 'unknown', 'unknown'


  def get_coverage_data(sample, coverage_files):
      # Extract coverage data for a sample (ONT or Illumina)
      coverage_data = {}
      
      for coverage_file in coverage_files:
          # Extract sample name from file path
          file_sample = os.path.basename(coverage_file).replace('_coverage_summary.tsv', '')
          if sample == file_sample:
              try:
                  with open(coverage_file, 'r') as f:
                      lines = f.readlines()
                      read_platform = None
                      
                      # Parse the main coverage statistics
                      if len(lines) > 1:  # Skip header
                          # Get the first data line
                          data_line = lines[1].strip().split('\\t')
                          if len(data_line) >= 6:
                              read_platform = data_line[1]
                              coverage_data.update({
                                  f'{read_platform.lower()}_total_bases': data_line[2],
                                  f'{read_platform.lower()}_covered_bases': data_line[3], 
                                  f'{read_platform.lower()}_coverage_breadth_percent': data_line[4],
                                  f'{read_platform.lower()}_max_depth': data_line[5],
                                  'read_platform': read_platform
                              })
                      
                      # If we found data, return it
                      if coverage_data:
                          return coverage_data
                          
              except Exception as e:
                  print("Error reading {}: {}".format(coverage_file, e))
                  continue
      return coverage_data

  def get_human_contamination_data(sample, human_contamination_files):
      # Extract human contamination data for a sample - return sourmash containment instead of contamination level
      for human_file in human_contamination_files:
          # Extract sample name from file path
          file_sample = os.path.basename(human_file).replace('_human_summary.txt', '')
          if sample == file_sample:
              try:
                  with open(human_file, 'r') as f:
                      line = f.readline().strip()
                      if line:
                          parts = line.split('\\t')
                          if len(parts) >= 2:
                              # Return the sourmash containment instead of contamination level
                              return parts[1]  # containment column
              except Exception as e:
                  print("Error reading {}: {}".format(human_file, e))
                  continue
      return 'unknown'

  def get_consensus_filepath(sample, consensus_files, outdir):
      # Get the consensus genome filepath for a sample
      # Prioritize the most processed version available
      # Return the full path where the file will be published
      
      # First priority: pypolca (final polishing step)
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_pypolca.fasta' in consensus_file:
              return os.path.join(outdir, sample, 'pypolca', os.path.basename(consensus_file))
      
      # Second priority: polypolish
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_polypolish.fasta' in consensus_file:
              return os.path.join(outdir, sample, 'polypolish', os.path.basename(consensus_file))
      
      # Third priority: clair3
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_clair3.fasta' in consensus_file:
              return os.path.join(outdir, sample, 'clair3', os.path.basename(consensus_file))
      
      # Fourth priority: reoriented
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_reoriented.fasta' in consensus_file:
              return os.path.join(outdir, sample, 'dnaapler', os.path.basename(consensus_file))
      
      # Fifth priority: unicycler (hybrid assembly)
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_unicycler.fasta' in consensus_file:
              return os.path.join(outdir, sample, 'unicycler', os.path.basename(consensus_file))
      
      # Fallback: any fasta file for the sample
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and consensus_file.endswith('.fasta'):
              return os.path.join(outdir, sample, os.path.basename(consensus_file))
      
      return ''

  def get_sub_fasta_filepath(sample, consensus_files, outdir):
      # Get the sub_fasta filepath if it exists
      # Return the full path where the file will be published
      
      # First priority: Look for files with 'sub_' prefix (circular assemblies with special headers)
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).replace('sub_', '').split('_')[0]
          if sample == file_sample and 'sub_' in consensus_file and consensus_file.endswith('.fasta'):
              return os.path.join(outdir, sample, 'consensus', os.path.basename(consensus_file))
      
      # Second priority: Look for reoriented files as they are submission-ready
      for consensus_file in consensus_files:
          file_sample = os.path.basename(consensus_file).split('_')[0]
          if sample == file_sample and '_reoriented' in consensus_file and consensus_file.endswith('.fasta'):
              return os.path.join(outdir, sample, 'dnaapler', os.path.basename(consensus_file))
      
      # If neither sub_ files nor reoriented files are found, return empty string
      return ''

  def main():
      try:
          # Read the original donut falls summary
          summary_data = read_tsv_file('${donut_falls_summary}')
          
          # Get all sourmash taxa files
          mash_files = glob.glob('sourmash_taxa/*_taxa.txt')
          
          # Get all consensus files
          consensus_files = glob.glob('consensus/*.fasta')

          # Get all coverage analysis files
          coverage_files = glob.glob('coverage/*_coverage_summary.tsv')

          # Get all human contamination files
          human_contamination_files = glob.glob('human_contamination/*_human_summary.txt')
          
          print("Found {} samples in summary".format(len(summary_data)))
          print("Found {} sourmash taxa files: {}".format(len(mash_files), mash_files))
          print("Found {} consensus files: {}".format(len(consensus_files), consensus_files))
          print("Found {} coverage analysis files: {}".format(len(coverage_files), coverage_files))
          print("Found {} human contamination files: {}".format(len(human_contamination_files), human_contamination_files))
          
          # Process each sample (enumerate to modify in place)
          for idx, row in enumerate(summary_data):
              sample = row['sample']
              print("Processing sample: {}".format(sample))
              
              # Rename pypolca columns first
              row = rename_pypolca_columns(row)
              
              # Get sourmash taxa and match containment ANI
              sourmash_taxa, match_ani = get_mash_taxa_and_distance(sample, mash_files)
              row['sourmash_taxa'] = sourmash_taxa
              row['match_containment_ani'] = match_ani

              # Get coverage data (ONT or Illumina)
              coverage_data = get_coverage_data(sample, coverage_files)
              for key, value in coverage_data.items():
                  row[key] = value

              # Get human contamination data (MASH distance)
              human_contamination = get_human_contamination_data(sample, human_contamination_files)
              row['human_contamination'] = human_contamination
              
              # Determine assembly type
              row['assembly_type'] = determine_assembly_type(row)
              
              # Get consensus filepath
              row['consensus_filepath'] = get_consensus_filepath(sample, consensus_files, '${params.outdir}')
              
              # Get sub_fasta filepath
              row['reoriented_assembly_filepath'] = get_sub_fasta_filepath(sample, consensus_files, '${params.outdir}')
              
              # Perform QC checks and add qc_check and reason_for_failure columns
              qc_status, qc_reason = perform_qc_checks(row)
              row['qc_check'] = qc_status
              row['reason_for_failure'] = qc_reason
              
              # Update the row in the list
              summary_data[idx] = row
          
          # Get all fieldnames (original + new columns)
          # Ensure qc_check and reason_for_failure come right after sample
          if summary_data:
              all_keys = list(summary_data[0].keys())
              # Reorder to put qc_check and reason_for_failure after sample
              fieldnames = ['sample']
              if 'qc_check' in all_keys:
                  fieldnames.append('qc_check')
              if 'reason_for_failure' in all_keys:
                  fieldnames.append('reason_for_failure')
              # Add remaining keys
              for key in all_keys:
                  if key not in fieldnames:
                      fieldnames.append(key)
          else:
              fieldnames = ['sample', 'qc_check', 'reason_for_failure']
          
          # Write the enhanced summary
          write_tsv_file('waphl_final_summary.tsv', summary_data, fieldnames)
          
          print("Enhanced summary created with {} samples".format(len(summary_data)))
          print("Found {} mash taxa files".format(len(mash_files)))
          print("Found {} consensus files".format(len(consensus_files)))
          print("Found {} coverage analysis files".format(len(coverage_files)))
          print("Found {} human contamination files".format(len(human_contamination_files)))
          
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