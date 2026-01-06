process pypolca {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/pypolca:0.4.0'
  time          '30m'
  
  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file("pypolca/*_pypolca.fasta"), optional: true, emit: fasta
  path "pypolca/*pypolca_summary.tsv", optional: true, emit: summary
  path "pypolca/*_pypolca_failed_contigs.fasta", optional: true, emit: failed_contigs
  path "pypolca/*/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--careful'
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_polypolish','')}"
  """
  mkdir -p pypolca/${prefix}
  mkdir -p pypolca_temp
  
  # Replace spaces in contig names to avoid issues with pypolca
  sed "s/ /_.._/g" ${fasta} > input.fasta
  
  # Split the input fasta into individual contigs
  awk '/^>/ {if(filename) close(filename); filename="pypolca_temp/contig_"++count".fasta"} {print > filename}' input.fasta
  
  # Process each contig with pypolca
  for contig_file in pypolca_temp/contig_*.fasta; do
    contig_name=\$(basename \${contig_file} .fasta)
    
    pypolca run ${args} \
      -a \${contig_file} \
      -1 ${fastq[0]} \
      -2 ${fastq[1]} \
      -t ${task.cpus} \
      -o pypolca/${prefix}/\${contig_name}
    
    # Collect the corrected contig if it exists, otherwise mark as failed
    if [ -f "pypolca/${prefix}/\${contig_name}/pypolca_corrected.fasta" ]; then
      cat pypolca/${prefix}/\${contig_name}/pypolca_corrected.fasta >> pypolca/${prefix}_pypolca.fasta
    else
      # If pypolca didn't produce output, add to failed contigs file
      cat \${contig_file} >> pypolca/${prefix}_pypolca_failed_contigs.fasta
    fi
  done
  
  # Restore original spacing in contig names for successful polishes
  if [ -f "pypolca/${prefix}_pypolca.fasta" ]; then
    sed -i "s/_.._/ /g" pypolca/${prefix}_pypolca.fasta
  fi
  
  # Restore original spacing in contig names for failed polishes
  if [ -f "pypolca/${prefix}_pypolca_failed_contigs.fasta" ]; then
    sed -i "s/_.._/ /g" pypolca/${prefix}_pypolca_failed_contigs.fasta
  fi
  
  # Combine all pypolca reports into a summary
  first_report=true
  for report in pypolca/${prefix}/*/pypolca.report; do
    if [ -f "\${report}" ]; then
      if [ "\${first_report}" = true ]; then
        # Create header from first report
        cut -f 1 -d : \${report} | \
          sed 's/ /_/g' | \
          tr "\\n" "\\t" | \
          awk '{print "sample\\t" \$0 }' \
          > pypolca/${prefix}_pypolca_summary.tsv
        first_report=false
      fi
      
      # Add data row for this contig
      cut -f 2 -d : \${report} | \
        awk '{( \$1 = \$1 ) ; print \$0 }' | \
        sed 's/ /_/g' | \
        tr "\\n" "\\t" | \
        awk '{print "${prefix}\\t" \$0 }' \
        >> pypolca/${prefix}_pypolca_summary.tsv
    fi
  done

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    pypolca: \$(pypolca --version | head -n 1 | awk '{print \$NF}')
  END_VERSIONS
  """
}