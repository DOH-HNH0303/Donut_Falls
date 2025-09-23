process ONT_COVERAGE {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "ont_coverage/*"
  container     'staphb/circulocov:0.1.20240104'
  time          '30m'

  input:
  tuple val(meta), file(fasta), file(nanopore)

  output:
  tuple val(meta), file("ont_coverage/*_ont_coverage_summary.tsv"), emit: summary
  path "ont_coverage/*", emit: everything
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p ont_coverage

  # Index the reference genome
  samtools faidx ${fasta}

  # Map ONT reads to the assembly
  minimap2 -ax map-ont -t ${task.cpus} ${fasta} ${nanopore} | \\
    samtools sort -@ ${task.cpus} -o ont_coverage/${prefix}_ont_mapped.bam -

  # Index the BAM file
  samtools index ont_coverage/${prefix}_ont_mapped.bam

  # Calculate coverage statistics
  samtools depth ont_coverage/${prefix}_ont_mapped.bam > ont_coverage/${prefix}_depth.txt

  # Calculate basic mapping statistics
  samtools flagstat ont_coverage/${prefix}_ont_mapped.bam > ont_coverage/${prefix}_flagstat.txt

  # Calculate coverage summary using awk
  awk '
  BEGIN {
    total_bases = 0
    covered_bases = 0
    sum_depth = 0
    max_depth = 0
  }
  {
    total_bases++
    if (\$3 > 0) {
      covered_bases++
      sum_depth += \$3
      if (\$3 > max_depth) max_depth = \$3
    }
  }
  END {
    mean_depth = (covered_bases > 0) ? sum_depth / covered_bases : 0
    coverage_breadth = (total_bases > 0) ? (covered_bases / total_bases) * 100 : 0
    print "sample\\ttotal_bases\\tcovered_bases\\tcoverage_breadth_percent\\tmean_depth\\tmax_depth"
    print "${prefix}\\t" total_bases "\\t" covered_bases "\\t" coverage_breadth "\\t" mean_depth "\\t" max_depth
  }' ont_coverage/${prefix}_depth.txt > ont_coverage/${prefix}_ont_coverage_summary.tsv

  # Get genome size from fasta index
  genome_size=\$(awk '{sum+=\$2} END {print sum}' ${fasta}.fai)

  # Calculate theoretical coverage from read data
  total_bases_reads=\$(zcat ${nanopore} | awk 'NR%4==2 {sum+=length(\$0)} END {print sum}')
  theoretical_coverage=\$(echo "scale=2; \$total_bases_reads / \$genome_size" | bc -l)

  # Add theoretical coverage to summary
  echo -e "\\n# Theoretical coverage based on read data" >> ont_coverage/${prefix}_ont_coverage_summary.tsv
  echo -e "genome_size\\ttotal_read_bases\\ttheoretical_coverage" >> ont_coverage/${prefix}_ont_coverage_summary.tsv
  echo -e "\$genome_size\\t\$total_bases_reads\\t\$theoretical_coverage" >> ont_coverage/${prefix}_ont_coverage_summary.tsv

  # Create a coverage histogram
  awk '{print \$3}' ont_coverage/${prefix}_depth.txt | \\
    sort -n | \\
    uniq -c | \\
    awk '{print \$2 "\\t" \$1}' > ont_coverage/${prefix}_coverage_histogram.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    minimap2: \$(minimap2 --version)
  END_VERSIONS
  """
}