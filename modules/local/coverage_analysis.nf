process COVERAGE_ANALYSIS {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "coverage_analysis/*"
  container     'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
  time          '30m'

  input:
  tuple val(meta), file(fasta), file(reads)

  output:
  tuple val(meta), file("coverage_analysis/*_coverage_summary.tsv"), emit: summary
  path "coverage_analysis/*", emit: everything
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  
  // Determine read type based on file structure
  def is_paired = reads instanceof List && reads.size() == 2
  def is_ont = !is_paired && reads.toString().endsWith('.fastq.gz')
  def read_type = is_paired ? "illumina" : (is_ont ? "ont" : "unknown")
  
  """
  mkdir -p coverage_analysis

  # Index the reference genome
  samtools faidx ${fasta}

  # Determine read type and map accordingly
  if [ "${read_type}" == "illumina" ]; then
    echo "Processing Illumina paired-end reads"
    
    # Map Illumina reads to the assembly using BWA
    bwa index ${fasta}
    bwa mem -t ${task.cpus} ${fasta} ${reads[0]} ${reads[1]} | \\
      samtools sort -@ ${task.cpus} -o coverage_analysis/${prefix}_mapped.bam -
    
    # Calculate read statistics for Illumina
    total_bases_reads=\$(zcat ${reads[0]} ${reads[1]} | awk 'NR%4==2 {sum+=length(\$0)} END {print sum}')
    read_platform="Illumina"
    
  elif [ "${read_type}" == "ont" ]; then
    echo "Processing ONT long reads"
    
    # Map ONT reads to the assembly using minimap2
    minimap2 -ax map-ont -t ${task.cpus} ${fasta} ${reads} | \\
      samtools sort -@ ${task.cpus} -o coverage_analysis/${prefix}_mapped.bam -
    
    # Calculate read statistics for ONT
    total_bases_reads=\$(zcat ${reads} | awk 'NR%4==2 {sum+=length(\$0)} END {print sum}')
    read_platform="ONT"
    
  else
    echo "Error: Unable to determine read type"
    exit 1
  fi

  # Index the BAM file
  samtools index coverage_analysis/${prefix}_mapped.bam

  # Calculate coverage statistics
  samtools depth coverage_analysis/${prefix}_mapped.bam > coverage_analysis/${prefix}_depth.txt

  # Calculate basic mapping statistics
  samtools flagstat coverage_analysis/${prefix}_mapped.bam > coverage_analysis/${prefix}_flagstat.txt

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
    print "sample\\tread_platform\\ttotal_bases\\tcovered_bases\\tcoverage_breadth_percent\\tmean_depth\\tmax_depth"
    print "${prefix}\\t" read_platform "\\t" total_bases "\\t" covered_bases "\\t" coverage_breadth "\\t" mean_depth "\\t" max_depth
  }' read_platform="\$read_platform" coverage_analysis/${prefix}_depth.txt > coverage_analysis/${prefix}_coverage_summary.tsv

  # Get genome size from fasta index
  genome_size=\$(awk '{sum+=\$2} END {print sum}' ${fasta}.fai)

  # Calculate theoretical coverage from read data
  theoretical_coverage=\$(awk -v reads=\$total_bases_reads -v genome=\$genome_size 'BEGIN {printf "%.2f", reads/genome}')

  # Add theoretical coverage to summary
  echo -e "\\n# Theoretical coverage based on read data" >> coverage_analysis/${prefix}_coverage_summary.tsv
  echo -e "genome_size\\ttotal_read_bases\\ttheoretical_coverage\\tread_platform" >> coverage_analysis/${prefix}_coverage_summary.tsv
  echo -e "\$genome_size\\t\$total_bases_reads\\t\$theoretical_coverage\\t\$read_platform" >> coverage_analysis/${prefix}_coverage_summary.tsv

  # Create a coverage histogram
  awk '{print \$3}' coverage_analysis/${prefix}_depth.txt | \\
    sort -n | \\
    uniq -c | \\
    awk '{print \$2 "\\t" \$1}' > coverage_analysis/${prefix}_coverage_histogram.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    minimap2: \$(minimap2 --version 2>/dev/null || echo "not available")
    bwa: \$(bwa 2>&1 | grep -i version | awk '{print \$NF}' || echo "not available")
  END_VERSIONS
  """
}