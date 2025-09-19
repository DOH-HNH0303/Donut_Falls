process MASH_TAXA {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "mash_taxa/*"
  container     'staphb/mash:2.3'
  time          '30m'
  
  // Note: Memory requirements depend on RefSeq database size
  // Typically requires 8-16 GB RAM for RefSeq database operations
  
  input:
  tuple val(meta), file(fasta)

  output:
  tuple val(meta), file("mash_taxa/*_taxa.txt"), emit: taxa
  path "mash_taxa/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def db     = task.ext.db     ?: '/mash_db/RefSeqSketchesDefaults.msh'
  """
  mkdir -p mash_taxa

  # Sketch the assembly
  mash sketch ${args} -o mash_taxa/${prefix} ${fasta}

  # Find closest matches in reference database
  mash dist ${db} mash_taxa/${prefix}.msh | \
    sort -gk3 | \
    head -10 | \
    awk -v sample=${prefix} 'BEGIN{OFS="\\t"; print "sample", "reference", "distance", "p_value", "shared_hashes", "genus_species"} 
    {
      # Extract genus and species from reference name
      split(\$1, ref_parts, "_")
      genus_species = ref_parts[1] "_" ref_parts[2]
      gsub(/.*\\//, "", genus_species)
      gsub(/\\.fna.*/, "", genus_species)
      print sample, \$1, \$3, \$4, \$5, genus_species
    }' > mash_taxa/${prefix}_taxa.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}