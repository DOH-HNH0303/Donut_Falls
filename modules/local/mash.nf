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
  def db     = task.ext.db     ?: '/db/RefSeqSketchesDefaults.msh'
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
      ref_name = \$1
      genus_species = "Unknown"
      
      # Remove path if present
      gsub(/.*\\//, "", ref_name)
      
      # Remove file extensions
      gsub(/\\.(fna|fasta|fa).*\$/, "", ref_name)
      gsub(/_genomic\$/, "", ref_name)
      
      # Try to extract genus species from different RefSeq formats
      # Format 1: GCF_000005825.2_ASM582v2 -> look for genus_species pattern after GCF_
      if (match(ref_name, /^GCF_[0-9]+\\.[0-9]+_/)) {
        # Remove the GCF prefix
        gsub(/^GCF_[0-9]+\\.[0-9]+_/, "", ref_name)
        
        # Look for genus_species pattern at start
        if (match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
          split(ref_name, parts, "_")
          if (length(parts) >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
            genus_species = parts[1] " " parts[2]
          }
        }
      }
      
      # Format 2: Direct genus_species format
      if (genus_species == "Unknown" && match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
        split(ref_name, parts, "_")
        if (length(parts) >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
          genus_species = parts[1] " " parts[2]
        }
      }
      
      # Format 3: Look for any genus species pattern in the name
      if (genus_species == "Unknown") {
        split(ref_name, parts, "_")
        for (i = 1; i <= length(parts)-1; i++) {
          if (match(parts[i], /^[A-Z][a-z]+\$/) && match(parts[i+1], /^[a-z]+\$/) && length(parts[i+1]) > 2) {
            genus_species = parts[i] " " parts[i+1]
            break
          }
        }
      }
      
      # Clean up and fallback
      if (genus_species == "Unknown" || genus_species == "") {
        genus_species = "Unknown organism"
      }
      
      # Remove any problematic characters
      gsub(/[^a-zA-Z0-9 -]/, "", genus_species)
      
      print sample, \$1, \$3, \$4, \$5, genus_species
    }' > mash_taxa/${prefix}_taxa.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}