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
      # Extract genus and species from RefSeq reference names
      ref_name = \$1
      
      # Remove path if present
      gsub(/.*\\//, "", ref_name)
      
      # Remove file extensions (.fna, .fasta, .fa, .genomic, etc.)
      gsub(/\\.(fna|fasta|fa).*\$/, "", ref_name)
      gsub(/_genomic\$/, "", ref_name)
      
      # Initialize genus_species
      genus_species = "Unknown"
      
      # Handle RefSeq GCF format: GCF_000005825.2_ASM582v2 -> extract from assembly name after GCF_XXXXXX.X_
      if (match(ref_name, /^GCF_[0-9]+\\.[0-9]+_(.+)\$/, arr)) {
        assembly_part = arr[1]
        # Many assembly names start with genus_species pattern
        if (match(assembly_part, /^([A-Z][a-z]+)_([a-z]+)/, arr2)) {
          genus_species = arr2[1] " " arr2[2]
        } else {
          # Some assembly names are more complex, try to find recognizable patterns
          # Look for two consecutive words that look like genus species
          split(assembly_part, parts, "_")
          for (i = 1; i <= length(parts)-1; i++) {
            if (match(parts[i], /^[A-Z][a-z]+\$/) && match(parts[i+1], /^[a-z]+\$/) && length(parts[i+1]) > 2) {
              genus_species = parts[i] " " parts[i+1]
              break
            }
          }
          # If still no match, use the assembly name with underscores replaced by spaces
          if (genus_species == "Unknown") {
            genus_species = assembly_part
            gsub(/_/, " ", genus_species)
            # Limit to first two words to avoid long strings
            if (split(genus_species, words, " ") >= 2) {
              genus_species = words[1] " " words[2]
            }
          }
        }
      }
      # Handle direct genus_species format (common in older RefSeq entries)
      else if (match(ref_name, /^([A-Z][a-z]+)_([a-z]+)/, arr)) {
        genus_species = arr[1] " " arr[2]
      }
      # Handle other formats - try to find genus species pattern
      else {
        split(ref_name, parts, "_")
        for (i = 1; i <= length(parts)-1; i++) {
          if (match(parts[i], /^[A-Z][a-z]+\$/) && match(parts[i+1], /^[a-z]+\$/) && length(parts[i+1]) > 2) {
            genus_species = parts[i] " " parts[i+1]
            break
          }
        }
      }
      
      # Clean up genus_species - remove any remaining problematic characters
      gsub(/[^a-zA-Z0-9 -]/, "", genus_species)
      
      # If we still have Unknown, try one more approach - look for any capitalized word followed by lowercase
      if (genus_species == "Unknown" || genus_species == "") {
        if (match(ref_name, /([A-Z][a-zA-Z]*)[_\\s]+([a-z][a-zA-Z]*)/, arr)) {
          genus_species = arr[1] " " arr[2]
        } else {
          genus_species = "Unknown organism"
        }
      }
      
      print sample, \$1, \$3, \$4, \$5, genus_species
    }' > mash_taxa/${prefix}_taxa.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}