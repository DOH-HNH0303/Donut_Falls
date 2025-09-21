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
      
      # IMPROVED: Handle complex RefSeq formats with embedded genus_species
      # Look for genus_species pattern anywhere in the name (not just at start)
      # Pattern: [A-Z][a-z]+_[a-z]+ (genus_species with underscore)
      if (match(ref_name, /[A-Z][a-z]+_[a-z][a-z]+/)) {
        # Extract the matched portion
        matched_part = substr(ref_name, RSTART, RLENGTH)
        
        # Split the matched part and validate
        split(matched_part, genus_parts, "_")
        if (length(genus_parts) >= 2 && match(genus_parts[1], /^[A-Z][a-z]+\$/) && match(genus_parts[2], /^[a-z]+\$/) && length(genus_parts[2]) > 2) {
          genus_species = genus_parts[1] " " genus_parts[2]
        }
      }
      
      # If still unknown, try standard GCF format
      if (genus_species == "Unknown" && match(ref_name, /GCF_[0-9]+\\.[0-9]+/)) {
        # Remove everything up to and including GCF_XXXXXX.X
        gsub(/.*GCF_[0-9]+\\.[0-9]+[^A-Za-z]*/, "", ref_name)
        
        # Look for genus_species pattern at start of remaining string
        if (match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
          split(ref_name, parts, "_")
          if (length(parts) >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
            genus_species = parts[1] " " parts[2]
          }
        }
      }
      
      # If still unknown, try direct genus_species format
      if (genus_species == "Unknown" && match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
        split(ref_name, parts, "_")
        if (length(parts) >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
          genus_species = parts[1] " " parts[2]
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