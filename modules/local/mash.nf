process MASH_TAXA {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "mash_taxa/*"   
    container     'staphb/mash:2.3'
    
    input:
    tuple val(meta), path(fasta)
    
    when:
    fasta_file != null && !fasta_file.name.contains('reoriented')

    output:
    tuple val(meta), path("mash_taxa/*.txt"), emit: taxa
    path "versions.yml", emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p mash_taxa
    
    # Sketch the assembly
    mash sketch -o mash_taxa/${prefix} ${fasta}
    
    # Find closest matches in reference database
    # Fixed AWK script that works with different AWK implementations
    mash dist /db/RefSeqSketchesDefaults.msh mash_taxa/${prefix}.msh | \\
        sort -gk3 | \\
        head -10 | \\
        awk -v sample=${prefix} '
        BEGIN{
            OFS="\\t"
            print "sample", "reference", "distance", "p_value", "shared_hashes", "genus_species"
        } 
        {
            # Extract genus and species from reference name
            ref_name = \$1
            genus_species = "Unknown"

            # Remove path if present
            gsub(/.*\\//, "", ref_name)

            # Remove file extensions
            gsub(/\\.(fna|fasta|fa).*\$/, "", ref_name)
            gsub(/_genomic\$/, "", ref_name)

            # Handle complex RefSeq formats with embedded genus_species
            if (match(ref_name, /[A-Z][a-z]+_[a-z][a-z]+/)) {
                matched_part = substr(ref_name, RSTART, RLENGTH)
                
                # Use split with return value check
                n = split(matched_part, genus_parts, "_")
                if (n >= 2 && match(genus_parts[1], /^[A-Z][a-z]+\$/) && match(genus_parts[2], /^[a-z]+\$/) && length(genus_parts[2]) > 2) {
                    genus_species = genus_parts[1] " " genus_parts[2]
                }
                delete genus_parts
            }

            # If still unknown, try standard GCF format
            if (genus_species == "Unknown" && match(ref_name, /GCF_[0-9]+\\.[0-9]+/)) {
                gsub(/.*GCF_[0-9]+\\.[0-9]+[^A-Za-z]*/, "", ref_name)
                
                if (match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
                    n = split(ref_name, parts, "_")
                    if (n >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
                        genus_species = parts[1] " " parts[2]
                    }
                    delete parts
                }
            }

            # If still unknown, try direct genus_species format
            if (genus_species == "Unknown" && match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
                n = split(ref_name, parts, "_")
                if (n >= 2 && match(parts[1], /^[A-Z][a-z]+\$/) && match(parts[2], /^[a-z]+\$/) && length(parts[2]) > 2) {
                    genus_species = parts[1] " " parts[2]
                }
                delete parts
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
    "WAPHL_ANALYSIS:MASH_TAXA":
        mash: \$( mash --version )
    END_VERSIONS
    """
}