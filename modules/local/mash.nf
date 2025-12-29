process MASH_TAXA {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "mash_taxa/*"   
    container     'staphb/mash:2.3'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("mash_taxa/*.txt"), emit: taxa
    path "versions.yml", emit: versions
    
    when:
    fasta != null && !fasta.name.contains('reoriented')
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def database = params.mash_db ?: ''
    """
    mkdir -p mash_taxa
    
    # Sketch the assembly
    mash sketch -o mash_taxa/${prefix} ${fasta}
    
    # Check if database is provided
    if [ -z "${database}" ]; then
        echo "WARNING: No Mash database specified!" >&2
        echo "Provide --mash_db to enable taxonomic classification" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tdistance\\tp_value\\tshared_hashes\\tgenus_species" > mash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tno_database_provided\\t1.0\\t1.0\\t0/1000\\tUnknown" >> mash_taxa/${prefix}_taxa.txt
        exit 0
    fi
    
    # Check if database file exists
    if [ ! -f "${database}" ]; then
        echo "WARNING: Mash database not found at: ${database}" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tdistance\\tp_value\\tshared_hashes\\tgenus_species" > mash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tdatabase_not_found\\t1.0\\t1.0\\t0/1000\\tUnknown" >> mash_taxa/${prefix}_taxa.txt
        exit 0
    fi
    
    # Find closest matches in reference database
    # Fixed AWK script that works with different AWK implementations
    mash dist ${database} mash_taxa/${prefix}.msh | \\
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

process MASH_HUMAN_CONTAMINATION {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "mash_human/*"   
    container     'staphb/mash:2.3'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("mash_human/*_human_contamination.txt"), emit: human_contamination
    tuple val(meta), path("mash_human/*_human_summary.txt"), emit: human_summary
    path "versions.yml", emit: versions
    
    when:
    fasta != null
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def database_human = params.mash_db_human ?: params.mash_db ?: ''
    """
    mkdir -p mash_human
    
    # Sketch the assembly
    mash sketch -o mash_human/${prefix} ${fasta}
    
    # Check if database is provided
    if [ -z "${database_human}" ]; then
        echo "WARNING: No Mash database specified for human contamination check!" >&2
        echo "Provide --mash_db_human (or --mash_db) to enable contamination detection" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tdistance\\tp_value\\tshared_hashes\\tcontamination_level" > mash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\tno_database_provided\\t1.0\\t1.0\\t0/1000\\tunknown" >> mash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\t1.0\\tunknown" > mash_human/${prefix}_human_summary.txt
        exit 0
    fi
    
    # Check if database file exists
    if [ ! -f "${database_human}" ]; then
        echo "WARNING: Mash database not found at: ${database_human}" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tdistance\\tp_value\\tshared_hashes\\tcontamination_level" > mash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\tdatabase_not_found\\t1.0\\t1.0\\t0/1000\\tunknown" >> mash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\t1.0\\tunknown" > mash_human/${prefix}_human_summary.txt
        exit 0
    fi
    
    # Check against database for human sequences
    mash dist ${database_human} mash_human/${prefix}.msh | \\
        grep -i "homo.*sapiens\\|human\\|GCF_000001405" | \\
        sort -gk3 | \\
        head -5 | \\
        awk -v sample=${prefix} '
        BEGIN{
            OFS="\\t"
            print "sample", "human_reference", "distance", "p_value", "shared_hashes", "contamination_level"
            min_distance = 1.0
            contamination_status = "none"
        } 
        {
            distance = \$3
            if (distance < min_distance) {
                min_distance = distance
            }
            
            # Determine contamination level based on distance
            if (distance < 0.01) {
                contamination_level = "high"
            } else if (distance < 0.05) {
                contamination_level = "moderate" 
            } else if (distance < 0.1) {
                contamination_level = "low"
            } else {
                contamination_level = "minimal"
            }
            
            print sample, \$1, distance, \$4, \$5, contamination_level
        }
        END {
            # If no human sequences found, create a default entry
            if (NR == 0) {
                print sample, "no_human_reference_found", "1.0", "1.0", "0/1000", "none"
            }
        }' > mash_human/${prefix}_human_contamination.txt

    # Create a summary line with overall contamination assessment
    awk -v sample=${prefix} '
    BEGIN {
        OFS="\\t"
        min_distance = 1.0
        overall_contamination = "none"
    }
    NR > 1 {  # Skip header
        if (\$3 < min_distance) {
            min_distance = \$3
            overall_contamination = \$6
        }
    }
    END {
        print sample, min_distance, overall_contamination > "mash_human/" sample "_human_summary.txt"
    }' mash_human/${prefix}_human_contamination.txt

    cat <<-END_VERSIONS > versions.yml
    "WAPHL_ANALYSIS:MASH_HUMAN_CONTAMINATION":
        mash: \$( mash --version )
    END_VERSIONS
    """
}