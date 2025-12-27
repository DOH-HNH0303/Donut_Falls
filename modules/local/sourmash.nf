process SOURMASH_TAXA {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "sourmash_taxa/*"   
    container     'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("sourmash_taxa/*.txt"), emit: taxa
    path "versions.yml", emit: versions
    
    when:
    fasta != null && !fasta.name.contains('reoriented')
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ksize = task.ext.ksize ?: '31'
    def scaled = task.ext.scaled ?: '1000'
    """
    mkdir -p sourmash_taxa
    
    # Sketch the assembly with specified k-mer size and scaled value
    sourmash sketch dna \\
        -p k=${ksize},scaled=${scaled},abund \\
        --name ${prefix} \\
        -o sourmash_taxa/${prefix}.sig \\
        ${fasta}
    
    # Gather against reference database to find closest matches
    # Note: This assumes a sourmash database is available at /db/sourmash_db
    # You may need to adjust the database path based on your setup
    if [ -f "/db/gtdb-rs214.genomic.k31.zip" ]; then
        DATABASE="/db/gtdb-rs214.genomic.k31.zip"
    elif [ -f "/db/genbank-2022.03.genomic.k31.zip" ]; then
        DATABASE="/db/genbank-2022.03.genomic.k31.zip"
    else
        echo "Warning: No sourmash database found at expected locations" >&2
        echo "Skipping gather step - creating empty output" >&2
        echo -e "sample\\treference\\tcontainment\\tf_match\\tANI\\tgenus_species" > sourmash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tno_database_found\\t0\\t0\\t0\\tUnknown organism" >> sourmash_taxa/${prefix}_taxa.txt
    fi
    
    if [ -n "\${DATABASE}" ]; then
        # Run gather to find best matches
        sourmash gather \\
            sourmash_taxa/${prefix}.sig \\
            \${DATABASE} \\
            -o sourmash_taxa/${prefix}_gather.csv \\
            --threshold-bp 50000 \\
            -k ${ksize}
        
        # Parse gather results to create taxa output compatible with downstream processing
        # Extract top 10 matches and format similar to MASH output
        awk -F',' -v sample=${prefix} '
        BEGIN {
            OFS="\\t"
            print "sample", "reference", "containment", "f_match", "ANI", "genus_species"
        }
        NR > 1 && NR <= 11 {
            # Extract genus and species from name field (column 2)
            ref_name = \$2
            genus_species = "Unknown"
            
            # Remove quotes
            gsub(/"/, "", ref_name)
            
            # Try to extract genus_species pattern
            if (match(ref_name, /[A-Z][a-z]+_[a-z][a-z]+/)) {
                matched_part = substr(ref_name, RSTART, RLENGTH)
                n = split(matched_part, parts, "_")
                if (n >= 2 && length(parts[2]) > 2) {
                    genus_species = parts[1] " " parts[2]
                }
            }
            
            # Try GCF/GCA format if still unknown
            if (genus_species == "Unknown" && match(ref_name, /GC[FA]_[0-9]+/)) {
                # Extract everything after the accession
                sub(/.*GC[FA]_[0-9]+\\.[0-9]+[^A-Za-z]*/, "", ref_name)
                if (match(ref_name, /^[A-Z][a-z]+_[a-z]+/)) {
                    n = split(ref_name, parts, "_")
                    if (n >= 2 && length(parts[2]) > 2) {
                        genus_species = parts[1] " " parts[2]
                    }
                }
            }
            
            # Clean up genus_species
            if (genus_species == "Unknown" || genus_species == "") {
                genus_species = "Unknown organism"
            }
            gsub(/[^a-zA-Z0-9 -]/, "", genus_species)
            
            # Output fields: sample, reference, containment, f_match, ANI, genus_species
            # Columns from gather CSV: 
            # intersect_bp, f_orig_query, f_match, f_unique_to_query, f_unique_weighted, 
            # average_abund, median_abund, std_abund, name, filename, md5, f_match_orig, 
            # unique_intersect_bp, gather_result_rank, remaining_bp, query_filename, 
            # query_name, query_md5, query_bp, ksize, moltype, scaled, query_n_hashes, 
            # query_abundance, query_containment_ani, match_containment_ani, 
            # average_containment_ani, max_containment_ani, potential_false_negative
            
            # Get containment (f_orig_query is column 3), f_match (column 4), and ANI
            containment = \$3
            f_match = \$4
            ani = \$25  # query_containment_ani
            
            # Remove quotes from fields
            gsub(/"/, "", containment)
            gsub(/"/, "", f_match)
            gsub(/"/, "", ani)
            
            print sample, \$2, containment, f_match, ani, genus_species
        }
        END {
            if (NR == 1) {
                # No matches found
                print sample, "no_matches_found", "0", "0", "0", "Unknown organism"
            }
        }' sourmash_taxa/${prefix}_gather.csv > sourmash_taxa/${prefix}_taxa.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "WAPHL_ANALYSIS:SOURMASH_TAXA":
        sourmash: \$( sourmash --version 2>&1 | sed 's/sourmash //' )
    END_VERSIONS
    """
}

process SOURMASH_HUMAN_CONTAMINATION {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "sourmash_human/*"   
    container     'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("sourmash_human/*_human_contamination.txt"), emit: human_contamination
    tuple val(meta), path("sourmash_human/*_human_summary.txt"), emit: human_summary
    path "versions.yml", emit: versions
    
    when:
    fasta != null
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ksize = task.ext.ksize ?: '31'
    def scaled = task.ext.scaled ?: '1000'
    """
    mkdir -p sourmash_human
    
    # Sketch the assembly
    sourmash sketch dna \\
        -p k=${ksize},scaled=${scaled},abund \\
        --name ${prefix} \\
        -o sourmash_human/${prefix}.sig \\
        ${fasta}
    
    # Check for human contamination using gather
    # This would ideally use a human genome database
    # For now, we'll search against GenBank and filter for human sequences
    
    if [ -f "/db/genbank-2022.03.genomic.k31.zip" ]; then
        # Run gather against GenBank database
        sourmash gather \\
            sourmash_human/${prefix}.sig \\
            /db/genbank-2022.03.genomic.k31.zip \\
            -o sourmash_human/${prefix}_gather_all.csv \\
            --threshold-bp 50000 \\
            -k ${ksize} || true
        
        # Filter for human sequences and create contamination report
        if [ -f "sourmash_human/${prefix}_gather_all.csv" ]; then
            awk -F',' -v sample=${prefix} '
            BEGIN {
                OFS="\\t"
                print "sample", "human_reference", "containment", "f_match", "ANI", "contamination_level"
                min_containment = 1.0
                contamination_status = "none"
            }
            NR > 1 {
                ref_name = tolower(\$2)
                if (ref_name ~ /homo.*sapiens|human|hg[0-9]+|grch[0-9]+/) {
                    containment = \$3
                    f_match = \$4
                    ani = \$25
                    
                    gsub(/"/, "", containment)
                    gsub(/"/, "", f_match)
                    gsub(/"/, "", ani)
                    
                    containment_val = containment + 0
                    if (containment_val < min_containment) {
                        min_containment = containment_val
                    }
                    
                    # Determine contamination level based on containment
                    # Higher containment = more human contamination
                    if (containment_val > 0.1) {
                        contamination_level = "high"
                    } else if (containment_val > 0.05) {
                        contamination_level = "moderate"
                    } else if (containment_val > 0.01) {
                        contamination_level = "low"
                    } else {
                        contamination_level = "minimal"
                    }
                    
                    print sample, \$2, containment, f_match, ani, contamination_level
                }
            }
            END {
                if (NR == 1) {
                    print sample, "no_human_reference_found", "0", "0", "0", "none"
                }
            }' sourmash_human/${prefix}_gather_all.csv > sourmash_human/${prefix}_human_contamination.txt
            
            # Create summary with lowest containment value (best human match)
            awk -v sample=${prefix} '
            BEGIN {
                OFS="\\t"
                min_containment = 0
                overall_contamination = "none"
            }
            NR > 1 {
                if (\$3 + 0 > min_containment) {
                    min_containment = \$3
                    overall_contamination = \$6
                }
            }
            END {
                print sample, min_containment, overall_contamination
            }' sourmash_human/${prefix}_human_contamination.txt > sourmash_human/${prefix}_human_summary.txt
        else
            # No gather results, create default output
            echo -e "sample\\thuman_reference\\tcontainment\\tf_match\\tANI\\tcontamination_level" > sourmash_human/${prefix}_human_contamination.txt
            echo -e "${prefix}\\tno_human_reference_found\\t0\\t0\\t0\\tnone" >> sourmash_human/${prefix}_human_contamination.txt
            echo -e "${prefix}\\t0\\tnone" > sourmash_human/${prefix}_human_summary.txt
        fi
    else
        echo "Warning: No GenBank database found" >&2
        echo -e "sample\\thuman_reference\\tcontainment\\tf_match\\tANI\\tcontamination_level" > sourmash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\tno_database_found\\t0\\t0\\t0\\tunknown" >> sourmash_human/${prefix}_human_contamination.txt
        echo -e "${prefix}\\t0\\tunknown" > sourmash_human/${prefix}_human_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "WAPHL_ANALYSIS:SOURMASH_HUMAN_CONTAMINATION":
        sourmash: \$( sourmash --version 2>&1 | sed 's/sourmash //' )
    END_VERSIONS
    """
}
