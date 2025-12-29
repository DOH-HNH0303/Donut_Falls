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
    def database_taxa = params.sourmash_db_taxa ?: params.sourmash_db ?: ''
    """
    mkdir -p sourmash_taxa
    
    # Sketch the assembly with specified k-mer size and scaled value
    sourmash sketch dna \\
        -p k=${ksize},scaled=${scaled},abund \\
        --name ${prefix} \\
        -o sourmash_taxa/${prefix}.sig \\
        ${fasta}
    
    # Check if database is provided
    if [ -z "${database_taxa}" ]; then
        echo "WARNING: No sourmash database specified for taxonomic identification!" >&2
        echo "Provide --sourmash_db_taxa (or --sourmash_db) to enable species identification" >&2
        echo "Example databases: GTDB (gtdb-rs214.genomic.k31.zip) or GenBank (genbank-2022.03.genomic.k31.zip)" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tcontainment\\tf_match\\tANI\\tgenus_species" > sourmash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tno_database_provided\\t0\\t0\\t0\\tUnknown organism" >> sourmash_taxa/${prefix}_taxa.txt
        exit 0
    fi
    
    # Check if database file exists
    if [ ! -f "${database_taxa}" ]; then
        echo "WARNING: Sourmash database not found at: ${database_taxa}" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tcontainment\\tf_match\\tANI\\tgenus_species" > sourmash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tdatabase_not_found\\t0\\t0\\t0\\tUnknown organism" >> sourmash_taxa/${prefix}_taxa.txt
        exit 0
    fi
    
    # Database exists, proceed with gather
    sourmash gather \\
        sourmash_taxa/${prefix}.sig \\
        ${database_taxa} \\
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

    cat <<-END_VERSIONS > versions.yml
    "WAPHL_ANALYSIS:SOURMASH_TAXA":
        sourmash: \$( sourmash --version 2>&1 | sed 's/sourmash //' )
    END_VERSIONS
    """
}
