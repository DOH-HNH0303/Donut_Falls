process SOURMASH_TAXA {
    tag           "${meta.id}"
    label         "process_medium"
    publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "sourmash_taxa/*"   
    container     'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'
    
    input:
    tuple val(meta), path(fasta)
    path database_taxa, stageAs: 'sourmash_db/*'

    output:
    tuple val(meta), path("sourmash_taxa/*.txt"), emit: taxa
    path "versions.yml", emit: versions
    
    when:
    fasta != null
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ksize = task.ext.ksize ?: '31'
    def scaled = task.ext.scaled ?: '1000'
    // database_taxa is now a staged file path from the input
    """
    mkdir -p sourmash_taxa
    
    # Sketch the assembly with specified k-mer size and scaled value
    sourmash sketch dna \\
        -p k=${ksize},scaled=${scaled},abund \\
        --name ${prefix} \\
        -o sourmash_taxa/${prefix}.sig \\
        ${fasta}
    
    # Check if database file is empty (happens when [] is passed as input)
    if [ "${database_taxa}" == "[]" ] || [ ! -s "${database_taxa}" ]; then
        echo "WARNING: No sourmash database specified for taxonomic identification!" >&2
        echo "Provide --sourmash_db_taxa (or --sourmash_db) to enable species identification" >&2
        echo "Example databases: GTDB (gtdb-rs214.genomic.k31.zip) or GenBank (genbank-2022.03.genomic.k31.zip)" >&2
        echo "Creating empty output..." >&2
        echo -e "sample\\treference\\tcontainment\\tf_match\\tANI\\tgenus_species" > sourmash_taxa/${prefix}_taxa.txt
        echo -e "${prefix}\\tno_database_provided\\t0\\t0\\t0\\tUnknown organism" >> sourmash_taxa/${prefix}_taxa.txt
        exit 0
    fi
    
    # Database exists and is staged, proceed with gather
    sourmash gather \\
        sourmash_taxa/${prefix}.sig \\
        ${database_taxa} \\
        -o sourmash_taxa/${prefix}_gather.csv \\
        --threshold-bp 50000 \\
        -k ${ksize}
    
    # Parse gather results to create taxa output compatible with downstream processing
    # Filter for match_containment_ani >= 0.95 and extract taxa (without accession)
    awk -F',' -v sample=${prefix} '
    BEGIN {
        OFS="\\t"
    }
    NR == 1 {
        # Print header
        print "sample", "gather_rank", "taxa", "match_containment_ani", "f_unique_to_query", "average_abund"
        next
    }
    NR > 1 {
        # Columns from gather CSV (1-indexed in awk):
        # 1=intersect_bp, 2=f_orig_query, 3=f_match, 4=f_unique_to_query, 5=f_unique_weighted, 
        # 6=average_abund, 7=median_abund, 8=std_abund, 9=filename, 10=name, 11=md5, 12=f_match_orig, 
        # 13=unique_intersect_bp, 14=gather_result_rank, 15=remaining_bp, 16=query_filename, 
        # 17=query_name, 18=query_md5, 19=query_bp, 20=ksize, 21=moltype, 22=scaled, 23=query_n_hashes, 
        # 24=query_abundance, 25=query_containment_ani, 26=match_containment_ani, 
        # 27=average_containment_ani, 28=max_containment_ani, 29=potential_false_negative
        
        # Extract fields
        full_name = \$10
        gather_rank = \$14
        f_unique = \$4
        avg_abund = \$6
        match_ani = \$26
        
        # Remove quotes
        gsub(/"/, "", full_name)
        gsub(/"/, "", gather_rank)
        gsub(/"/, "", f_unique)
        gsub(/"/, "", avg_abund)
        gsub(/"/, "", match_ani)
        
        # Filter for match_containment_ani >= 0.95
        if (match_ani + 0 >= 0.95) {
            # Extract taxa name (everything after first space, removing accession)
            taxa = full_name
            # Split on first space and take rest
            if (match(full_name, / /)) {
                taxa = substr(full_name, RSTART + 1)
            }
            
            # Clean up taxa name
            if (taxa == "" || taxa == full_name) {
                taxa = "Unknown organism"
            }
            
            # Output directly - gather CSV is already sorted by rank
            print sample, gather_rank, taxa, match_ani, f_unique, avg_abund
        }
    }' sourmash_taxa/${prefix}_gather.csv > sourmash_taxa/${prefix}_taxa_tmp.txt
    
    # Check if we got any matches (file has more than just header)
    if [ \$(wc -l < sourmash_taxa/${prefix}_taxa_tmp.txt) -le 1 ]; then
        # No matches above threshold - add a placeholder line
        echo -e "${prefix}\\t0\\tNo matches >= 95% ANI\\t0\\t0\\t0" >> sourmash_taxa/${prefix}_taxa_tmp.txt
    fi
    
    # Move temp file to final location
    mv sourmash_taxa/${prefix}_taxa_tmp.txt sourmash_taxa/${prefix}_taxa.txt

    cat <<-END_VERSIONS > versions.yml
    "WAPHL_ANALYSIS:SOURMASH_TAXA":
        sourmash: \$( sourmash --version 2>&1 | sed 's/sourmash //' )
    END_VERSIONS
    """
}
