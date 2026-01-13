process HOSTILE {
    tag "${meta.id}"
    label 'process_high'
    
    conda (params.enable_conda ? 'bioconda::hostile=2.0.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0' :
        'biocontainers/hostile:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val index_name

    output:
    tuple val(meta), path("*.clean.fastq.gz"), emit: reads
    path "*_hostile_stats.txt"                , emit: stats
    path "*_hostile_log.json"                 , emit: log
    path "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_arg = index_name ? "--index ${index_name}" : ""
    def threads_arg = task.cpus ? "--threads ${task.cpus}" : ""
    
    """
    # Run hostile and capture JSON log
    hostile clean \\
        --fastq1 ${reads} \\
        ${index_arg} \\
        ${threads_arg} \\
        --output . \\
        > ${prefix}_hostile_log.json
    
    # Parse JSON log to create human-readable stats file
    python3 <<CODE
import json
import sys

try:
    with open('${prefix}_hostile_log.json', 'r') as f:
        # Skip any INFO lines and find the JSON array
        lines = f.readlines()
        json_start = -1
        for i, line in enumerate(lines):
            if line.strip().startswith('['):
                json_start = i
                break
        
        if json_start >= 0:
            json_text = ''.join(lines[json_start:])
            data = json.loads(json_text)
            
            # Hostile returns a list with one entry for single FASTQ input
            if data and len(data) > 0:
                result = data[0]
                
                reads_in = result.get('reads_in', 0)
                reads_out = result.get('reads_out', 0)
                reads_removed = result.get('reads_removed', 0)
                reads_removed_prop = result.get('reads_removed_proportion', 0.0)
                
                # Write stats file
                with open('${prefix}_hostile_stats.txt', 'w') as out:
                    out.write(f"sample\\t${meta.id}\\n")
                    out.write(f"reads_in\\t{reads_in}\\n")
                    out.write(f"reads_out\\t{reads_out}\\n")
                    out.write(f"reads_removed\\t{reads_removed}\\n")
                    out.write(f"reads_removed_proportion\\t{reads_removed_prop:.6f}\\n")
                    out.write(f"human_contamination_pct\\t{reads_removed_prop * 100:.4f}\\n")
            else:
                sys.stderr.write("Warning: No data found in hostile JSON output\\n")
                with open('${prefix}_hostile_stats.txt', 'w') as out:
                    out.write(f"sample\\t${meta.id}\\n")
                    out.write(f"reads_in\\t0\\n")
                    out.write(f"reads_out\\t0\\n")
                    out.write(f"reads_removed\\t0\\n")
                    out.write(f"reads_removed_proportion\\t0\\n")
                    out.write(f"human_contamination_pct\\t0\\n")
        else:
            sys.stderr.write("Warning: Could not find JSON array in hostile log\\n")
            with open('${prefix}_hostile_stats.txt', 'w') as out:
                out.write(f"sample\\t${meta.id}\\n")
                out.write(f"reads_in\\t0\\n")
                out.write(f"reads_out\\t0\\n")
                out.write(f"reads_removed\\t0\\n")
                out.write(f"reads_removed_proportion\\t0\\n")
                out.write(f"human_contamination_pct\\t0\\n")
                
except Exception as e:
    sys.stderr.write(f"Error parsing hostile log: {str(e)}\\n")
    with open('${prefix}_hostile_stats.txt', 'w') as out:
        out.write(f"sample\\t${meta.id}\\n")
        out.write(f"reads_in\\t0\\n")
        out.write(f"reads_out\\t0\\n")
        out.write(f"reads_removed\\t0\\n")
        out.write(f"reads_removed_proportion\\t0\\n")
        out.write(f"human_contamination_pct\\t0\\n")
CODE

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version 2>&1 | sed 's/hostile //; s/ .*//')
    END_VERSIONS
    """
}
