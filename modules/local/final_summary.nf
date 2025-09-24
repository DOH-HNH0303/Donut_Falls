#!/usr/bin/env nextflow


process FINAL_SUMMARY {
    tag "Creating final summary"
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path(contamination_files)
    path(taxa_files) 
    path(coverage_files)
    path(assembly_files)
    
    output:
    path("waphl_final_summary.tsv"), emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    
    import glob
    import pandas as pd
    import os
    from pathlib import Path
    
    # Create output directory
    os.makedirs("human_contamination", exist_ok=True)
    
    # Process contamination files with unique naming
    contamination_files = glob.glob("*human_contamination.txt") + glob.glob("*human_summary.txt")
    
    # Group files by sample ID and type to avoid naming conflicts
    processed_files = {}
    
    for file_path in contamination_files:
        filename = os.path.basename(file_path)
        
        # Extract sample ID from filename
        if "_human_contamination.txt" in filename:
            sample_id = filename.replace("_human_contamination.txt", "")
            file_type = "contamination"
        elif "_human_summary.txt" in filename:
            sample_id = filename.replace("_human_summary.txt", "")
            file_type = "summary"
        else:
            continue
            
        # Create unique identifier
        unique_key = f"{sample_id}_{file_type}"
        
        # Only process if we haven't seen this combination
        if unique_key not in processed_files:
            processed_files[unique_key] = file_path
            
            # Copy with unique name to avoid conflicts
            new_name = f"human_contamination/{sample_id}_{file_type}.txt"
            os.system(f"cp '{file_path}' '{new_name}'")
    
    # Process other input files similarly
    taxa_files = glob.glob("*taxa*")
    coverage_files = glob.glob("*coverage*") 
    assembly_files = glob.glob("*assembly*") + glob.glob("*stats*")
    
    # Create the final summary
    summary_data = []
    
    # Get unique sample IDs
    sample_ids = set()
    for key in processed_files.keys():
        sample_id = key.split("_")[0]  # Get the part before first underscore
        sample_ids.add(sample_id)
    
    # Process each sample
    for sample_id in sorted(sample_ids):
        sample_data = {"sample": sample_id}
        
        # Add contamination data
        contamination_file = f"human_contamination/{sample_id}_contamination.txt"
        if os.path.exists(contamination_file):
            with open(contamination_file, 'r') as f:
                content = f.read().strip()
                sample_data["human_contamination"] = content
        
        # Add summary data  
        summary_file = f"human_contamination/{sample_id}_summary.txt"
        if os.path.exists(summary_file):
            with open(summary_file, 'r') as f:
                content = f.read().strip()
                sample_data["human_summary"] = content
        
        # Add taxa information
        taxa_file = [f for f in taxa_files if sample_id in f]
        if taxa_file:
            with open(taxa_file[0], 'r') as f:
                content = f.read().strip()
                sample_data["taxa"] = content
        
        # Add coverage information
        coverage_file = [f for f in coverage_files if sample_id in f]
        if coverage_file:
            with open(coverage_file[0], 'r') as f:
                content = f.read().strip()
                sample_data["coverage"] = content
        
        summary_data.append(sample_data)
    
    # Convert to DataFrame and save
    df = pd.DataFrame(summary_data)
    df.to_csv("waphl_final_summary.tsv", sep="\\t", index=False)
    
    print(f"Processed {len(summary_data)} samples")
    print("Files processed:")
    for key, file_path in processed_files.items():
        print(f"  {key}: {file_path}")
    """
}

/*
 * Alternative fix: Modify the input channel handling to avoid name collisions
 * This approach fixes the issue at the channel level before it reaches FINAL_SUMMARY
 */

workflow WAPHL_ANALYSIS_FIXED {
    take:
    contamination_ch
    taxa_ch
    coverage_ch
    assembly_ch
    
    main:
    // Fix naming conflicts by adding unique prefixes
    contamination_fixed = contamination_ch
        .map { sample, files ->
            def fixed_files = []
            files.each { file ->
                def basename = file.getBaseName()
                def extension = file.getExtension()
                def newName = "${sample}_${basename}.${extension}"
                fixed_files.add(file.copyTo(newName))
            }
            [sample, fixed_files]
        }
        .flatten()
        .filter { it.toString().endsWith('.txt') }
        .collect()
    
    // Similarly fix other channels
    taxa_fixed = taxa_ch.collect()
    coverage_fixed = coverage_ch.collect() 
    assembly_fixed = assembly_ch.collect()
    
    // Run the fixed final summary process
    FINAL_SUMMARY(
        contamination_fixed,
        taxa_fixed,
        coverage_fixed, 
        assembly_fixed
    )
    
    emit:
    summary = FINAL_SUMMARY.out.summary
}

/*
 * Quick fix option: Modify the existing process to handle duplicates
 * Add this to your main.nf file or as a separate module
 */

process FINAL_SUMMARY_DEDUP {
    tag "Creating final summary with deduplication"
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path(input_files, stageAs: "input/*")
    
    output:
    path("waphl_final_summary.tsv"), emit: summary
    
    script:
    """
    #!/bin/bash
    
    # Create unique file mapping to avoid name conflicts
    mkdir -p processed_files
    
    # Process each input file and rename to avoid conflicts
    for file in input/*; do
        if [ -f "\$file" ]; then
            filename=\$(basename "\$file")
            # Add timestamp or hash to make unique if duplicate names exist
            if [ -f "processed_files/\$filename" ]; then
                timestamp=\$(date +%s%N)
                new_name="processed_files/\${timestamp}_\$filename"
                cp "\$file" "\$new_name"
            else
                cp "\$file" "processed_files/\$filename"
            fi
        fi
    done
    
    # Now run the original summary script with deduplicated files
    cd processed_files
    
    # Your original summary generation code here
    python3 ${projectDir}/bin/create_summary.py
    
    # Move result to expected location
    mv summary.tsv ../waphl_final_summary.tsv
    """
}