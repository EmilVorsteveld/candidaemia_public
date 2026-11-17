#!/bin/bash

# This is code to run isocall on the bulk kinnex samples

# Isoforms have been prepared as follows:
#isocall-v0.14.0 prep-isoforms --gtf ./gencode.v49.annotation.gtf.gz --output ./gencode.v49.isoforms.gz

# ==========================================
# INPUT & OUTPUT HANDLING
# ==========================================

# $1 captures the first argument (BAM input directory)
# If empty, defaults to the current directory "."
BAM_DIR="${1:-.}"
BAM_DIR="${BAM_DIR%/}" # Remove trailing slash

# $2 captures the second argument (Output directory)
# If empty, defaults to the current directory "."
OUTPUT_DIR="${2:-.}"
OUTPUT_DIR="${OUTPUT_DIR%/}" # Remove trailing slash

# Safety Check 1: Verify the input directory actually exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: Input directory '$BAM_DIR' does not exist."
    echo "Usage: $0 <input_bam_dir> [output_dir]"
    exit 1
fi

# Safety Check 2: If the output directory doesn't exist, create it automatically
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory '$OUTPUT_DIR' does not exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

echo "Input directory  : $BAM_DIR"
echo "Output directory : $OUTPUT_DIR"
echo "--------------------------------------------------"

# ==========================================
# RUN PROCESSING
# ==========================================

# Count how many BAM files are present
bam_count=$(ls -1 "$BAM_DIR"/*.bam 2>/dev/null | wc -l)

if [ "$bam_count" -eq 0 ]; then
    echo "Error: No .bam files found in '$BAM_DIR'"
    exit 1
fi

echo "Found $bam_count BAM file(s). Starting processing..."


# Loop through all .bam files in the directory
for bam_path in "$INPUT_DIR"/*.bam; do
    
    # Check if any BAM files actually exist to avoid loop errors
    [ -e "$bam_path" ] || continue
    
    # Extract just the filename (e.g., RNA25-00042.isoseq-refine.mapped...)
    file_name=$(basename "$bam_path")
    
    # Extract the sample ID (everything before the first dot '.')
    sample_id="${file_name%%.*}"
    
    # Define the output name
    output_name="${sample_id}.gz"
    
    echo "--------------------------------------------------"
    echo "Processing sample: $sample_id"
    echo "--------------------------------------------------"
    
    # Run the isocall command
    isocall-v0.14.0 profile \
        --reads "$bam_path" \
        --sample "$sample_id" \
        --output "$output_name"

done

echo "All files processed successfully!"


# Merge profiles
echo "Merging the profiles..."
isocall-v0.14.0 merge --profiles *.gz --output merged.gz
echo "Profiles have been merged"

# Call isoforms
echo "Calling isoforms..."
isocall-v0.14.0 call --merged-profile merged.gz --known-isoforms gencode.v49.isoforms.gz --reference /mnt/data/inbox/hg38/hg38/hg38.fa --output-prefix merged
echo "Isoforms have been called"

# Prepare gtf files for Oarfish
echo "Prepare GTF for Oarfish..."
## Generate a fasta file for the matched isoforms
gunzip merged.isoforms.gtf.gz
gffread -w merged.isoforms.fa -g /mnt/data/inbox/hg38/hg83/hg38.fa merged.isoforms.gtf

## Run Oarfish
echo "Counting with Oarfish..."
for bam_path in "$BAM_DIR"/*.bam; do
    
    [ -e "$bam_path" ] || continue
    
    # Extract filename and sample ID
    file_name=$(basename "$bam_path")
    sample_id="${file_name%%.*}"
    
    # Define sample-specific output folder inside your configured output directory
    sample_output_dir="$OUTPUT_DIR/${sample_id}"
    
    echo "--------------------------------------------------"
    echo "Running oarfish for sample: $sample_id"
    echo "Input:  $bam_path"
    echo "Output: $sample_output_dir"
    echo "--------------------------------------------------"
    
    oarfish -j 15 \
        --reads "$bam_path" \
        --annotated merged.isoforms.fa \
        --filter-group no-filters \
        --seq-tech pac-bio-hifi \
        --output "$sample_output_dir"
done

echo "Done running Oarfish for all samples!"