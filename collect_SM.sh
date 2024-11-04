#!/bin/bash

# Define the base directory containing the BAM files
BASE_DIR="/gatk/source/raw_data"

# Output file for the table
OUTPUT_DIR="/gatk/project/sample_id_table"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"


OUTPUT_FILE="$OUTPUT_DIR/sample_id_table.txt"

# Print the header for the table in the output file
echo -e "Sample_ID\tBAM_File_Path" > "$OUTPUT_FILE"

# Loop through all subdirectories and find BAM files
find "$BASE_DIR" -type f -name "*.bam" | while read -r BAM_FILE; do
    # Extract the sample ID (SM tag) using samtools
    SAMPLE_ID=$(samtools view -H "$BAM_FILE" | grep '^@RG' | awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}' | sed 's/SM://')

    # Check if SAMPLE_ID was found
    if [ -n "$SAMPLE_ID" ]; then
        # Append the sample ID and file path to the output file
        echo -e "$SAMPLE_ID\t$BAM_FILE" >> "$OUTPUT_FILE"
    else
        echo "No sample ID found for BAM file: $BAM_FILE"
    fi
done

# Print a message indicating the table has been created
echo "Sample ID table created at $OUTPUT_FILE"
