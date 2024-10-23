#!/bin/bash

# Set GATK options
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
INPUT_BAM="/source/raw_data/1d93cae2-d8e3-4e9d-a097-240d3c06e6a1/81e85e60-a0b6-4898-8cb0-02ac1a70f157.rna_seq.genomic.gdc_realn.bam"  
OUTPUT_DIR="/gatk/data/elin_output"
FINAL_OUTPUT_BAM="$OUTPUT_DIR/marked_duplicates.bam"  # Output from MarkDuplicates
SPLIT_OUTPUT_BAM="$OUTPUT_DIR/split_ncigar_reads.bam"  # Final output file after SplitNCigarReads
METRICS_FILE="$OUTPUT_DIR/marked_dup_metrics.txt"  # Metrics file (optional)

# Set reference genome path to the new file
REF_GENOME="/gatk/project/GRCh38.d1.vd1.fa"  # Update with the new reference genome file

# Create a temporary directory for samtools
TMP_DIR="/gatk/project/tmp"
mkdir -p "$TMP_DIR"

# Function to log the current time and message
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Check for the index (.fai) in the same directory as the reference genome
if [ ! -f "${REF_GENOME}.fai" ]; then
    log "Indexing the reference genome..."
    samtools faidx "$REF_GENOME"  # Index file (.fai) will be created here
    log "Reference genome indexing completed."
else
    log "Reference genome index already exists."
fi

# Run samtools sort (skip if already sorted)
log "Starting samtools sort..."
samtools sort -n -o "$TMP_DIR/sorted.bam" "$INPUT_BAM" --threads 4

# Check if samtools sort completed successfully
if [ $? -ne 0 ]; then
    log "samtools sort failed."
    exit 1
else
    log "samtools sort completed successfully."
fi

# Run MarkDuplicates
log "Running MarkDuplicates..."
java -Dsamjdk.use_async_io_read_samtools=false \
    -Dsamjdk.use_async_io_write_samtools=true \
    -Dsamjdk.use_async_io_write_tribble=false \
    -Dsamjdk.compression_level=2 \
    -jar "$GATK_JAR" MarkDuplicates \
    -I "$TMP_DIR/sorted.bam" \  # Input from the temporary sorted file
    -O "$FINAL_OUTPUT_BAM" \  # Output to final output file
    -M "$METRICS_FILE" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true \
    --VERBOSITY DEBUG

# Check if MarkDuplicates completed successfully
if [ $? -ne 0 ]; then
    log "MarkDuplicates failed."
    exit 1
else
    log "MarkDuplicates completed successfully."
fi

# Use SplitNCigarReads with stack trace enabled
log "Running SplitNCigarReads..."
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SplitNCigarReads \
    -R "$REF_GENOME" \  # Use the correct reference genome path
    -I "$FINAL_OUTPUT_BAM" \
    -O "$SPLIT_OUTPUT_BAM"  # Output the result to a new file

# Check if SplitNCigarReads completed successfully
if [ $? -eq 0 ]; then
    log "SplitNCigarReads completed successfully."
else
    log "SplitNCigarReads failed."
    exit 1
fi

# Cleanup temporary files
log "Cleaning up temporary files..."
rm -f "$TMP_DIR/sorted.bam"
rm -rf "$TMP_DIR"
log "Temporary files cleaned up."

# Final success message
log "Processing completed successfully."
