#!/bin/bash
export PATH="/home/qluelda/miniconda3/envs/vardict_env/bin:$PATH"

# Paths and configurations
REF_GENOME="GRCh38.d1.vd1.fa"
AF_THRESHOLD="0.01"
OUTPUT_ROOT="/data/elin_output/VARDICT"

echo "PATH is: $PATH"

# Function to log messages
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Check if VarDict is installed
if ! command -v vardict-java &>/dev/null; then
    echo "VarDict is not installed or accessible in your PATH."
    exit 1
fi

# Check if reference genome is indexed
if [ ! -f "${REF_GENOME}.fai" ]; then
    log "Indexing the reference genome..."
    samtools faidx "$REF_GENOME"
    log "Reference genome indexing completed."
else
    log "Reference genome index already exists."
fi

# Loop over each sample directory in /data/output/EGFR_output
for sample_dir in /data/output/EGFR_output/*/; do
    # Find the recalibrated BAM file in each directory
    BAM_INPUT="${sample_dir}recalibrated_output_reads.bam"
    
    # Check if BAM file exists, skip if not
    if [ ! -f "$BAM_INPUT" ]; then
        log "BAM file not found in $sample_dir, skipping..."
        continue
    fi

    # Extract sample ID and prepare output directory
    ID_PREFIX=$(basename "$sample_dir")
    OUTPUT_DIR="${OUTPUT_ROOT}/${ID_PREFIX}"
    mkdir -p "$OUTPUT_DIR"

    # Define output VCF file
    OUTPUT_VCF="${OUTPUT_DIR}/${ID_PREFIX}_variants.vardict.vcf"

    log "Processing $BAM_INPUT with VarDict, output to $OUTPUT_VCF"

    # Run VarDict on the recalibrated BAM file
    vardict-java -G "$REF_GENOME" -f "$AF_THRESHOLD" -N "$ID_PREFIX" -b "$BAM_INPUT" -z -F 0 -c 1 -S 2 -E 3 -g 4 -R chr7:55019017-55211628:EGFR |
        teststrandbias.R |
        var2vcf_valid.pl -N "$ID_PREFIX" -E -f "$AF_THRESHOLD" >"$OUTPUT_VCF" || log "Error in processing $ID_PREFIX"

    log "Completed processing for $ID_PREFIX, results saved to $OUTPUT_VCF"
done

log "All BAM files processed!"