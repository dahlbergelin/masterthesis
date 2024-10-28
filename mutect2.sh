#!/bin/bash

# Set GATK options
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
INPUT_BAM="/gatk/source/raw_data/1d93cae2-d8e3-4e9d-a097-240d3c06e6a1/81e85e60-a0b6-4898-8cb0-02ac1a70f157.rna_seq.genomic.gdc_realn.bam"

# Reference genomes
REF_GENOME="/gatk/project/GRCh38.d1.vd1.fa"
REF_DICT="/gatk/project/GRCh38.d1.vd1.dict"
REF_FAI="/gatk/project/GRCh38.d1.vd1.fa.fai"

# VCF files
KNOWN_SITES="/gatk/project/updated_common_all_20180418.vcf.gz"
PON_VCF="/gatk/project/1000g_pon.hg38.vcf.gz"  # Panel of Normals
GERMLINE_RESOURCE="/gatk/project/af-only-gnomad.hg38.vcf.gz"
CONT_VCF="/gatk/project/small_exac_common_3.hg38.vcf.gz"

# Set unique identifier for each run, e.g., timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Output files
ID_PREFIX=$(basename "$INPUT_BAM" | cut -c1-7) # ID: first seven characters from the BAM file name
OUTPUT_DIR="/gatk/data/elin_output/${ID_PREFIX}" # Unique output directory for each input file
mkdir -p "$OUTPUT_DIR"
FINAL_OUTPUT_BAM="$OUTPUT_DIR/marked_duplicates_$TIMESTAMP.bam"  # Output from MarkDuplicates
SPLIT_OUTPUT_BAM="$OUTPUT_DIR/split_ncigar_reads_$TIMESTAMP.bam"  # Final output file after SplitNCigarReads
METRICS_FILE="$OUTPUT_DIR/marked_dup_metrics_$TIMESTAMP.txt"  # Metrics file
SORTED_BAM="$OUTPUT_DIR/sorted_$TIMESTAMP_.bam"  # Path to the already sorted BAM file, was in /tmp before 
OUTPUT_TABLE="$OUTPUT_DIR/recal_data_$TIMESTAMP.table"
OUTPUT_BAM="$OUTPUT_DIR/recalibrated_output_$TIMESTAMP.bam"  # Final output BAM after applying BQSR
SOMATIC_VCF="$OUTPUT_DIR/somatic_$TIMESTAMP.vcf.gz"  # Output VCF file for somatic variants
OUTPUT_SEG="$OUTPUT_DIR/segments_$TIMESTAMP.table"
PILEUP_TAB="$OUTPUT_DIR/pileups_$TIMESTAMP.table"  # Output for GetPileupSummaries
CONT_TAB="$OUTPUT_DIR/contamination_$TIMESTAMP.table"
FILT_MUTCAL="$OUTPUT_DIR/filtered_without_seg_$TIMESTAMP.vcf.gz"
INTERVAL_LIST="$OUTPUT_DIR/hg38_exons_$TIMESTAMP.interval_list"
ALLELIC_COUNTS="$OUTPUT_DIR/allelic_counts_hets_$TIMESTAMP.tsv"

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
samtools sort "$INPUT_BAM" -o "$SORTED_BAM"

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
    -I "$SORTED_BAM" \
    -O "$FINAL_OUTPUT_BAM" \
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
    -R "$REF_GENOME" \
    -I "$FINAL_OUTPUT_BAM" \
    -O "$SPLIT_OUTPUT_BAM"

# Check if SplitNCigarReads completed successfully
if [ $? -eq 0 ]; then
    log "SplitNCigarReads completed successfully."
else
    log "SplitNCigarReads failed."
    exit 1
fi


# Run Base Quality Recalibration on the output from SplitNCigarReads
log "Running BaseRecalibrator..."
gatk --java-options '-Xmx4g' BaseRecalibrator \
    -I "$SPLIT_OUTPUT_BAM" \
    -R "$REF_GENOME" \
    --known-sites "$KNOWN_SITES" \
    -O "$OUTPUT_TABLE"

if [ $? -ne 0 ]; then
    log "BaseRecalibrator failed."
    exit 1
else
    log "BaseRecalibrator completed successfully."
fi

# Apply Base Quality Score Recalibration
log "Applying Base Quality Score Recalibration..."
gatk --java-options '-Xmx4g' ApplyBQSR \
    -R "$REF_GENOME" \
    -I "$SPLIT_OUTPUT_BAM" \
    --bqsr-recal-file "$OUTPUT_TABLE" \
    -O "$OUTPUT_BAM"  # Ensure OUTPUT_BAM includes the file name

if [ $? -ne 0 ]; then
    log "ApplyBQSR failed."
    exit 1
else
    log "ApplyBQSR completed successfully."
fi

# Variant Calling with Mutect2
log "Running Mutect2 for variant calling..."
gatk Mutect2 \
    -R "$REF_GENOME" \
    -I "$OUTPUT_BAM" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --panel-of-normals "$PON_VCF" \ 
    -O "$SOMATIC_VCF"

if [ $? -ne 0 ]; then
    log "Mutect2 variant calling failed."
    exit 1
else
    log "Mutect2 variant calling completed successfully."
fi


# Index the VCF file using tabix
log "Indexing the VCF file with tabix..."
tabix -p vcf "$CONT_VCF"

# Calculate contamination
log "Running GetPileupSummaries..."
gatk GetPileupSummaries \
    -I "$SPLIT_OUTPUT_BAM" \
    -V "$CONT_VCF" \
    -L "$CONT_VCF" \
    -O "$PILEUP_TAB"  

if [ $? -ne 0 ]; then
    log "GetPileupSummaries failed."
    exit 1
else
    log "GetPileupSummaries completed successfully."
fi

# Calculate contamination
log "Running CalculateContamination..."
gatk CalculateContamination \
    -I "$PILEUP_TAB" \
    -O "$CONT_TAB" \
    -tumor-segmentation "$OUTPUT_SEG"

if [ $? -ne 0 ]; then
    log "CalculateContamination failed."
    exit 1
else
    log "CalculateContamination completed successfully."
fi

# Filter Variants
log "Running FilterMutectCalls..."
gatk FilterMutectCalls \
    -R "$REF_GENOME" \
    -V "$SOMATIC_VCF" \
    --contamination-table "$CONT_TAB" \
    --tumor-segmentation "$OUTPUT_SEG" \
    -O "$FILT_MUTCAL"

if [ $? -ne 0 ]; then
    log "FilterMutectCalls failed."
    exit 1
else
    log "FilterMutectCalls completed successfully."
fi

# Completion message
log "Pipeline finished!"




:'

# Learn Orientation Bias Artifacts (if needed)
# Add command here if applicable

# Filter Variants
echo "Running FilterMutectCalls..."
gatk FilterMutectCalls \
    -R $REFERENCE \
    -V somatic.vcf.gz \
    --contamination-table contamination.table \
    --tumor-segmentation segments.tsv \
    -O filtered.vcf.gz

# Annotate Variants (if needed)
# Add command here if applicable
'
#this wil run if the code succed
echo "Pipeline finished!"


