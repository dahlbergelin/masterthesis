#!/bin/bash

# Paths and variables
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
INPUT_BAM="/gatk/source/raw_data/1d93cae2-d8e3-4e9d-a097-240d3c06e6a1/81e85e60-a0b6-4898-8cb0-02ac1a70f157.rna_seq.genomic.gdc_realn.bam"
FUNCOTATOR_DATA_SOURCES="/gatk/project/funcotator_dataSources.v1.8.hg38.20230908s"

# Reference genomes
REF_GENOME="/gatk/project/GRCh38.d1.vd1.fa"
REF_DICT="/gatk/project/GRCh38.d1.vd1.dict"
REF_FAI="/gatk/project/GRCh38.d1.vd1.fa.fai"

# VCF files
KNOWN_SITES="/gatk/project/updated_common_all_20180418.vcf.gz" # BaseRecalibrator
PON_VCF="/gatk/project/1000g_pon.hg38.vcf.gz"  # Mutect2: Panel of Normals 1000 genomes
GERMLINE_RESOURCE="/gatk/project/af-only-gnomad.hg38.vcf.gz" # Mutect2: gmomAD  
CONT_VCF="/gatk/project/small_exac_common_3.hg38.vcf.gz" # GetPileUpSummeries

ID_PREFIX=$(basename "$INPUT_BAM" | cut -c1-7)
OUTPUT_DIR="/gatk/data/elin_output/${ID_PREFIX}"
mkdir -p "$OUTPUT_DIR"

# Output file paths
SORTED_BAM="$OUTPUT_DIR/sorted_${ID_PREFIX}.bam"
FINAL_OUTPUT_BAM="$OUTPUT_DIR/marked_duplicates_${ID_PREFIX}.bam"
METRICS_FILE="$OUTPUT_DIR/marked_dup_metrics_${ID_PREFIX}.txt"
SPLIT_OUTPUT_BAM="$OUTPUT_DIR/split_ncigar_reads_${ID_PREFIX}.bam"
OUTPUT_TABLE="$OUTPUT_DIR/recal_data_${ID_PREFIX}.table"
OUTPUT_BAM="$OUTPUT_DIR/recalibrated_output_${ID_PREFIX}.bam"
SOMATIC_VCF="$OUTPUT_DIR/somatic_${ID_PREFIX}.vcf.gz"
PILEUP_TAB="$OUTPUT_DIR/pileups_${ID_PREFIX}.table"
CONT_TAB="$OUTPUT_DIR/contamination_${ID_PREFIX}.table"
FILT_MUTCAL="$OUTPUT_DIR/filtered_without_seg_${ID_PREFIX}.vcf.gz"
F1R2_FILE="$OUTPUT_DIR/f1r2_${ID_PREFIX}.tar.gz"         # Intermediate F1R2 file for orientation model
ORIENTATION_MODEL="$OUTPUT_DIR/read-orientation-model_${ID_PREFIX}.tar.gz"  # Output from LearnReadOrientationModel
FUNCOTATED_VCF="$OUTPUT_DIR/funcotated_output_${ID_PREFIX}.vcf"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Check for existing index
if [ ! -f "${REF_GENOME}.fai" ]; then
    log "Indexing the reference genome..."
    samtools faidx "$REF_GENOME"
    log "Reference genome indexing completed."
else
    log "Reference genome index already exists."
fi

# Sort BAM file
log "Starting samtools sort..."
samtools sort "$INPUT_BAM" -o "$SORTED_BAM"
if [ $? -ne 0 ]; then log "samtools sort failed."; exit 1; fi

# Mark Duplicates
log "Running MarkDuplicates..."
java -jar "$GATK_JAR" MarkDuplicates \
    -I "$SORTED_BAM" \
    -O "$FINAL_OUTPUT_BAM" \
    -M "$METRICS_FILE" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true
if [ $? -ne 0 ]; then log "MarkDuplicates failed."; exit 1; fi

# Split NCigar Reads
log "Running SplitNCigarReads..."
gatk SplitNCigarReads \
    -R "$REF_GENOME" \
    -I "$FINAL_OUTPUT_BAM" \
    -O "$SPLIT_OUTPUT_BAM"
if [ $? -ne 0 ]; then log "SplitNCigarReads failed."; exit 1; fi

# Base Quality Recalibration
log "Running BaseRecalibrator..."
gatk BaseRecalibrator \
    -I "$SPLIT_OUTPUT_BAM" \
    -R "$REF_GENOME" \
    --known-sites "$KNOWN_SITES" \
    -O "$OUTPUT_TABLE"
if [ $? -ne 0 ]; then log "BaseRecalibrator failed."; exit 1; fi

# Apply BQSR
log "Applying Base Quality Score Recalibration..."
gatk ApplyBQSR \
    -R "$REF_GENOME" \
    -I "$SPLIT_OUTPUT_BAM" \
    --bqsr-recal-file "$OUTPUT_TABLE" \
    -O "$OUTPUT_BAM"
if [ $? -ne 0 ]; then log "ApplyBQSR failed."; exit 1; fi

# Mutect2 Variant Calling
log "Running Mutect2 for variant calling..."
gatk Mutect2 \
    -R "$REF_GENOME" \
    -I "$OUTPUT_BAM" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --panel-of-normals "$PON_VCF" \
    -O "$SOMATIC_VCF" \
    --f1r2-tar-gz "$F1R2_FILE"
if [ $? -ne 0 ]; then log "Mutect2 variant calling failed."; exit 1; fi

# Run LearnReadOrientationModel to process the F1R2 file
log "Running LearnReadOrientationModel..."
gatk LearnReadOrientationModel \
    -I "$F1R2_FILE" \
    -O "$ORIENTATION_MODEL"
if [ $? -ne 0 ]; then log "LearnReadOrientationModel failed."; exit 1; fi

# Indexing VCF
log "Indexing the VCF file with tabix..."
tabix -p vcf "$SOMATIC_VCF"

# Get Pileup Summaries
log "Running GetPileupSummaries..."
gatk GetPileupSummaries \
    -I "$OUTPUT_BAM" \
    -V "$CONT_VCF" \
    -L "$CONT_VCF" \
    -O "$PILEUP_TAB"
if [ $? -ne 0 ]; then log "GetPileupSummaries failed."; exit 1; fi

# Calculate Contamination
log "Running CalculateContamination..."
gatk CalculateContamination \
    -I "$PILEUP_TAB" \
    -O "$CONT_TAB"
if [ $? -ne 0 ]; then log "CalculateContamination failed."; exit 1; fi


# Filter Mutect Calls
log "Running FilterMutectCalls..."
gatk FilterMutectCalls \
    -R "$REF_GENOME" \
    -V "$SOMATIC_VCF" \
    --contamination-table "$CONT_TAB" \
    --orientation-bias-artifact-priors "$ORIENTATION_MODEL" \
    -O "$FILT_MUTCAL"
if [ $? -ne 0 ]; then log "FilterMutectCalls failed."; exit 1; fi


gatk --java-options "-Xmx8g" Funcotator \
    -R "$REF_GENOME" \
    -V "filtered.vcf" \
    --ref-version hg38 \
    --data-sources-path "$FUNCOTATOR_DATA_SOURCES" \
    -O "$FUNCOTATED_VCF" \
    --output-file-format VCF \
    --transcript-selection-mode CANONICAL
if [ $? -ne 0 ]; then log "Funcotator VCF generation failed."; exit 1; fi

:'
gatk VariantsToTable \
    -v "$FUNCOTATED_VCF" -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O "/gatk/data/elin_output/${ID_PREFIX}/output-funct.table"
'
# Completion message
log "Pipeline finished!"

