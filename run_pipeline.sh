#!/bin/bash

# Set GATK options
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
INPUT_BAM="/gatk/source/raw_data/1d93cae2-d8e3-4e9d-a097-240d3c06e6a1/81e85e60-a0b6-4898-8cb0-02ac1a70f157.rna_seq.genomic.gdc_realn.bam"  
OUTPUT_DIR="/gatk/data/elin_output"
SORTED_BAM="/gatk/project/tmp/sorted.bam"  # Path to the already sorted BAM file
FINAL_OUTPUT_BAM="$OUTPUT_DIR/marked_duplicates.bam"  # Output from MarkDuplicates
SPLIT_OUTPUT_BAM="$OUTPUT_DIR/split_ncigar_reads.bam"  # Final output file after SplitNCigarReads
METRICS_FILE="$OUTPUT_DIR/marked_dup_metrics.txt"  # Metrics file (optional)

# Set reference genome path to the new file
REF_GENOME="/gatk/project/GRCh38.d1.vd1.fa" 
REF_DICT="/gatk/project/GRCh38.d1.vd1.dict"
REF_FAI="/gatk/project/GRCh38.d1.vd1.fa.fai"


#  VCF files 
KNOWN_SITES="/gatk/project/common_all_20180418.vcf"  

# Output table
OUTPUT_TABLE="$OUTPUT_DIR/recal_data.table"  

# Function to log the current time and message
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}   


:'
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
'

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



#!/bin/bash

# Set GATK options
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
FINAL_OUTPUT_BAM="/gatk/data/output/marked_duplicates.bam"  # Existing output from MarkDuplicates
SPLIT_OUTPUT_BAM="/gatk/data/output/split_ncigar_reads.bam"  # Final output file after SplitNCigarReads
REF_GENOME="/gatk/project/GRCh38.p13.genome.fa"  # Reference genome

# Use SplitNCigarReads with stack trace enabled
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SplitNCigarReads \
    -R "$REF_GENOME" \  # Use the correct reference genome path
    -I "$FINAL_OUTPUT_BAM" \
    -O "$SPLIT_OUTPUT_BAM"  # Output the result to a new file

# Check if SplitNCigarReads completed successfully
if [ $? -eq 0 ]; then
    echo "SplitNCigarReads completed successfully."
else
    echo "SplitNCigarReads failed."
    exit 1
fi

# Final success message
echo "Processing completed successfully."


#!/bin/bash

# Set GATK options
GATK_JAR="/gatk/gatk-package-4.6.0.0-local.jar"
INPUT_BAM="/gatk/source/raw_data/bef72dd6-39e3-4b46-85c3-064af1d1979a/e653f78b-82dc-4068-938e-c59cfae00c68.rna_seq.transcriptome.gdc_realn.bam"
OUTPUT_DIR="/gatk/data/output"
FINAL_OUTPUT_BAM="$OUTPUT_DIR/marked_duplicates.bam"  # Intermediate output file after MarkDuplicates
SPLIT_OUTPUT_BAM="$OUTPUT_DIR/split_ncigar_reads.bam"  # Final output file after SplitNCigarReads
METRICS_FILE="$OUTPUT_DIR/marked_dup_metrics.txt"  # Metrics file (optional)

# Set reference genome path to the new file
REF_GENOME="/gatk/project/GRCh38.p13.genome.fa"


# Sort the input BAM file
samtools sort "$INPUT_BAM" -o "$SORTED_BAM"
if [ $? -ne 0 ]; then
    echo "Error sorting BAM file."
    exit 1
fi

# Run MarkDuplicates
java -Dsamjdk.use_async_io_read_samtools=false \
    -Dsamjdk.use_async_io_write_samtools=true \
    -Dsamjdk.use_async_io_write_tribble=false \
    -Dsamjdk.compression_level=2 \
    -jar "$GATK_JAR" MarkDuplicates \
    -I "$SORTED_BAM" \
    -O "$MARKED_DUP_BAM" \
    -M "$METRICS_FILE" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true \
    --VERBOSITY DEBUG

# Check if the process completed successfully
if [ $? -eq 0 ]; then
    echo "MarkDuplicates completed successfully."
else
    echo "MarkDuplicates failed."
    exit 1
fi


# Use SplitNCigarReads with stack trace enabled
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SplitNCigarReads \
    -R "$REF_GENOME" \  # Use the correct reference genome path
    -I "$FINAL_OUTPUT_BAM" \
    -O "$SPLIT_OUTPUT_BAM"  # Output the result to a new file

# Check if SplitNCigarReads completed successfully
if [ $? -eq 0 ]; then
    echo "SplitNCigarReads completed successfully."
else
    echo "SplitNCigarReads failed."
    exit 1
fi

# Final success message
echo "Processing completed successfully."


# Data Cleanup
echo "Running MarkDuplicates..."
gatk MarkDuplicates \
    I=$BAM_FILE \
    O="/gatk/data/output/output.bam" \
    M="/gatk/data/output/marked_dup_metrics.txt"\
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true


# Data Cleanup
 echo "Running MergeBamAlignment..."
 java -jar /gatk/picard.jar MergeBamAlignment \
     R=$REFERENCE \
     I=unaligned.bam \
     P=649e5ee8-e10c-4268-a1d7-77b38f5154ea.rna_seq.transcriptome.gdc_realn.bam \
     O=merged.bam \
     CREATE_INDEX=true
'

:'
# Split NCigar Reads
echo "Running SplitNCigarReads..."
gatk SplitNCigarReads \
    -R $REFERENCE \
    -I input.bam \
    -O output.bam
    
# Base Quality Recalibration
echo "Running BaseRecalibrator..."
gatk BaseRecalibrator \
    -I my_reads.bam \
    -R $REFERENCE \
    --known-sites sites_of_variation.vcf \
    --known-sites another/optional/setOfSitesToMask.vcf \
    -O recal_data.table

echo "Applying BQSR..."
gatk ApplyBQSR \
    -R $REFERENCE \
    -I input.bam \
    --bqsr-recal-file recal_data.table \
    -O output.bam

# Variant Calling
echo "Running Mutect2..."
gatk Mutect2 \
    -R $REFERENCE \
    -I tumor.bam \
    -I normal.bam \
    -normal normal_sample_name \
    --germline-resource af-only-gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O somatic.vcf.gz

# Calculate contamination
echo "Running GetPileupSummaries..."
gatk GetPileupSummaries \
    -I tumor.bam \
    -V common_biallelic.vcf.gz \
    -L common_biallelic.vcf.gz \
    -O pileups.table

echo "Running CalculateContamination..."
gatk CalculateContamination \
    -I pileups.table \
    -O contamination.table

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


OUT_PUT_BAM= "/gatk/data/output/output.bam"
'
#this wil run if the code succed
echo "Pipeline finished!"

