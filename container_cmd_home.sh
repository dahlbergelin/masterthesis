sudo docker ps # check running containers 
sudo docker ps -a # check docker all containers 
docker inspect <container_id_or_name> # Get Detailed Information
exit # How to exit a container 
docker start <container_name_or_id> # Step 1: Start the Container
sudo docker exec -it adoring_carson /bin/bash  #execute an already excisting container

# run a container 
sudo docker run -v ~/project:/gatk/project \
               -v /source:/gatk/source \
               -v /data:/gatk/data \
               -it broadinstitute/gatk:4.6.0.0

/gatk/project/sample_id_table/sample_id_table.txt
# to update version outside my container
scp /Users/qluelda/Desktop/Master_thesis/run_pipeline.sh qluelda@bender:~/my_scripts
scp /Users/qluelda/Desktop/Master_thesis/final_mutect2.sh qluelda@bender:~/my_scripts
scp /Users/qluelda/Desktop/Master_thesis/final_mutect2.sh qluelda@bender:~/my_scripts
scp /Users/qluelda/Desktop/Master_thesis/final_vardict.sh qluelda@bender:~/my_scripts
scp /Users/qluelda/Desktop/Master_thesis/target_genes.bed qluelda@bender:~/project/
scp /Users/qluelda/Desktop/masterthesis/vardict_test.sh qluelda@bender:~/my_scripts
scp /Users/qluelda/Desktop/elin_bin.sh qluelda@bender:~/project/my_scripts

# run the script in the container 
# to run on windows
sed -i 's/\r$//' run_pipeline.sh 
sed -i 's/\r$//' elin_bin.sh
sed -i 's/\r$//' vardict_test.sh
:'
./gatk/project/run_pipeline.sh
./run_pipeline.sh
./final_mutect2.sh
'
# curl to use 
sudo curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.12.3/picard.jar
# dowload file with curl
curl -L -o GRCh38.d1.vd1_GATK_indices.tar.gz "https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5"

# make sure a script i executable
chmod +x run_pipeline.sh
chmod +x final_mutect2.sh
chmod +x elin_bin.sh
chmod +x collectSM.sh
chmod +x vardict_test.sh

# <<<<<<<< VIEW RESULUTS >>>>>>>>>
    grep -v '^#' filtered_without_seg_5158a03.vcf.gz | cut -f 1-7 | head -n 10
    grep -v '^#' funcotated_output_5158a03.vcf | cut -f 1-8 | head -n 10
    gzcat  HG38_f1r1.tar.gz | head
    gzcat  somatic_variants.vcf.gz | less
head -n 10 /gatk/data/elin_output/contamination.table # get the ten first heads in the file
cat /gatk/project/run_pipeline.sh #view script
zless /gatk/data//elin_output/filtered_without_seg.vcf.gz # view a xipped file. press key: q to exit 
bcftools view -f .,PASS -H -Ov filtered_without_seg.vcf.gz | less -S # view filtered variants n .vcf file

# Remove all files that looks like this, whithout askking before files deleted
rm -f samtools.846.67.tmp.*.bam
rm -f * # alla i den mappen :) 

# -----COLORS of files in Terminal-------
# Executable file in container : -rwxr-xr-x (green)
# Directory (blue)
# not executable file in container : -rw-r--r-- (black)


#<<<<<<<< start a screen session>>>>>>>>>
# keep processes running on a remote server even if you lose your connection or need to close the terminal
screen -S session_name
Ctrl + A, then D # Detach from a Screen Session (Keep it Running in the Background)
screen -r session_name #3. Reattach to a Screen Session
screen -ls #List All Active Screen Sessions
screen -S session_name -X quit #kill a screen session

# Baserecalibration
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
# index the .vcf file 
gatk IndexFeatureFile -I /gatk/project/common_all_20180418.vcf # output: /gatk/project/common_all_20180418.vcf.idx
# update the .vcf file with the right names 1 -> chr1
awk '{if($1 ~ /^[0-9]+$/ || $1 ~ /^X$/ || $1 ~ /^Y$/) $1="chr"$1; print}' common_all_20180418.vcf | \
bgzip -c > updated_common_all_20180418.vcf.gz
tabix -p vcf /gatk/project/updated_common_all_20180418.vcf.gz # index the new file
# now update the script with updated_common_all_20180418.vcf.gz

# Get a .GTF file for the specific intervals
curl -L -o gencode.v47.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
# convert it to a .BED file
zcat gencode.v47.annotation.gtf | awk '$3 == "exon"' | awk '{print $1"\t"$4-1"\t"$5"\t"$9}' > gencode.v47.annotation.exons.bed
# Convert BED to Interval List using Picard's BedToIntervalList tool
java -jar picard.jar BedToIntervalList \
    I=gencode.v47.annotation.exons.bed \
    O=hg38_exons.interval_list \
    SD=GRCh38.d1.vd1.dict

gatk VariantsToTable \
    -v "$FUNCOTATED_VCF" -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O "/gatk/data/elin_output/${ID_PREFIX}/output-funct.table"

# <<<<<<<<   WORKFLOW MUTECT2   >>>>>>>>>
# Reference genome - dowload file with curl
curl -L -o GRCh38.d1.vd1_GATK_indices.tar.gz "https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5"
REF_GENOME="/gatk/project/GRCh38.d1.vd1.fa"
REF_DICT="/gatk/project/GRCh38.d1.vd1.dict"
REF_FAI="/gatk/project/GRCh38.d1.vd1.fa.fai"

# Baserecalibration
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
# Youtube - somatic pipeline with Mutect2 
# download known sites files for BQSR from GATK resource bundle
wget -P ~/gatk/project/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/gatk/project/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# MUTECT2 Variant calling: dowload a panel of normals 
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz /gatk/project/
PON_VCF="/gatk/project/1000g_pon.hg38.vcf.gz"  # Mutect2: Panel of Normals 1000 genomes

# index the .vcf file for PoN 
gatk IndexFeatureFile \
    -I /gatk/project/1000g_pon.hg38.vcf.gz

# download the germline-resource file
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz /gatk/project/
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi /gatk/project/
GERMLINE_RESOURCE="/gatk/project/af-only-gnomad.hg38.vcf.gz"

# GetPileUpSummeries: dowload common somatic variant sites .vcf
gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz /gatk/project/
CONT_VCF="/gatk/project/small_exac_common_3.hg38.vcf.gz" 

# Get file name for segmentation in filter Mutec Calls 
samtools view - H 
# Get the RG information from the .BAM file 
samtools view -H Your_BAM_File.bam | grep '^@RG'

# FUNCOTATOR DATA SOURCES
# get the latest version for hg:38
FUNCOTATOR_DATA_SOURCES="/gatk/project/funcotator_dataSources.v1.8.hg38.20230908s"
https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
# Unzip the file 
sudo tar -xzvf funcotator_dataSources.v1.7.20200521s.tar.gz -C /gatk/project/
# unzip the gnomad files in the directory 
sudo tar -xzvf gnomAD_exome.tar.gz
sudo tar -xzvf gnomAD_genome.tar.gz

# view how many mutations that passed the filter
view -f PASS filtered_without_seg_5158a03.vcf.gz | grep -v '^#' | wc -l
bcftools view -f PASS filtered_without_seg_5158a03.vcf.gz | grep -v '^#' | head -n 10
bcftools view -f PASS funcotated_output_5158a03.vcf | grep -v '^#' | head