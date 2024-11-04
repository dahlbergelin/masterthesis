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

# # GetPileUpSummeries: dowload common somatic variant sites .vcf
gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz /gatk/project/
CONT_VCF="/gatk/project/small_exac_common_3.hg38.vcf.gz" 

# Get file name for segmentation in filter Mutec Calls 
samtools view - H 
# get the RG information from the .BAM file 
samtools view -H Your_BAM_File.bam | grep '^@RG'

# Fucotator data sources
# get the latest version for hg:38
FUNCOTATOR_DATA_SOURCES="/gatk/project/funcotator_dataSources.v1.8.hg38.20230908s"
https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
# Unzip the file 
sudo tar -xzvf funcotator_dataSources.v1.7.20200521s.tar.gz -C /gatk/project/
# unzip the gnomad files in the directory 
sudo tar -xzvf gnomAD_exome.tar.gz
sudo tar -xzvf gnomAD_genome.tar.gz

