
sudo docker ps # check running containers 
sudo docker ps -a # check docker all containers 
docker inspect <container_id_or_name> # Get Detailed Information
exit # How to exit a container 
docker start <container_name_or_id> # Step 1: Start the Container
sudo docker exec -it adoring_carson /bin/bash  #execute an already excisting container

#view script 
cat /gatk/project/run_pipeline.sh

# to update version outside my container
scp /Users/qluelda/Desktop/Master_thesis/run_pipeline.sh qluelda@bender:~/project/
scp /Users/qluelda/Desktop/Master_thesis/final_mutect2.sh qluelda@bender:~/project/
scp /Users/qluelda/Desktop/elin_bin.sh qluelda@bender:~/project/
rsync -avz /Users/qluelda/Desktop/Master_thesis/run_pipeline.sh qluelda@bender:~/project/

# run a container 
sudo docker run -v ~/project:/gatk/project \
               -v /source:/gatk/source \
               -v /data:/gatk/data \
               -it broadinstitute/gatk:4.6.0.0

# curl to use 
sudo curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.12.3/picard.jar
# dowload file with curl
curl -L -o GRCh38.d1.vd1_GATK_indices.tar.gz "https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5"

# make sure a script i executable
chmod +x run_pipeline.sh
chmod +x final_mutect2.sh
chmod +x elin_bin.sh

# run the script in the container 
sed -i 's/\r$//' run_pipeline.sh # to run on windows
:'
./gatk/project/run_pipeline.sh
./run_pipeline.sh
./final_mutect2.sh
'


# -----COLORS of files in Terminal-------
# Executable file in container : -rwxr-xr-x (green)
# Directory (blue)
# not executable file in container : -rw-r--r-- (black)

# Remove all files that looks like this, whithout askking before files deleted
rm -f samtools.846.67.tmp.*.bam
rm -f * # alla i den mappen :) 

# Reference genome

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

# dowload a panel of normals 
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz /gatk/project/

# index the .vcf file for PoN 
gatk IndexFeatureFile \
    -I /gatk/project/1000g_pon.hg38.vcf.gz

# download the germline-resource file
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz /gatk/project/
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi /gatk/project/
GERMLINE_RESOURCE="/gatk/project/af-only-gnomad.hg38.vcf.gz"

# dowload common somatic variant sites .vcf
gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz /gatk/project/

# Get file name for segmentation in filter Mutec Calls 
samtools view - 

#<<<<<<<< start a screen session>>>>>>>>>
# keep processes running on a remote server even if you lose your connection or need to close the terminal
screen -S session_name
Ctrl + A, then D # Detach from a Screen Session (Keep it Running in the Background)
screen -r session_name #3. Reattach to a Screen Session
screen -ls #List All Active Screen Sessions
screen -S session_name -X quit #kill a screen session
