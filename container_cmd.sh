# check running containers 
sudo docker ps

# check docker all containers 
sudo docker ps -a

# Get Detailed Information
docker inspect <container_id_or_name>

# How to exit a container 
exit

# curl to use 
sudo curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.12.3/picard.jar

# Step 1: Start the Container
docker start <container_name_or_id>

# dowload file with curl
curl -L -o GRCh38.d1.vd1_GATK_indices.tar.gz "https://api.gdc.cancer.gov/data/2c5730fb-0909-4e2a-8a7a-c9a7f8b2dad5"

#execute an already excisting container
sudo docker exec -it adoring_carson /bin/bash

# make sure a script i executable
chmod +x run_pipeline.sh
chmod +x final_mutect2.sh

# run the script
:'
./gatk/project/run_pipeline.sh
./run_pipeline.sh
./final_mutect2.sh
'

#view script 
cat /gatk/project/run_pipeline.sh


# to update version outside my container
scp /Users/qluelda/Desktop/Master_thesis/run_pipeline.sh qluelda@bender:~/project/
scp /Users/qluelda/Desktop/Master_thesis/final_mutect2.sh qluelda@bender:~/project/
rsync -avz /Users/qluelda/Desktop/Master_thesis/run_pipeline.sh qluelda@bender:~/project/

# run a container 
sudo docker run -v ~/project:/gatk/project \
               -v /source:/gatk/source \
               -v /data:/gatk/data \
               -it broadinstitute/gatk:4.6.0.0

# -----COLORS of files in Terminal-------
# Executable file in container : -rwxr-xr-x (green)
# Directory (blue)
# not executable file in container : -rw-r--r-- (black)

# Remove all files that looks like this, whithout askking before files deleted
rm -f samtools.846.67.tmp.*.bam
