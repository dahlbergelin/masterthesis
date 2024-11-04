#<<<<<<<<<<<<<<< INSTALL CONDA >>>>>>>>>>>
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#Run the installer:
bash Miniconda3-latest-Linux-x86_64.sh
# done --> restart the terminal 
conda --version # see that it was correctly installed

#<<<<<<<<<<<<<<< INSTALL VARDCIT >>>>>>>>>>>
# STEP 1: create an enviroment 
conda create --name vardict_env #environment location: /home/qluelda/miniconda3/envs/vardict_env
# activate envrioment
conda activate vardict_env
# Step 2: Add Bioconda Channel 
conda config --add channels bioconda
# - add the conda-forge channel:
conda config --add channels conda-forge
# Step 3: Install vardict-java
conda install vardict-java
# verify that Vardict was installed 
vardict-java -h

# ex on running vardcict:
vardict-java -G reference.fa -f 0.01 -N sample_name -b input.bam -R region.bed > output.vardict

# Deactivate the Environment (If Needed)
# When you're done working in this environment, you can deactivate it with:
conda deactivate