![ProteoSeeker Logo](images/ProteoSeeker_Logo.png)

<br>
<br>

## [<ins>Usage Details</ins>](./README.md)

<br>
<br>

# 1. Installation
## 1.1 Docker
<p align="justify">To install ProteoSeeker from Docker Hub as a Docker image, Docker must be installed in your system. To install Docker in Ubuntu, follow the instructions provided by the link below:</p>

Docker engine for Ubuntu: https://docs.docker.com/engine/install/ubuntu/

<p align="justify">Then, download the image of ProteoSeeker from Docker Hub. There are two versions.</p>

The main_v1.0.0 version has a download size of  **13.16 GB** and decompressed has a size of **30.4 GB**. To install the main_v1.0.0 version use **one** of the following commands:
~~~bash
sudo docker image pull skretaslab/proteoseeker
or
sudo docker image pull skretaslab/proteoseeker:latest
or
sudo docker image pull skretaslab/proteoseeker:main_v1.0.0
~~~

The light_v1.0.0 version has a download size of  **7.66 GB** and decompressed has a size of **22.3 GB**. To install the light_v1.0.0 version use the following command:
~~~bash
sudo docker image pull skretaslab/proteoseeker:light_v1.0.0
~~~

## 2.2 Source code
### 2.2.1 Prerequisites
#### Anaconda
<p align="justify">To install ProteoSeeker from source code, conda, from Anaconda, must be installed and activated in your system. Instructions for the installation of Anaconda in Linux are provided through the following link:</p>

Anaconda for Linux: https://docs.anaconda.com/free/anaconda/install/linux/

#### git
<p align="justify">Necessary to download the ProteoSeeker repository.</p>

### 2.2.2 Installation
<p align="justify">Execute the commands below to perform the following installation steps.</p>

~~~bash
git clone https://github.com/SkretasLab/ProteoSeeker.git
cd ProteoSeeker-main
sudo chmod -R 755 installation
cd installation
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible
./install.sh
~~~

### 2.2.3 Parameter files
<p align="justify">You can create a set of "template" parameter files which can be used to run the seek or the taxonomy mode or both modes of ProteoSeeker by running the following script in the installation directory from the same directory. This set of files is generated in the main directory of ProteoSeeker.</p>

~~~bash
./parameter_files.sh
~~~

### 2.2.4 COMEBin - GPU
<p align="justify">It should be noted that COMEBin can also be installed and run on a GPU. Instructions are available at: https://github.com/ziyewang/COMEBin and also below:</p>

~~~bash
conda create -n ps_comebin_gpu
conda activate ps_comebin_gpu
conda install -c conda-forge -c bioconda comebin
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia -c conda-forge
conda deactivate
~~~

<p align="justify">We have observed that running COMEBin on a GPU offers a great improvement in the running time of COMEBin. To use COMEBin in an environment which allows the usage of a GPU, provide that environments's name to the default name of the environment for COMEBin and also modify the path to the COMEBin directory of the environment. These paths refer to the following options of ProteoSeeker:</p>

~~~bash
   -sen/--comebin-env             Str -Opt: ps_comebin- The conda environment for COMEBin. 'None/none': To not use an environment at all.

   -cfp/--comebin-folder-path     Str -Opt- The path to the parent directory of "run_comebin.sh" of COMEBin.
~~~

<p align="justify">For example, the second option in our system has the following value: "/home/compteam/anaconda3/envs/ps_comebin_gpu/bin/COMEBin"</p>

### 2.2.5 Removing installation environments, files and directories
<p align="justify">To remove the environments, all their files and the directories that were created during the installation of ProteoSeeker (by running the "./install.sh" script), run the Bash script below, in the installation directory from the installation directory. You can then delete the main directory of ProteoSeeker and all environments and files associated with installing ProteoSeeker will have been removed by your system. The Bash script below will also remove the environment ("ps_result_analysis") created by running the commands provided at section "4.1" and used to run the analysis of the taxonomy evaluation results.</p>

~~~bash
./remove.sh
~~~

## 2.3 Phobius
<p align="justify">For either case of installation process, in order to use the topology and signal peptide predictions provided by Phobius you must download Phobius from https://phobius.sbc.su.se/data.html. As described in section "3.4" to utilize Phobius when running ProteoSeeker through the command-line you should also provide the path to the Phobius directory in the parameter file or as a parameter through the corresponding option of "proteoseeker.py". The default path for the Phobius installation in a Docker container from the proteoseeker Docker image is already set to the phobius directory of the shared directory and you should download and copy the Phobius installation files in that directory as explained in section "3.3" below. In any other case, ProteoSeeker will run without performing topology and signal peptide predictions in its seek functionality.</p>

## 2.4 UniprotKB/Swiss-Prot and Pfam preprocessed datasets
<p align="justify">In case the user needs to update the preprocessed files generated from processing information from the UniprotKB/Swiss-Prot flat file or the Pfam database the scripts available in the ps_scripts can be utilized.</p>

<p align="justify">For the Swiss-Prot/UniprotKB datasets the user should download the flat file for UniprotKB/Swiss-Prot and then run the command provided below. This command will generate a series of files as output in the same directory. The files to be updated are to be replaced in the "profile_protein_dbs" directory. These files are: "prfamilies_numbered.tsv", "prfamilies_pfamdomains.tsv" and "prfamilies_length.tsv". For more information about the script the user can run the script with the "-h" option.</p>

~~~bash
python uniprot_data_acs.py -i <uniprotkb/swiss-prot_flat_file_path>
~~~

<p align="justify">For the Pfam datasets the Pfam HMM database is needed. This database should already be installed in the "pfam_database" directory of the main "proteoseeker" directory. Run the command provided below. This command will generate a series of files as output in the same directory. The files to be updated are to be replaced in the "profile_protein_dbs" directory. These files are: "pfam_accs_names.tsv" and "profiles_lengths.tsv". For more information about the script the user can run the script with the "-h" option.</p>

~~~bash
python pfam_info.py -i <pfam_hmm_file_path>
~~~

# Publication
[https://doi.org/10.1002/advs.202414877](https://doi.org/10.1002/advs.202414877)

# Citation
G. Filis, D. Bezantakou, K. Rigkos, D. Noti, P. Saridis, D. Zarafeta, G. Skretas, ProteoSeeker: A Feature-Rich Metagenomic Analysis Tool for Accessible and Comprehensive Metagenomic Exploration. Adv. Sci. 2025, 2414877. https://doi.org/10.1002/advs.202414877

# Website
<p align="justify">General information about ProteoSeeker, including part of the installation and usage instructions found in this repository, can be found at the link below:</p>

https://skretaslab.gr/proteoseeker

# Contacts and bug reports
Feel free to send questions or bug reports by using one of the following email addresses:
1. ProteoSeeker Team: proteoseeker@fleming.gr
2. Georgios Filis: filis@fleming.gr
3. Dimitra Bezantakou: bezantakou@fleming.gr
4. Dimitra Zarafeta: zarafeta@fleming.gr
5. Georgios Skretas: skretas@fleming.gr

# License
This project is licensed under the GPLv3 License. See the [LICENSE](LICENSE) file for details.
