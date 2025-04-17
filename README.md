![ProteoSeeker Logo](images/ProteoSeeker_Logo.png)

# 1. Installation
## 1.1 Source code
### 1.1.1 Prerequisites
#### Anaconda
<p align="justify">To install ProteoSeeker from source code, conda, from Anaconda, must be installed and activated in your system. Instructions for the installation of Anaconda in Linux are provided through the following link:</p>

Anaconda for Linux: [https://www.anaconda.com/docs/getting-started/anaconda/install](https://www.anaconda.com/docs/getting-started/anaconda/install)

#### git
<p align="justify">Necessary to download the ProteoSeeker repository.</p>

### 1.1.2 Installation process
<p align="justify">Execute the commands below to perform the following installation steps. The installation process requires approximately <strong>20 GB</strong> of disk space and <strong>30 minutes</strong> to be completed.</p>

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

### 1.1.3 Removing installation environments, files and directories
<p align="justify">To remove the environments, all their files and the directories that were created during the installation of ProteoSeeker (by running the "./install.sh" script), run the Bash script below, in the installation directory from the installation directory. You can then delete the main directory of ProteoSeeker and all environments and files associated with installing ProteoSeeker will have been removed by your system. The Bash script below will also remove the environment ("ps_result_analysis") created by running the commands provided at section "4.1" and used to run the analysis of the taxonomy evaluation results.</p>

~~~bash
./remove.sh
~~~

### 1.1.4 Parameter files
<p align="justify">You can create a set of "template" parameter files which can be used to run the seek or the taxonomy mode or both modes of ProteoSeeker by running the following script in the installation directory from the same directory. This set of files is generated in the "parameter_files" directory of ProteoSeeker.</p>

~~~bash
./parameter_files.sh
~~~

<p align="justify">In addition, to facilitate further the process of running ProteoSeeker, two parameter files, "par_DRR083188_sra" and "par_DRR083188_run", are already present in the "parameter_files" directory. The latter files can easily be modified to run ProteoSeeker based on the specifications of your system and analysis.</p>

### 1.1.5 COMEBin - GPU
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

## 1.2 Docker
<p align="justify">To install ProteoSeeker from Docker Hub as a Docker image, Docker must be installed in your system. To install Docker in Ubuntu, follow the instructions provided by the link below. Then, download the image of ProteoSeeker from Docker Hub. There are two versions.</p>

Docker engine for Ubuntu: https://docs.docker.com/engine/install/ubuntu/

The main_v1.0.0 version has a download size of  **13.16 GB** and decompressed has a size of **30.4 GB**:
~~~bash
sudo docker image pull skretaslab/proteoseeker
or
sudo docker image pull skretaslab/proteoseeker:latest
or
sudo docker image pull skretaslab/proteoseeker:main_v1.0.0
~~~

The light_v1.0.0 version has a download size of  **7.66 GB** and decompressed has a size of **22.3 GB**:
~~~bash
sudo docker image pull skretaslab/proteoseeker:light_v1.0.0
~~~

## 1.3 Phobius
<p align="justify">For either case of installation process, in order to use the topology and signal peptide predictions provided by Phobius you must download Phobius from https://phobius.sbc.su.se/data.html. As described in section "3.4" to utilize Phobius when running ProteoSeeker through the command-line you should also provide the path to the Phobius directory in the parameter file or as a parameter through the corresponding option of "proteoseeker.py". The default path for the Phobius installation in a Docker container from the proteoseeker Docker image is already set to the phobius directory of the shared directory and you should download and copy the Phobius installation files in that directory as explained in section "3.3" below. In any other case, ProteoSeeker will run without performing topology and signal peptide predictions in its seek functionality.</p>

# 2. Run ProteoSeeker
<p align="justify"> For detailed information regarding the use, options and parameter files of ProteoSeeker check its Wiki page.</p>
   
<p align="justify">In general, the most basic options in a parameters file that control the main functionalities of ProteoSeeker are:
- seek_mode=True/False
   The mode is applied when True and not applied when False.
- taxonomy_mode=True/False
   The mode is applied when True and not applied when False.
- seek_route=1/2/3
   Seek route 1 is based on screening each protein by searching for at least one domain corresponding to a profile from the selected protein families. Seek route 2 is based on screening each protein against a filtered protein database by searching for at least one hit with an e-value higher than a specific threshold. The filtered protein database is based on filtering the provided protein database based on the protein names corresponding to the selected protein families.
- taxonomy_route=1/2
   Taxonomy route 1 is the Kraken2 taxonomy route. Taxonomy route 2 is the COMEBin/MetaBinner taxonomy route.
</p>
   
<p align="justify">Examples of running ProteoSeeker can be based on the template files "par_DRR083188.txt", "par_DRR083188_sra.txt" and "par_DRR083188_run.txt" which are already present in the "parameter_files" directory. To use these files you should change certain paths based on the specifications of your system. These paths are described below. It should not be necessary to change the relative paths when running ProteoSeeker from the main directory.</p>

<p align="justify">
   
- protein_db_path: Can be left empty if not applying the seek mode with the seek routes (analysis types) 2 or 3 nor applying the taxonomy mode with the taxonomy route 2 (COMEBin/MetaBinner taxonomoy route).
- kraken_db_path: Can be left empty if not applying the taxonomy mode with the taxonomoy route 1 (Kraken taxonomy route).
- metabinner_bin_path: Can be left empty if not applying binning tool 1 (MetaBinner) in the taxonomy route 2 of the taxonomy mode.
- comebin_bin_path: Can be left empty if not applying binning tool 2 (COMEBin) in the taxonomy route 2 of the taxonomy mode.
- profiles_broad_path, swissprot_path, fraggenescanrs_path, conda_bin, conda_sh: Should already be set properly.

</p>

~~~bash
# Relative paths:
protein_db_path="" or protein_db_path="other/swissprot_proteins_13_04_2025/uniprot_sprot.fasta"
kraken_db_path="ps_tools/kraken2/kraken2_databases/kraken2_8st_db"
profiles_broad_path="pfam_database/Pfam-A.hmm"
swissprot_path="swissprot_database/swissprot"
fraggenescanrs_path="ps_tools/fgsrs/FragGeneScanRs"

# Full paths:
conda_bin="/home/gfilis/anaconda3_2024_02_1"
conda_sh="/home/gfilis/anaconda3_2024_02_1/etc/profile.d/conda.sh"
metabinner_bin_path="/home/gfilis/anaconda3_2024_02_1/envs/ps_metabinner/bin"
comebin_bin_path="/home/gfilis/anaconda3_2024_02_1/envs/ps_comebin/bin/COMEBin"
~~~

<p align="justify">You can use these files by running ProteoSeeker from the main directory as shown below.</p>

<p align="justify">
   
- par_DRR083188_run.txt: ProteoSeeker downloads and processes the SRA sample DRR083188. ProteoSeeker applies the seek mode and seek route 1 and the taxonomy mode and taxonomy route 1 (Kraken2) on sample DRR083188.
- par_DRR083188_sra.txt: ProteoSeeker downloads and processes the SRA sample DRR083188 (it does not apply any kind of analysis).
- par_DRR083188_protein_database.txt: ProteoSeeker downloads and processes the SRA sample DRR083188. ProteoSeeker applies the seek mode and seek route 3 and the taxonomy mode and taxonomy route 2 (COMEBin/MetaBinner) on sample DRR083188.
- par_DRR083188_tool_params.txt: ProteoSeeker downloads and processes the SRA sample DRR083188. ProteoSeeker applies the seek mode and seek route 1 and the taxonomy mode and taxonomy route 1 (Kraken2). Custom parameter values have been set for BBDuk and Megahit.

</p>

~~~bash
conda activate ps_env
python proteoseeker -pfp parameter_files/par_DRR083188_run.txt
python proteoseeker -pfp parameter_files/par_DRR083188_sra.txt
python proteoseeker -pfp parameter_files/par_DRR083188_protein_database.txt
python proteoseeker -pfp parameter_files/par_DRR083188_tool_params.txt
conda deactivate
~~~

# Wiki
<p align="justify">Extensive documentation on the functionality and usage of ProteoSeeker is available in its [<ins>Wiki</ins>](https://github.com/SkretasLab/ProteoSeeker/wiki).</p>

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
