![ProteoSeeker Logo](images/ProteoSeeker_Logo.png)

<br>
<br>

## [<ins>Usage Details</ins>](./README.md)

<br>
<br>

# 1. Installation
## 1.1 Source code
### 1.1.1 Prerequisites
#### Anaconda
<p align="justify">To install ProteoSeeker from source code, conda, from Anaconda, must be installed and activated in your system. Instructions for the installation of Anaconda in Linux are provided through the following link:</p>

Anaconda for Linux: https://docs.anaconda.com/free/anaconda/install/linux/

#### git
<p align="justify">Necessary to download the ProteoSeeker repository.</p>

### 1.1.2 Installation process
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

### 1.1.3 Removing installation environments, files and directories
<p align="justify">To remove the environments, all their files and the directories that were created during the installation of ProteoSeeker (by running the "./install.sh" script), run the Bash script below, in the installation directory from the installation directory. You can then delete the main directory of ProteoSeeker and all environments and files associated with installing ProteoSeeker will have been removed by your system. The Bash script below will also remove the environment ("ps_result_analysis") created by running the commands provided at section "4.1" and used to run the analysis of the taxonomy evaluation results.</p>

~~~bash
./remove.sh
~~~

### 1.2 Parameter files
<p align="justify">You can create a set of "template" parameter files which can be used to run the seek or the taxonomy mode or both modes of ProteoSeeker by running the following script in the installation directory from the same directory. This set of files is generated in the main directory of ProteoSeeker.</p>

~~~bash
./parameter_files.sh
~~~

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

### 1.3 COMEBin - GPU
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

## 1.4 Phobius
<p align="justify">For either case of installation process, in order to use the topology and signal peptide predictions provided by Phobius you must download Phobius from https://phobius.sbc.su.se/data.html. As described in section "3.4" to utilize Phobius when running ProteoSeeker through the command-line you should also provide the path to the Phobius directory in the parameter file or as a parameter through the corresponding option of "proteoseeker.py". The default path for the Phobius installation in a Docker container from the proteoseeker Docker image is already set to the phobius directory of the shared directory and you should download and copy the Phobius installation files in that directory as explained in section "3.3" below. In any other case, ProteoSeeker will run without performing topology and signal peptide predictions in its seek functionality.</p>

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
