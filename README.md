# 1. Overiview
ProteoSeeker is feature-rich metagenomic analysis tool which allows for the identification of novel proteins originating from metagenomes. These proteins may be part of specified protein families or/and may be subjected to a taxonomy analysis.

# 2. Installation
## 2.1 Source code
To install ProteoSeeker from source code, conda (from Anaconda or Miniconda) must be installed and activated in your system. Instructions for installing Anaconda and Miniconda in Linux are provided in the following links:
Anaconda for Linux: https://docs.anaconda.com/free/anaconda/install/linux/
Miniconda for Linux: https://docs.anaconda.com/free/miniconda/miniconda-install/

Then download the repository, extract it, move to the installation folder and run the installation script. Follow the steps below:
~~~bash
git clone https://github.com/SkretasLab/ProteoSeeker.git
unzip ProteoSeeker.zip
cd ProteoSeeker/installation
bash instal.sh
~~~

## 2.2 Docker
To install ProteoSeeker from DockerHub as a docker image docker must be installed in your system. To install docker in Linux follow the instructions provided by the link below
Docker engine for Linux: https://docs.docker.com/engine/install/ubuntu/

Then simply download the image from dockerhub:
~~~bash
pull proteoseeker
~~~

### Phobius
For either case of installation process followed in order to use the topology and signal peptide predictions provided by Phobius you must download Phobius from .... If not downloaded ProteoSeeker will run normaly without including topology and signal peptide predictions in its seek functionality.

# 3. Use


