# 1. Overiview
ProteoSeeker is feature-rich metagenomic analysis tool which allows for the identification of novel proteins originating from metagenomes. These proteins may be part of specified protein families or/and may be subjected to a taxonomy analysis.

![ProteoSeeker Overview](images/Figure_1.png)

1. Specific characteristics for the environment of sample collection are documented.
2. The sample is collected by the environment.
3. The DNA of the sample is isolated and prepared for sequencing.
4. The DNA is sequenced through Next-Generation Sequencing (NGS).
5. Files containing reads are the result of NGS.
6. Such files can be provided directly to ProteoSeeker for analysis. Alternatively, these files may be uploaded in an open-access database alongside their metadata. The data of such databases form the exploration ground for ProteoSeeker in discovering novel enzymes functioning in predefined conditions enriching the scientific community's capacity to explore microbial ecosystems. A user may download such data and provide it to ProteoSeeker or in the case of SRA of NCBI use the code (SRA accession) of the sample from the online database directly to ProteoSeeker.
7. ProteoSeeker identifies proteins originating from the input reads.
8. ProteoSeeker offers two main functionalities. The first one is the “seek” functionality. Certain proteins are grouped in protein families, initially selected by the user, based on profiles from the Pfam database. This functionality is used to discover novel enzymes with specific functionality able to perform under predefined conditions, as documented for the sample analyzed.
9. The second one is the “taxonomy” functionality. One or more organisms are assigned to certain proteins.*

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


