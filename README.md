# 1. Overiview
ProteoSeeker is a feature-rich metagenomic analysis tool which allows for the identification of novel proteins originating from metagenomes. These proteins may be part of specified protein families or/and may be subjected to a taxonomy analysis.

![ProteoSeeker Overview](images/Figure_1.png)

1. Sample Documentation: Specific characteristics of the sample's environmental source, including factors such as location, habitat, sampling conditions and collection method are documented.
2. Sample Collection: The metagenomic material is collected from the enviromnental nieche of interest
3. DNA Isolation and Preparation: Following DNA isolation, DNA extraction of the collected sample is performed and preped for sequencing.
4. Next-Generation Sequencing (NGS):  NGS technology is performed to collect the metagenomic dataset of the sample
5. NGS Data Processing: The sequence files (reads) resulting from NGS are generated and  quality control is applied to the data.
6. Data Sharing and Accessibility: Datasets and metadata are shared in open-access databases, facilitating collaborative research and data reuse. Such files can be provided directly to ProteoSeeker for analysis, forming the exploration ground for the tool. Proteoseeker aspires to provide a coprehensive, user-friendly platform for the discovery of novel proteins/enzymes originating from envoronmetns of interest, enriching the scientific community's capacity to explore microbial ecosystems. A user may download such data and provide it to ProteoSeeker or in the case of SRA of NCBI use the code (SRA accession) of the sample from the online database directly to ProteoSeeker.
7. ProteoSeeker Analysis: The tool identifies putative proteins derived from the input reads.
8. Functional Analysis with ProteoSeeker: Detailed functionalities offered by ProteoSeeker, includes "seek" and "taxonomy" functionalities, and their respective purposes in protein/enzyme discovery and taxonomic assignment.
9. Protein Family Profiling: Protein family profiles from databases like Pfam groups proteins and facilitate the discovery of novel proteins/enzymes with specific functionalities.
10. Taxonomic Assignment: The tool expands on the process of assigning one or more organisms to identified proteins, aiding in the understanding of microbial community composition.


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


