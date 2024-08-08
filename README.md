# 1. Overview
## 1.1 Scope of usage
<p align="justify">ProteoSeeker is a feature-rich metagenomic analysis tool for accessible and comprehensive metagenomic exploration. It allows for the identification of novel proteins originating from metagenomes. These proteins may be part of specified protein families or/and may be subjected to a taxonomy analysis.</p>

![ProteoSeeker Overview](images/Figure_1.png)

1. **Sampling Site Documentation:** Specific characteristics of the sample's environmental source, including factors such as location, habitat, sampling conditions and collection method are documented.
2. **Sample Collection:** The metagenomic material is collected from the enviromnental nieche of interest.
3. **DNA Isolation and Preparation:** Following DNA extraction, the metagenomic material is preped for sequencing.
4. **Next-Generation Sequencing (NGS):**  NGS is performed to collect the metagenomic dataset of the sample
5. **NGS Data Processing:** The sequencing files (reads) resulting from NGS are generated and data quality control is performed. Datasets and metadata are shared in open-access databases, facilitating collaborative research and data reuse.
   
<p align="justify">*Such files can be provided directly to ProteoSeeker for analysis, forming the exploration ground for the tool. Proteoseeker aspires to provide a coprehensive, user-friendly platform for the discovery of novel proteins/enzymes originating from envoronmetns of interest, enriching the scientific community's capacity to explore microbial ecosystems. A user may download such data and provide it to ProteoSeeker or in the case of SRA of NCBI use the code (SRA accession) of the sample from the online database directly to ProteoSeeker.*</p>

6. **ProteoSeeker Analysis:** The selected dataset is uploaded to Proteoseeker. The tool identifies putative proteins derived from the input reads.
7. **Functional Analysis:** Functionalities offered by ProteoSeeker include "seek" and "taxonomy" functionalities, and their respective purposes in protein/enzyme discovery and taxonomic assignment.
8. **Protein Family Profiling:** Protein family profiles from databases like Pfam, groups proteins and facilitates the discovery of novel proteins/enzymes with specific functionalities.
9. **Taxonomic Assignment:** The tool expands on the process of assigning one or more organisms to identified proteins, aiding in the understanding of microbial community composition.

## 1.2 Pipeline
<p align="justify">The overall pipeline implemented by ProteoSeeker can be seen in the image below. ProteoSeeker offers two main functionalities with a multitude of options for users not accustomed to metagenomic analysis tools and more advanced users that may want to modify the behaviour of certain tools included in ProteoSeeker's pipeline. The first functionality is to **seek** proteins that may be part of selected protein families and the second funcitonality is to find the **taxonomy** of the proteins discovered from the analysis of a sample. The second functionality is based either on the taxonomy of the reads (kraken2 route) or on binning the contigs and searching for the taxonomy of the proteins through the "taxonomy filtered protein database" (COMEBin/MetaBinner route).</p>

The stages of the “seek” mode of ProteoSeeker. ProteoSeeker offers two main functionalities applied through the seek mode (blue) and the “taxonomy” mode (green). Each stage is colored based on the mode it belongs to. The possible types of input for ProteoSeeker include an SRA code, reads in FASTQ files, contigs or genomes or proteins if FASTA format. If an SRA code is provided the corresponding SRA file and FASTQ files are generated.

![ProteoSeeker Seek Functionality](images/Figure_2.png)

1. The protein families selected are determined based on their input codes.
2. The profiles and protein names associated with the selected families are collected and the spd is created.
3. The protein database is filtered based on the collected protein names the sfpd is created.
4. The reads of the FASTQ files undergo several quality control checks by FastQC.
5. The reads are preprocessed by BBDuk and then are analyzed by FastQC for a second time.
6. The preprocessed reads are assembled into contigs by Megahit.
7. Protein coding regions (pcdrs) are predicted in the contigs by FragGeneScanRs.
8. CD-HIT is used to reduce the redundancy of the pcdrs.
9. The pcdrs are screened against the spd with HMMER. Any pcdr with at least one hit based on the latter screening is retained (set 1).
10. The rest of the pcdrs are screened against the sfpd and only those with a hit of low enough e-value are retained (set 2). In addition, set 1 is screened against the Swiss-Prot protein database.
11. Both sets are screened against all the profiles of the Pfam database with HMMER.
12. Topology prediction are performed by Phobius.
13. Motifs provided by the user are screened against each protein.
14. The protein family of each protein is predicted. 15. Annotation files are written.

The stages of the “taxonomy” mode of ProteoSeeker. ProteoSeeker offers two main functionalities applied through the “seek” mode (blue) and the “taxonomy” mode (green). The taxonomic analysis can be performed by either of two”routes” of analysis. The fist one is based on binning based on the taxonomic analysis of reads from Kraken2 (orange) and the second one is performed by COMEBin or MetaBinner (purple). The possible types of input for ProteoSeeker include an SRA code, reads in FASTQ files, contigs or genomes or proteins if FASTA format. If an SRA code is provided the corresponding SRA file and FASTQ files are generated.

![ProteoSeeker Taxonomy Functionality](images/Figure_3.png)

1. The protein families selected are determined based on their input codes.
2. The profiles and protein names associated with the selected families are collected and the tpd is created.
3. The protein database is filtered based on the collected protein names the tfpd is created.
4. The reads of the FASTQ files undergo several quality control checks by FastQC.
5. The reads are preprocessed by BBDuk and then are analyzed by FastQC for a second time.
6. The preprocessed reads are assembled into contigs by Megahit.
7. Protein coding regions (pcdrs) are predicted in the contigs by FragGeneScanRs.
8. CD-HIT is used to reduce the redundancy of the pcdrs.
9. Bowtie2 maps the reads to the contigs.
10. The contigs are binned based on MetaBinner or COMEBin.
11. The pcdrs are screened against the tpd with HMMER.
12. Any pcdr with at least one hit against the tpd is screened against the tfpd.
13. Species are assigned to the proteins based on the species assignment to the bins.
14. Each bin, along with any taxa assigned to it, is quantified based on the reads mapped to its contigs.
15. Annotation files are written.
16. Species are assigned to the reads based on Kraken2.
17. Through the read-contig mapping each species is quantified for each contig. Species are assigned to the contigs.
18. The contigs are binned based on their species.
19. Species are assigned to the bins.
20. Species are assigned to the proteins of the latter bins.
21. Annotation files are written.

# 2. Installation
## 2.1 Source code
### Prerequisites
#### 1. Anaconda:
<p align="justify">To install ProteoSeeker from source code, conda (from Anaconda or Miniconda) must be installed and activated in your system. Instructions for installing Anaconda and Miniconda in Linux are provided in the following links:</p>

Anaconda for Linux: https://docs.anaconda.com/free/anaconda/install/linux/

#### 2. git:
<p align="justify">Necessary to download the ProteoSeeker repository.</p>

### Dependencies
<p align="justify">All dependencies, except for the protein database, are automatically installed by the installation process of ProteoSeeker. You can skip this part if you want to go straight to the installation instructions but do read the information related to the protein database which is not installed automatically by the installation process of ProteoSeeker.</p>

#### 1. git, wget, gzip, tar:
<p align="justify">These packages are handled by the installation process of ProteoSeeker and are installed in a enviroment from which the installation of ProteoSeeker will be initiated. A new enviroment called "ps_install" is created, used to download and collect the latter packages. If one wants to use a custom conda enviroment to base upon it the installation of ProteoSeeker, the latter enviroment should contain the following:</p>
   
~~~bash
conda install anaconda::git
conda install anaconda::wget
conda install conda-forge::gzip
conda install conda-forge::tar
~~~

<p align="justify">In addition, the name of the custom enviroment must replace the name of the enviroment ("ps_install") in line 14 ("conda activate ps_install") in the file "install.sh". It is suggested not to perform these changes if possible and let the installation process handle the creation and usage of the installation enviroment for these packages.</p>

#### 2. Pipeline tools
<p align="justify">These tools are all aumatically installed by the installation process of ProteoSeeker or have already been set in the Docker image of ProteoSeeker. The spcific versions of these tools which are during ProteoSeeker's installation are the ones also used for the evaluation of the "seek" and "taxonomy" modes of ProteoSeeker. These tools and their versions are the following:</p>

1. bbmap: 39.01
2. bowtie2: 2.5.3
3. cd-hit: 4.8.1
4. comebin: 1.0.4
5. diamond: 2.1.9
6. fastqc 0.12.1
7. hmmer: 3.4
8. kraken2: 2.1.3
9. megahit: 1.2.9
10. metabinner: 1.4.2
11. sra-tools: 3.1.0
12. taxonkit: 0.16.0
13. csvtk: 0.30.0
14. FragGeneScanRs: 1.1.0

#### 3. Databases
<p align="justify">The latests versions of the following databases are installed automatically by ProteoSeeker. To see the datatime of their collection for the evaluation of ProteoSeeker check the evaluation section. Only the protein database should be installed by the user, in which case he can use one of the already prepared shell scripts for installing either the nr or the Uniref50 or Uniref90 or Uniref100 databases. The latter files are located in the "installation" directory. The user should be certain that the system has enough available memory to hold the decompressed protein database that he wants to download.</p>

1. Pfam database: Latest - Installed automatically
2. Swiss-Prot/UniprotKB database: Latest - Installed automatically
3. GTDB taxonomy taxdump files: Latest - Installed automatically
4. Kraken 2/Bracken Refseq indexes: Collection Standard-8 of 05/06/2024 - Installed automatically
5. nr database: Latest - Installed manually

### Installation
<p align="justify">Open a terminal, download the repository, extract it and move to the installation folder. Make sure that all shell scripts in the installation directory can be executed as programs. If not, change the permissions as to allow for their execution. The conda enviroment where the installation will take place should be activate (whether it is the base or another enviroment). Then run the installation script. Depending on your internet connection it might take from 5 to 30 minutes for the installation to be completed. Follow the steps below:</p>

~~~bash
git clone https://github.com/SkretasLab/ProteoSeeker.git
cd ProteoSeeker_main
sudo chmod -R 755 installation
cd installation
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible
./instal.sh
~~~

### Parameter files
<p align="justify">You can create a set of "template" parameter files which can be used to run the "seek" or the "taxonomy" mode or both modes of ProteoSeeker based on the taxonomy route of Kraken2 or COMEBin/MetaBinner by running the following script in the installation directory from the same directory. This set of files is generated in the main directory of ProteoSeeker.</p>

~~~bash
./parameter_files.sh
~~~

### COMEBin - GPU
<p align="justify">It should be noted that COMEBin can also be installed and run on a GPU. Instructions are available at: https://github.com/ziyewang/COMEBin and also below:</p>

~~~bash
conda create -n ps_comebin_gpu
conda activate ps_comebin_gpu
conda install -c conda-forge -c bioconda comebin
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia -c conda-forge
conda deactivate
~~~

<p align="justify">We have observed that running COMEBin with a GPU offers a great improvement in the running time of COMEBin. To use COMEBin in an enviroment which allows the usage of a GPU simply provide that enrioment's name to the default name of the enviroment for COMEBin and also modify the path to the COMEBin directory of the enviroment. These paths refer to the following options of ProteoSeeker:</p>

~~~bash
   -sen/--comebin-env             Str -Opt: ps_comebin- The conda enviroment for sra tools.
                                  'None/none': To not use an enviroment at all.

   -cfp/--comebin-folder-path     Str -Opt- The path to the bin folder of COMEBin.
~~~

<p align="justify">For example, the second option in our system has the following value: "/home/compteam/anaconda3/envs/ps_comebin_gpu/bin/COMEBin"</p>

### Removing installation enviroments, files and diretories
<p align="justify">To remove the enviroments, all their files and the directories which enviroments and directories were created during the installation of ProteoSeeker (by running the "./install.sh" script) you can run the script below, in the installation directory from the installation directory. You can then delete the main directory of ProteoSeeker and all enviroments and files associated with installing ProteoSeeker will have been removed by your system.</p>

~~~bash
./remove.sh
~~~

## 2.2 Docker
<p align="justify">To install ProteoSeeker from DockerHub as a docker image docker must be installed in your system. To install docker in Linux follow the instructions provided by the link below:</p>

Docker engine for Linux: https://docs.docker.com/engine/install/ubuntu/

<p align="justify">Then download the image of ProteoSeeker from dockerhub. There are two versions. The "main" version contains the "Kraken 2/Bracken Refseq indexe Collection Standard-8 database" while the "light" version does not. Hence, the "main" version can be used directly to run the "seek" or the "taxonomy" mode of ProteoSeeker, specificaly through the "Kraken2" route. The "light" version can be used directly to run only the "seek" mode of ProteoSeeker. Neither version contains a protein database. The process of using a protein database through Docker is explained below. Both versions can be modified to utilize a protein database and thus be used to run the "taxonomy" mode of ProteoSeeker through the "COMEBin/MetaBinner" route.</p>

The "main" version has a size of  **29.6 GB**. To install the "main" version use one of the following commands:
~~~bash
sudo docker image pull proteoseeker
sudo docker image pull proteoseeker:latest
sudo docker image pull proteoseeker:main
~~~

The "light" version has a size of  **15.6 GB**. To install the "light" version use one of the following commands:
~~~bash
sudo docker image pull proteoseeker:light
~~~

### Phobius
For either case of installation process followed in order to use the topology and signal peptide predictions provided by Phobius you must download Phobius from https://phobius.sbc.su.se/data.html. Otherwise, ProteoSeeker will run without performing topology and signal peptide predictions in its seek functionality.

# 3. Use
## 3.1 Parameter file
<p align="justify">In general the easiest way to run ProteoSeeker is to use a parameter file. A parameter file should at the very least contain the parameters for the options of ProteoSeeker which are to be modified from their default values. Paramteter files for differect case-scenarios may be downloaded from *here*. These files are also generated when running the "install.sh" in the installation folder and more specifically one can generate them at any time when running the "parameter_files.sh" in the installation folder. The latter files contains all options and have parameters specificilly for those that are to be modified. A parameter file is the proposed way of using ProteoSeeker eitehr when installing from source or when installing it as a Docker image. In addition, we advise that the paths used as input to ProteoSeeker (for files or databases) contain no whitespaces, although ProteoSeeker is designed to handle them properly. </p>

## 3.2 Options
<p align="justify">The options of ProteoSeeker as a command-line tool, their default parameters and descriptions are described below.</p>

~~~bash
ProteoSeeker Version 1.0.0

python enzann.py [options] -i <Input_Folder_Path>

Option description:
1. Parameter type
2. Req: Required, Opt: Optional
3. Default value (shown if not empty or none)
4. Description
Terminology:
spd: seek profile database
tpd: taxonomy profile database
sfpd: seek filtered protein database
tfpd: taxonomy filtered protein database

Options:
---------Input and output options---------
   -i/--input                     Str -Req- The path of the folder with the input files. The
                                  input files can either be single-end or paired-end FASTQ
                                  files but not both.

   -sc/--sra-code                 Str -Opt- A RUN accession code from the SRA database of
                                  NCBI.

   -c/--contigs                   True/False -Opt: False- Indicates whether the files in the
                                  input folder are (non-compressed) files with contigs or
                                  genomes in FASTA format.

   -pi/--protein-input            True/False -Opt: False- Indicates whether the input folder
                                  contains a file with protein sequences in FASTA format.

   -a/--adapters                  Str -Req: adapters.fa- Path to the file with the adapters.

   -pdp/--protein-database-path   Str -Opt- Path to the protein database file.

   -kdp/--kraken-database-path    Str -Opt- Path to the kraken database folder.

   -psp/--profiles-seek-path      Str -Opt- Path to the seek profile database with the
                                  profiles associated with one or more protein families.

   -pyp/--profiles-phylo-path     Str -Opt- Path to the phylo profile database with the
                                  profiles associated with one or more protein families.

   -pbp/--profiles-broad-path     Str -Opt- Path to the profile database with the profiles to
                                  be searched in the each protein identified with at least one
                                  profile from the seek profile database.

   -sp/--swissprot-path           Str -Opt- Path to the Swiss-Prot protein database.

   -mp/--motifs-path              Str -Opt- Path to the file with the motifs.

   -pfp/--parameters-file-path    Str -Opt- The path to the file with the parameters and their
                                  values.

   -o/--output                    Str -Opt- Path to the output folder.

---------Protein family options---------
   -fc/--family-code              Str -Opt- The seek protein family codes.

   -fct/--family-code-taxonomy    Str -Opt- The phylo protein family codes.

   -dn/--database-name            Str -Opt- The seek profile database (spd) and seek filtered
                                  protein database name (sfdp).

   -dnt/--database-name-taxonomy  Str -Opt- The taxonomy profile database (tpd) and taxonomy
                                  filtered protein database name (tfdp).

   -sns/--seek-names-status       True/False -Opt: False- Determines whether the protein names
                                  used to filter the protein database and create the sfpd will
                                  be determined solely based on the protein names provided by
                                  the user (True) or solely based on protein names
                                  automatically identified with or without the addition
                                  protein names provided by the user (False).

   -spn/--seek-protein-names      Str -Opt- Protein names, divided by commas, given as input
                                  from the user and used to filter the seek protein database.
                                  If such protein names are indeed provided and -sns is False,
                                  then these protein names are added to the automatically
                                  identified ones. If such protein names are not provided and
                                  -sns is False, then only the aumatically identified ones are
                                  used for the filtering.

   -tns/--taxonomy-names-status   True/False -Opt: False- Determines whether the protein names
                                  used to filter the protein database and create the tfpd will
                                  be determined solely based on the protein names provided by
                                  the user (True) or solely based on protein names
                                  automatically identified with or without the addition
                                  protein names provided by the user (False).

   -tpn/--taxonomy-protein-names  Str -Opt- Protein names, divided by commas, given as input
                                  from the user and used to filter the taxonomy protein
                                  database. If such protein names are indeed provided and -tns
                                  is False, then these protein names are added to the
                                  automatically identified ones. If such protein names are not
                                  provided and -tns is False, then only the aumatically
                                  identified ones are used for the filtering.

   -nt/--name-threshold           Float -Opt: 0.5- The threshold used to filter the protein
                                  names associated with each protein family. Any protein name
                                  with a frequency below this threshold is omitted.

---------General options: Pipeline---------
   -sr/--seek-route               Int -Opt: 3- There are three analysis modes. The analysis
                                  mode determines the type of analysis by the seek
                                  functionality. '1': The seek functionality will only search
                                  for proteins to be annotated which include at least one of
                                  the profiles of the spd. '2': The seek functionality will
                                  only search for proteins to be annotated that have at least
                                  one hit against the sfpd with a low enough e-value. '3': The
                                  seek functionality includes both types of analysis '1' and
                                  '2'.

   -p/--paired-end                True/False -Opt: True- Indicates whether the files in the
                                  input folder are paired-end (True) or single-end (False)
                                  files.

   -k/--compressed                True/False -Opt: True- Indicates whether the files in the
                                  input folder are compressed (True) or not (False).

   -fpd/--filter-protein-database True/False -Opt: False- Determines whether the protein
                                  database will be filtered based on protein names to create
                                  the sfpd and tfpd.

   -ps/--preftech-size            Int -Opt: 20- The maximum file size to download in KB (K for
                                  kilobytes, M for megabytes, G gigabytes).

   -as/--adapters-status          Str -Opt: 'pre'- The following options are available: ide:
                                  Adds the overrepresented sequences identified by FastQC in
                                  the file with the adapters. 'fas': The file with the
                                  adapters will include only the overrepresented sequences
                                  identified by FastQC. 'pre': The  file with the adapters is
                                  used without any modification.

   -asi/--add-seek-info           True/False -Opt: True- Determines whether the results in the
                                  TXT and the EXCEL file will only contain information for the
                                  proteins identified through the seek mode (True) or not
                                  (False). In case only the taxonomy mode is applied, this
                                  option has no effect on the results.

   -ati/--add-taxonomy-info       True/False -Opt: True- Determines whether the results in the
                                  TXT and the EXCEL file will only contain information for the
                                  proteins characterized through the taxonomy mode (True) or
                                  not (False). The latter proteins are the ones encoded by
                                  genes which are part of contigs that are grouped in bins,
                                  which bins have also been associated with at least one
                                  species. In case only the seek mode is applied, this option
                                  has no effect on the results.

   -h/--help                      None -- Displays the help message.

---------General options: FastQC---------
   -sf/--skip-fastqc              True/False -Opt: False- Determines whether the second time
                                  the FastQC analysis is applied, at the preprocessed reads,
                                  will be omitted (True) or not (False).

---------General options: BBDuk---------
   -umr/--bbduk-max-ram           Int -Opt 4: False- The maximum number of GBs of RAM that may
                                  be utilized by BBDuk.

   -cs/--clear-space              True/False -Opt: False- Determines whether the compressed
                                  preprocessed reads (if any) will be deleted (True) before
                                  the assembly or not (False).

---------General options: Megahit---------
   -kl/--k-list                   Str -Opt- A list of k-mers (e.g.,
                                  15,17,21,29,39,59,79,99,119,141) to be used by Megahit.

---------General options: Kraken---------
   -km/--kraken-mode              True/False -Opt: True- Determines whether the taxonomy
                                  functionality will be based on the taxonomy analysis applied
                                  by kraken2 (True) or not (False).

   -kt/--kraken-threshold         Int/Float -Opt -1: False- A list with read-filtering
                                  threshold for the species reported by kraken. The list
                                  should include integers of floats seperated by commads. An
                                  integer is used as an absolute read threshold and a float is
                                  used as a percentage threshold applied to the percentagies
                                  reported by kraken for each species (e.g., 100 to represent
                                  a threshold of 100 reads, 1 to represent a threshold of 1
                                  read, 1.0 to represent a threshold of 1%, 12.5 to represent
                                  a threshold of 12.5%). In addition, the values of -1 or -2
                                  can be provided, to automatically set the threshold. For the
                                  value of -1 the threshold is set specifically for non-gut
                                  metagenomes and for the value of -2 the threshold is set
                                  specifically for gut metagenomes. When kraken is selected a
                                  binning process takes place based on the filtered species
                                  from the results of kraken. The latter binning process is
                                  based on the filtering performed based on the first
                                  threshold value of the list (if not only one).

   -kmm/--kraken-memory-mapping   True/False -Opt: True- Determines whether kraken2 will use
                                  memory mapping (True) or not (False). With memory mapping
                                  the analysis performed by kraken2 is slower but is not
                                  limited by the size of the RAM available at the time of the
                                  analysis, rather than by the free memory space of the disk.
                                  Without memory mapping the analysis performed by kraken2 is
                                  faster but is limited by the size of the RAM available at
                                  the time of the analysis.

---------General options: Binning---------
   -bt/--binning-tool             Int -Opt: 1- Determines the binning tool to be used by the
                                  functionality of taxonomy, when kraken2 is set not to be
                                  used (-km False). '1': MetaBinner. '2': COMEBin.

   -bmr/--binning-max-ram         Int -Opt: 4- The maximum number of GBs of RAM that may be
                                  utilized by binning.

   -bc/--bin-contig-len           Int -Opt: 500- The threshold for filtering the contigs based
                                  on their lengths before binning. Any contig with length
                                  below the threshold is omitted from the binning process.

   -bk/--bin-kmer                 Int -Opt: 4- The number of kmers to be used by the binning
                                  tool (MetaBinner or COMEBin).

   -cbs/--comebin-batch-size      Int -Opt: 256- The batch size to be used for the analysis of
                                  COMEBin.

---------General options: CD-HIT---------
   -ct/--cdhit-threshold          Float -Opt: 0.99- The threshold used by CD-HIT. The value
                                  must be a float number between 0 and 1.

   -cmr/--cd-hit-max-ram          Int -Opt: 4000- The maximum number of MBs of RAM that may be
                                  utilized by CD-HIT.

---------General options: Gene prediction---------
   -ge/--gene-encoding            Int -Opt: 1- Determines whether the proteins will be
                                  provided directly from the output of FragGeneScanRs (1) or
                                  will be the output (2) from applying the genetic code
                                  indicated by option -gc to encode the predicted genes from
                                  FragGeneScanRs.

   -gc/--genetic-code             Int -Opt: 11- The genetic code to be used to encode the
                                  genes predicted by FragGeneScanRs to proteins, if such an
                                  action has been selected (-ge 2).

---------General options: HMMER---------
   -st/--score-type               Str -Opt cut_ga- The scoring method used by HMMER. 'cut_ga':
                                  HMMER will use the GA gathering cutoffs of the profile to set all
                                  thresholding. 'default': HMMER will use its default scoring
                                  method.

   -sds/--second-domain-search    True/False -Opt: True- Determines whether the screening of
                                  the proteins by hmmscan of HMMER against the Pfam database
                                  will be applied (True) or not (False) during the seek
                                  functionality.

   -ndt/--no-domains-thr          Int -Opt: 70- The negative of this number becomes the power
                                  of 10 and the result is the e-value threshold used to retain
                                  proteins during analysis mode 2 in the seek functionality
                                  for further annotation.

---------General options: Annotation---------
   -at/--add-type                 Str -Opt- A comma-seperated list which includes kinds of
                                  information related to the analysis. The items of the list
                                  are added for each protein in each of the annotation files.

   -ai/--add-info                 Str -Opt- A comma-seperated list which includes values for
                                  the corresponding kinds of information set by -at. The items
                                  of the list are added for each protein in each of the
                                  annotation files.

---------Threads---------
   -t/--threads                   Int -Opt: 4- The maximum number of threads to be used by any
                                  of the processes used by ProteoSeeker.

   -ft/--filtering-threads        Int -Opt: -t- The maximum number of threads to be used by
                                  the filtering process of the protein database by
                                  ProteoSeeker. If not modified, the default value is equal to
                                  the value given at -t. Otherwise, it overwrites the value of
                                  -t.

---------Processes performed after---------
   -afp/--after-peprocessing      True/False -Opt: False- The pipeline start after the
                                  preprocessing of the reads.

   -afa/--after-assembly          True/False -Opt: False- The pipeline starts after the
                                  assembly of the reads. The files generated by the previous
                                  steps should be present in the output folder provided by the
                                  user.

   -afg/--after-gene-pred         True/False -Opt: False- The pipeline starts after gene
                                  prediction. The files generated by the previous steps should
                                  be present in the output folder provided by the user.

   -afb/--after-binning           True/False -Opt: False- The pipeline starts after binning.
                                  The files generated by the previous steps should be present
                                  in the output folder provided by the user.

   -adb/--after-after-db          True/False -Opt: False- The pipeline starts after screening
                                  against the filtered protein database in both cases of the
                                  seek and taxonomy functionalities. The files generated by
                                  the previous steps should be present in the output folder
                                  provided by the user.

   -atp/--after-topology-prediction True/False -Opt: False- The pipeline starts after the
                                  topology prediction in the case of the seek functionality.
                                  The files generated by the previous steps should be present
                                  in the output folder provided by the user.

   -afr/--after-analysis-processes True/False -Opt: False- The pipeline starts after all the
                                  analysis processes and only writes the annotation files for
                                  the proteins, both cases of the seek and taxonomy
                                  functionalities. The files generated by the previous steps
                                  should be present in the output folder provided by the user.

---------Processes performed up to---------
   -uts/--up-to-sra               True/False -Opt: False- The pipeline ends after downloading
                                  and processing the sample corresponding to the SRA code.

   -utd/--up-to-databases         True/False -Opt: False- The pipeline ends after creating the
                                  seek and taxonomy profile databases and the seek and
                                  taxonomy filtered protein databases are created.

   -utpc/--up-to-preprocessing-com True/False -Opt: False- The pipeline ends after the
                                  preprocessing of the reads.

   -utpu/--up-to-preprocessing-uncom True/False -Opt: False- The pipeline ends after
                                  decompressing the compressed preprocessed reads.

   -uta/--up-to-assembly          True/False -Opt: False- The pipeline ends after the assembly
                                  of the reads.

---------Tool enviroments---------
   -sen/--sra-env                 Str -Opt: ps_sra_tools- The conda enviroment for sra tools.
                                  'None/none': To not use an enviroment at all.

   -fen/--fastqc-env              Str -Opt: ps_fastqc- The conda enviroment for FastQC.
                                  'None/none': To not use an enviroment at all.

   -uen/--bbtools-env             Str -Opt: ps_bbtools- The conda enviroment for bbtools.
                                  'None/none': To not use an enviroment at all.

   -men/--megahit-env             Str -Opt: ps_megahit- The conda enviroment for megahit.
                                  'None/none': To not use an enviroment at all.

   -ken/--kraken-env              Str -Opt: ps_kraken- The conda enviroment for kraken2.
                                  'None/none': To not use an enviroment at all.

   -nen/--metabinner-env          Str -Opt: ps_metabinner- The conda enviroment for
                                  MetaBinner. 'None/none': To not use an enviroment at all.

   -sen/--comebin-env             Str -Opt: ps_comebin- The conda enviroment for sra tools.
                                  'None/none': To not use an enviroment at all.

   -ien/--cdhit-env               Str -Opt: ps_cd_hit- The conda enviroment for CD-HIT.
                                  'None/none': To not use an enviroment at all.

   -gen/--genepred-env            Str -Opt- The conda enviroment for FragGeneScanRs.
                                  'None/none': To not use an enviroment at all.

   -hen/--hmmer-env               Str -Opt: ps_hmmer- The conda enviroment for HMMER.
                                  'None/none': To not use an enviroment at all.

   -den/--dimaond-env             Str -Opt: ps_diamond- The conda enviroment for DIMAOND
                                  BLASTP. 'None/none': To not use an enviroment at all.

   -ten/--taxonkit-env            Str -Opt: ps_taxonkit- The conda enviroment for taxonkit.
                                  'None/none': To not use an enviroment at all.

   -pen/--phobius-env             Str -Opt: ps_phobius- The conda enviroment for Phobius.
                                  'None/none': To not use an enviroment at all.

   -ben/--bowtie-env              Str -Opt: ps_bowtie- The conda enviroment for Bowtie2.
                                  'None/none': To not use an enviroment at all.

---------Tool paths---------
   -adp/--anaconda-dir-path       Str -Opt- The path to the anaconda installation directory.
                                  This directory includes directories like "bin" and "etc".

   -asp/--anaconda-sh-path        Str -Opt- The path to conda.sh. If provided the path to
                                  conda.sh will not be automatically determined by the path to
                                  the conda installation directory (-adp).

   -rfp/--prefetch-path           Str -Opt- The path to the prefetch executable.

   -vvp/--vdb_validate-path       Str -Opt- The path to the vdb-validate executable.

   -fdp/--fastq-dump-path         Str -Opt- The path to the fastq-dump executable.

   -fp/--fastqc-path              Str -Opt- The path to the fastqc executable.

   -gzp/--gzip-path               Str -Opt- The path to the gzip executable.

   -ctp/--cat-path                Str -Opt- The path to the cat executable.

   -bdp/--bbduk-path              Str -Opt- The path to the bbduk executable.

   -mp/--megahit-path             Str -Opt- The path to the megahit executable.

   -kp/--kraken-path              Str -Opt- The path to the kraken executable.

   -bfp/--binner-folder-path      Str -Opt- The path to the bin folder of MetaBiner.

   -cfp/--comebin-folder-path     Str -Opt- The path to the bin folder of COMEBin.

   -chp/--cd-hit-path             Str -Opt- The path to the CD-HIT executable.

   -fgp/--fraggenescars-path      Str -Opt- The path to the FragGeneScanRs executable.

   -hp/--hmmscan-path             Str -Opt- The path to the hmmscan executable.

   -hpp/--hmmpress-path           Str -Opt- The path to the hmmpress executable.

   -hfp/--hmmfetch-path           Str -Opt- The path to the hmmfetch executable.

   -dp/--diamond-path             Str -Opt- The path to the diamond executable.

   -tkp/--taxonkit-path           Str -Opt- The path to the taxonkit executable.

   -php/--phobius-folder-path     Str -Opt- The path to the folder of phobius.

   -bbp/--bowtie-build-path       Str -Opt- The path to the bowtie build executable.
~~~

## 3.3 Command-line
<p align="justify">To run ProteoSeeker through the command-line a parameter file facilitates the process greatly. By using one of the template parameter files one can easily customize the values for the options of ProteoSeeker and run it. To run ProteoSeeker, at first one should activate its enviroment. The commands below suffice to run ProteoSeeker through the command-line from the main directory of ProteoSeeker by its "taxonomy" mode through the Kraken2 route. To run ProteoSeeker by its seek mode one should set the path to the protein database in the parameters file and then select one of the template files that include the "seek" mode at its name.</p>

~~~bash
conda activate ps_env
python proteoseeker.py -pfp par_tax_k_p.txt
conda deactivate
~~~

## 3.4 Docker
<p align="justify">We present two ways to run ProteoSeeker through its image. Both ways depend on creating a volume or bind-mount and attaching it to the container running based on the image. The first way is running ProteoSeeker directly with creating a container. The second way is to start a container in interactive mode and then run ProteoSeeker. In addition, in either case the volume or bind-mount can be used to provide a tool or database to ProteoSeeker in the container. In both cases (volume and bind mount) the protein database provided as an example is a small part of the nr database with additions of proteins associated with RNA polymerase, as for the file to work to test the functionality of the "seek" mode through "type 2" analysis and that of the "taxonomy" mode of ProteoSeeker through the route of "COMEBin/MetaBinner". You should provide your own protein database, ideally that of the decompressed nr database, in order to use properly the "seek" mode through the "type 2" analysis and the "taxonomy" mode through the "COMEBin/MetaBinner" route of ProteoSeeker.</p>

### 3.4.1 Volume
<p align="justify">A volume is a directory inside Docker. Volumes can be found in the "volumes" directory of your Docker installation (e.g., /var/lib/docker/volumes). The data of the volume is stored in the "_data" directory of the volume. This data are retained in the volume after the container is stopped or exits, may be used by different containers and are also accessible by the local host. Any directory or file placed in the "_data" directory will be accessible from the local host and the container to which is has been added. From the **main direcotry** of ProteoSeeker run the Bash script with sudo (mandatory) below which performs the following actions:</p>

1. Creates a Docker volume.
2. Finds the full path of the Docker volume.
3. Creates directories in the volumes needed to run ProteoSeeker in the image and collect the results in the host system.
4. Copies the parameter files used as examples to run ProteoSeeker in the image.

~~~bash
sudo ./installation/docker_vol_setup.sh
~~~

<p align="justify">Then run ProteoSeeker in the Docker image based on one of the parameter files. The following script will run ProteoSeeker in the Docker image in interactive mode so the user can observe the stages of the pipeline being run. As an example you can select to run either a "seek" and "taxonomy" analysis with the Kraken2 route (selection: 1) or only a "taxonomy" analysis with the Kraken2 route (selection: 2).</p>

~~~bash
./docker_vol_run_proteoseeker.sh
~~~

<p align="justify">So in this case, the data shared between the local host and the container will be located at "/var/lib/docker/volumes/ps_vol/_data". To find the results of the run check in the directory "/var/lib/docker/volumes/ps_vol/_data/results". Similarly, you can run ProteoSeeker based on a custom parameters file. After ProteoSeeker has terminated the container stops and exits, hence one can not find the results in that container. One can also run the same command of the script by adding the "-d" option which runs the container in the background. He can then attach to the container to observe which stage the pipeline is currently runnning. In addition, one can use "/bin/bash/ as the command to run in the container and enter the container in interactive mode, then run ProteoSeeker from the container (as one would directly from the command-line of the host) and then transfer the results to the host through a shared volume or bind.</p>

### 3.4.2 Bind mount
<p align="justify">A volume is a directory located in the local host and not run by Docker. As for the volume, the data stored in the mount are reatined after the container is stopped or exits, may be used by different containers and are also accessible by the local host. From the **main direcotry** of ProteoSeeker run the Bash script below (should not be run with sudo, otherwise it might not work) which performs the following actions:</p>

1. Creates the "docker_mount_dir" directory in the home directory of the user.
2. Create directories in the volumes needed to run ProteoSeeker in the image and collect the results in the host system.
3. Copies the parameter files used as examples to run ProteoSeeker in the image.

~~~bash
./installation/docker_bindmount_setup.sh
~~~

<p align="justify">Then run ProteoSeeker in the Docker image based on one of the parameter files. The following script will run ProteoSeeker in the Docker image in interactive mode so the user can observe the stages of the pipeline being run. As an example you can select to run either a "seek" and "taxonomy" analysis with the Kraken2 route (selection: 1) or only a "taxonomy" analysis with the Kraken2 route (selection: 2).</p>

~~~bash
./docker_bind_run_proteoseeker.sh
~~~

<p align="justify">When running the Docker container and collecting the results in a directory located in the "${HOME}" path of the user, one should first change the permissions of the results directory by providing read and write permissions to all its contents, otherwise the directory generated as output from ProteoSeeker from the Docker container will be missing these permissions.</p>

~~~bash
sudo chmod -R a+rw "${HOME}/proteoseeker_bindmount/results"
~~~

<p align="justify">So in this case, the data shared between the local host and the container will be located at "/${home}/docker_mount_dir". To find the results of the run check in the directory "/${home}/docker_mount_dir/results". One can also run the same command of the script by adding the "-d" option which runs the container in the background. He can then attach to the container to observe which stage the pipeline is currently runnning.</p>

<p align="justify">In addition, in both cases of volumes and bind mounts, one can use "/bin/bash/ as the command to run in the container and enter the container in interactive mode, then run ProteoSeeker from the container (as one would directly from the command-line of the host) and then transfer the results to the host through the shared volume or bind mount.</p>

### 3.4.3 Utilizing Phobius in the docker container
<p align="justify">The "phobius" directory in the volume or bind mount can be used to store the files of the Phobius installation. If done so, ProteoSeeker will detect and use Phobius when running in the Docker container (except if set otherwise in the parameters file). Download and install Phobius based on the instructions from: https://phobius.sbc.su.se/data.html</p>

# 4. Test cases
<p align="justify">All tests for the evaluation were run based on the ProteoSeeker version 1.0.0 and the tool versions described in it, which refer to the current release of ProteoSeeker. The collection dates for the databases used in the evaluation are found below. In addition, we also mention the download date of the flat file for the reviewed (Swiss-Prot) proteins of the Uniprot database, which was used to collect information about protein families, protein Pfam domains, protein names, protein lengths and their associations.</p>

1. Pfam database: 29/05/2024
2. Swiss-Prot/UniprotKB database: 29/05/2024
3. GTDB taxonomy taxdump files: 29/05/2024
4. Kraken 2/Bracken Refseq indexes: Collection Standard-8: 05/06/2024 (prior to the update)
5. Kraken 2/Bracken Refseq indexes: Collection Standard-16: 05/06/2024 (prior to the update)
6. Kraken 2/Bracken Refseq indexes: Collection Standard: 05/06/2024 (prior to the update)
7. nr database: 27/06/2024
8. Reviewed (Swiss-Prot) flat file: 04/08/2023

## 4.1 Seek and taxonomy modes evaluation
<p align="justify">The seek and taxonomy mode of ProteoSeeker can be tested on specific samples which have been analyzed by our lab in the process of discovering novel enzymes with desirable characteristics. The following samples were used for these analyses: "+++", "+++" and "+++". More details about the datasets can be found in the [MANUSCRIPT - UNDER REVIEW]. To analyze the results from these analyses you can run the "analyze_seek_tax_results.sh" script from the "tests/seek_tax_results_analysis" directory which utilizes a Python script from the same directory. The latter script uses an input file which contains information about the species of the best hit of each of the evaluated proteins against the nr database through BLASTP. The best hit is identified based on the lowest e-vaalue amongst all the hits."</p>

## 4.2 Taxonomy mode evaluation
<p align="justify">To run ProteoSeeker on the 19 benchmark datasets which correspond to the gold standard populations and which were used to evaluate the taxonomy mode of ProteoSeeker, the folowing steps must be followed.</p>

1. Create the database for the taxonomy route of COMEBin/MetaBinner. Copy the parameter file "taxonomy_tests/dbs/par_phylo_dbs_nr_rna_pol.txt" to the ProteoSeeker directory. Modify the values for the parameters for the path of the protein database (nr database) and the output path. Then, activate "ps_env" and run:

~~~bash
python proteoseeker.py -pfp par_phylo_dbs_nr_rna_pol.txt
~~~

2. Modify the "par_sX_sra.txt" file for sample X based on the paths to the databases and tools installed in your system. For each sample to be downloaded and processed by ProteoSeeker separately and not during the analysis run, modify its corresponding parameter file and run from the parent directory of ProteoSeeker the script:

~~~bash
./taxonomy_tests/sra_process/run_sra_1.sh
~~~

3. In the sciprt "taxonomy_tests/setup.sh" modify the paths for PS_PATH and HDD_PATH. The first path is the location of the installation of ProteoSeeker. The second path is a directory where all of the results generated from ProteoSeeker for every run will be stored. Only the files with taxonomy related information are retained in the results directory in PS_PATH. This is helpfull because the size of all the results generated by each run for a sample might be very large and may be suitable to be stored elsewhere.

4. Run the Python script which creates the parameters files for ProteoSeeker and the Bash script files to be used for running the taxonomy evaluation.

~~~bash
python parameter_files/create_par_files.py
~~~

5. Once the Python script creates the parameter files and the Basch scripts there are four Bash scripts associated with running each sample for all methods, "run_X_all.sh", "run_X_kraken.sh", "run_X_metabinner.sh" and "run_X_comebin.sh". The first one runs all methods, the second one runs Kraken2 with every database, the third one runs MetaBinner with nr and the last one runs COMEBin with nr. The first (MetaBinner) and third (COMEBin) one can be run at any time. The MetaBinner and COMEBin methods require for the Kraken2 Basch script to have run beforehand. To run any of the Basch scripts you must be located at the installation directory of ProteoSeeker. For example, after activating ps_env:

~~~bash
./taxonomy_tests/run_1_kraken.sh
~~~

6. To extract the statistics from the results of ProteoSeeker for the taxonomy evaluation all the methods from at least one sample must have run. The Bash script will analyze and provide output (plots, statistics etc.) for any number of samples analyzed based on all methods. To run the script simply move to the directory "/tests/tax_results_analysis" and run the following:

~~~bash
./tests/tax_results_analysis/analyze_tax_results.sh
~~~

The "analyze_tax_results_spec.sh" can be used for secondry analysis of different subsets of the results. This script is located also in the directory "/tests/tax_results_analysis".
