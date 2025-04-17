def help_par(text, step):
    words = text.split(" ")
    lines = [words[0]]
    for word in words[1:]:
        if len(lines[-1]) + len(word) < step:
            lines[-1] += (" " + word)
        else:
            lines.append(word)
    return lines


def help_message():
    print("ProteoSeeker Version 1.0.0")
    print()
    print("Usage:")
    print("python proteoseeker.py -pfp <parameters_file_path> [other_options]")
    print("or")
    print("python proteoseeker.py [options]")
    print()
    print("Option description:")
    print("1. Parameter type")
    print("2. Req: Required, Opt: Optional")
    print("3. Default value (shown if not empty or none)")
    print("4. Description")
    print()
    print("Terminology:")
    print("SPD: seek profile database")
    print("TPD: taxonomy profile database")
    print("SFPD: seek filtered protein database")
    print("TFPD: taxonomy filtered protein database")
    print()
    print("Options:")
    print("---------Input and output options---------")
    help_notes_dict = {
        "-i/--input": "Str -Req- The path of the directory with the input files. The input files can either be single-end or paired-end FASTQ files but not both.",
        "-sc/--sra-code": "Str -Opt- A RUN accession code from the SRA database of NCBI.",
        "-c/--contigs": "True/False -Opt: False- Indicates whether the files in the input directory are (non-compressed) files with contigs or genomes in FASTA format.",
        "-pi/--protein-input": "True/False -Opt: False- Indicates whether the input directory contains a file with protein sequences in FASTA format.",
        "-a/--adapters": "Str -Req: adapters.fa- Path to the file with the adapters.",
        "-pdp/--protein-database-path": "Str -Opt- Path to the protein database file.",
        "-kdp/--kraken-database-path": "Str -Opt- Path to the kraken database directory.",
        "-psp/--profiles-seek-path": "Str -Opt- Path to the seek profile database with the profiles associated with one or more protein families.",
        "-ptp/--profiles-taxonomy-path": "Str -Opt- Path to the phylo profile database with the profiles associated with one or more protein families.",
        "-pbp/--profiles-broad-path": "Str -Opt- Path to the profile database with the profiles to be searched in the each protein identified with at least one profile from the seek profile database.",
        "-sp/--swissprot-path": "Str -Opt- Path to the Swiss-Prot protein database.",
        "-mp/--motifs-path": "Str -Opt- Path to the file with the motifs.",
        "-pfp/--parameters-file-path": "Str -Opt- The path to the file with the parameters and their values.",
        "-o/--output": "Str -Opt- Path to the output directory.",
        "-sm/--seek-mode": "True/False -Opt: True- Determines whether the seek mode will be applied (True) or not (False).",
        "-tm/--taxonomy-mode": "True/False -Opt: False- Determines whether the taxonomy mode will be applied (True) or not (False).",
        "-sr/--seek-route": "Int -Opt: 1- Determines the seek route to be applied in the seek mode. The seek route refers to the seek analysis mode. The first seek route (1) is based on screening the proteins based on the domains corresponding to the selected seek protein famlily codes. The second seek route (2) is based on screening the proteins against the seek filtered protein database (accoring to the seek protein names). The third seek route (3) is based on applying both the first and the second seek route.",
        "-tr/--taxonomy-route": "Int -Opt: 1- Determines the taxonomy route to be applied in the taxonomy mode. The taxonomy routes are the Kraken2 (1) and the COMEBin/MetaBinner (2) taxonomy routes.",
        "-sfc/--seek-family-code": "Str -Opt- The seek protein family codes.",
        "-tfc/--taxonomy-family-code": "Str -Opt- The taxonomy protein family codes.",
        "-dn/--database-name": "Str -Opt- The seek profile database (SPD) and seek filtered protein database name (SFPD).",
        "-dnt/--database-name-taxonomy": "Str -Opt- The taxonomy profile database (TPD) and taxonomy filtered protein database name (TFPD).",
        "-sns/--seek-names-status": "True/False -Opt: False- Determines whether the protein names used to filter the protein database and create the SFPD will be determined solely based on the protein names provided by the user (True) or solely based on protein names automatically identified with or without the addition protein names provided by the user (False).",
        "-spn/--seek-protein-names": "Str -Opt- Protein names, divided by commas, given as input from the user and used to filter the seek protein database. If such protein names are indeed provided and -sns is False, then these protein names are added to the automatically identified ones. If such protein names are not provided and -sns is False, then only the aumatically identified ones are used for the filtering.",
        "-tns/--taxonomy-names-status": "True/False -Opt: False- Determines whether the protein names used to filter the protein database and create the TFPD will be determined solely based on the protein names provided by the user (True) or solely based on protein names automatically identified with or without the addition protein names provided by the user (False).",
        "-tpn/--taxonomy-protein-names": "Str -Opt- Protein names, divided by commas, given as input from the user and used to filter the taxonomy protein database. If such protein names are indeed provided and -tns is False, then these protein names are added to the automatically identified ones. If such protein names are not provided and -tns is False, then only the aumatically identified ones are used for the filtering.",
        "-nt/--name-threshold": "Float -Opt: 0.5- The threshold used to filter the protein names associated with each protein family. Any protein name with a frequency below this threshold is omitted.",
        "-p/--paired-end": "True/False -Opt: True- Indicates whether the files in the input directory are paired-end (True) or single-end (False) files.",
        "-k/--compressed": "True/False -Opt: True- Indicates whether the files in the input directory are compressed (True) or not (False).",
        "-fpd/--filter-protein-database": "True/False -Opt: False- Determines whether the protein database will be filtered based on protein names to create the SFPD and TFPD.",
        "-ps/--preftech-size": "Int -Opt: 20- The maximum file size to download in KB (K for kilobytes, M for megabytes, G gigabytes).",
        "-as/--adapters-status": "Str -Opt: 'pre'- The following options are available: ide: Adds the overrepresented sequences identified by FastQC in the file with the adapters. 'fas': The file with the adapters will include only the overrepresented sequences identified by FastQC. 'pre': The  file with the adapters is used without any modification.",
        "-asi/--add-seek-info": "True/False -Opt: True- Determines whether the results in the TXT and the EXCEL file will only contain information for the proteins identified through the seek mode (True) or not (False). In case only the taxonomy mode is applied, this option has no effect on the results.",
        "-ati/--add-taxonomy-info": "True/False -Opt: True- Determines whether the results in the TXT and the EXCEL file will only contain information for the proteins characterized through the taxonomy mode (True) or not (False). The latter proteins are the ones encoded by genes which are part of contigs that are grouped in bins, which bins have also been associated with at least one species. In case only the seek mode is applied, this option has no effect on the results.",
        "-h/--help": "None -- Displays the help message.",
        "-sf/--skip-fastqc": "True/False -Opt: False- Determines whether the second time the FastQC analysis is applied, at the preprocessed reads, will be omitted (True) or not (False).",
        "-umr/--bbduk-max-ram": "Int -Opt 4: False- The maximum number of GBs of RAM that may be utilized by BBDuk.",
        "-cs/--clear-space": "True/False -Opt: False- Determines whether the compressed preprocessed reads (if any) will be deleted (True) before the assembly or not (False).",
        "-kl/--k-list": "Str -Opt- A list of k-mers (e.g., 15,17,21,29,39,59,79,99,119,141) to be used by Megahit.",
        "-kt/--kraken-threshold": "Float/Int -Opt: 0.1- A list with read-filtering threshold for the species reported by kraken. The list should include integers of floats seperated by commads. An integer is used as an absolute read threshold and a float is used as a percentage threshold applied to the percentagies reported by kraken for each species (e.g., 100 to represent a threshold of 100 reads, 1 to represent a threshold of 1 read, 1.0 to represent a threshold of 1%, 12.5 to represent a threshold of 12.5%). In addition, the values of -1 or -2 can be provided, to automatically set the threshold. For the value of -1 the threshold is set specifically for non-gut metagenomes and for the value of -2 the threshold is set specifically for gut metagenomes. When kraken is selected a binning process takes place based on the filtered species from the results of kraken. The latter binning process is based on the filtering performed based on the first threshold value of the list (if not only one).",
        "-kmm/--kraken-memory-mapping": "True/False -Opt: True- Determines whether kraken2 will use memory mapping (True) or not (False). With memory mapping the analysis performed by kraken2 is slower but is not limited by the size of the RAM available at the time of the analysis, rather than by the free memory space of the disk. Without memory mapping the analysis performed by kraken2 is faster but is limited by the size of the RAM available at the time of the analysis.",
        "-bl/--bracken-length": "Int -Opt: 100- The ideal length of reads in the sample. The read length based on which to get all classifications. A kmer distribution file should exist for the selected read length in the Kraken2 database directory provided for the Kraken2 analysis.",
        "-bl/--bracken-level": "Str -Opt: 'S'- Specifies the taxonomic rank to analyze. Each classification at this specified rank will receive an estimated number of reads belonging to that rank after abundance estimation. Available values for selection: 'D','P','C','O','F','G','S'",
        "-bh/--bracken-threshold": "Int -Opt: 10- Specifies the minimum number of reads required for a classification at the specified rank. Any classifications with less than the specified threshold will not receive additional reads from higher taxonomy levels when distributing reads for abundance estimation.",
        "-bt/--binning-tool": "Int -Opt: 1- Determines the binning tool to be used by the functionality of taxonomy, when kraken2 is set not to be used (-km False). '1': MetaBinner. '2': COMEBin.",
        "-bmr/--binning-max-ram": "Int -Opt: 4- The maximum number of GBs of RAM that may be utilized by binning.",
        "-bc/--bin-contig-len": "Int -Opt: 500- The threshold for filtering the contigs based on their lengths before binning. Any contig with length below the threshold is omitted from the binning process.",
        "-bk/--bin-kmer": "Int -Opt: 4- The number of kmers to be used by the binning tool (MetaBinner or COMEBin).",
        "-cbs/--comebin-batch-size": "Int -Opt: 256- The batch size to be used for the analysis of COMEBin.",
        "-ct/--cdhit-threshold": "Float -Opt: 0.99- The threshold used by CD-HIT. The value must be a float number between 0 and 1.",
        "-cmr/--cd-hit-max-ram": "Int -Opt: 4000- The maximum number of MBs of RAM that may be utilized by CD-HIT.",
        "-ge/--gene-encoding": "Int -Opt: 1- Determines whether the proteins will be provided directly from the output of FragGeneScanRs (1) or will be the output (2) from applying the genetic code indicated by option -gc to encode the predicted genes from FragGeneScanRs.",
        "-gc/--genetic-code": "Int -Opt: 11- The genetic code to be used to encode the genes predicted by FragGeneScanRs to proteins, if such an action has been selected (-ge 2).",
        "-st/--score-type": "Str -Opt cut_ga- The scoring method used by HMMER. 'cut_ga': HMMER will use the GA gathering cutoffs of the profile to set all thresholding. 'default': HMMER will use its default scoring method.",
        "-sds/--second-domain-search": "True/False -Opt: True- Determines whether the screening of the proteins by hmmscan of HMMER against the Pfam database will be applied (True) or not (False) during the seek functionality.",
        "-ndt/--no-domains-thr": "Int -Opt: 70- The negative of this number becomes the power of 10 and the result is the e-value threshold used to retain proteins during analysis mode 2 in the seek functionality for further annotation.",
        "-at/--add-type": "Str -Opt- A comma-seperated list which includes kinds of information related to the analysis. The items of the list are added for each protein in each of the annotation files.",
        "-ai/--add-info": "Str -Opt- A comma-seperated list which includes values for the corresponding kinds of information set by -at. The items of the list are added for each protein in each of the annotation files.",
        "-t/--threads": "Int -Opt: 4- The maximum number of threads to be used by any of the processes used by ProteoSeeker.",
        "-ft/--filtering-threads": "Int -Opt: -t- The maximum number of threads to be used by the filtering process of the protein database by ProteoSeeker. If not modified, the default value is equal to the value given at -t. Otherwise, it overwrites the value of -t.",
        "-fpar/--fastqc-parameters": "Str -Opt- A string which may contain any option and corresponding value for FastQC, except for the options related with input and output files.",
        "-bpar/--bbduk-parameters": "Str -Opt- A string which may contain any option and corresponding value for BBDuk, except for the options related with input and output files.",
        "-mpar/--megahit-parameters": "Str -Opt- A string which may contain any option and corresponding value for Megahit, except for the options related with input and output files.",
        "-rpar/--fraggenescanrs-parameters": "Str -Opt- A string which may contain any option and corresponding value for FragGeneScanRs, except for the options related with input and output files.",
        "-cpar/--cd-hit-parameters": "Str -Opt- A string which may contain any option and corresponding value for CD-HIT, except for the options related with input and output files.",
        "-hpar/--hmmerscan-parameters": "Str -Opt- A string which may contain any option and corresponding value for hmmerscan, except for the options related with input and output files.",
        "-dpar/--diamond-blastp-parameters": "Str -Opt- A string which may contain any option and corresponding value for DIAMOND BLASTP, except for the options related with input and output files.",
        "-kpar/--kraken-parameters": "Str -Opt- A string which may contain any option and corresponding value for Kraken2, except for the options related with input and output files.",
        "-npar/--bracken-parameters": "Str -Opt- A string which may contain any option and corresponding value for Bracken, except for the options related with input and output files.",
        "-epar/--metabinner-parameters": "Str -Opt- A string which may contain any option and corresponding value for MetaBinner, except for the options related with input and output files.",
        "-opar/--comebin-parameters": "Str -Opt- A string which may contain any option and corresponding value for COMEBin, except for the options related with input and output files.",
        "-wpar/--bowtie-parameters": "Str -Opt- A string which may contain any option and corresponding value for Bowtie2, except for the options related with input and output files.",
        "-afp/--after-peprocessing": "True/False -Opt: False- The pipeline starts after the preprocessing of the reads.",
        "-afa/--after-assembly": "True/False -Opt: False- The pipeline starts after the assembly of the reads. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-afg/--after-gene-pred": "True/False -Opt: False- The pipeline starts after gene prediction. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-afb/--after-binning": "True/False -Opt: False- The pipeline starts after binning. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-afm/--after-mapping": "True/False -Opt: False- The pipeline starts after mapping the reads to the contigs. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-adb/--after-after-db": "True/False -Opt: False- The pipeline starts after screening against the filtered protein database in both cases of the seek and taxonomy functionalities. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-atp/--after-topology-prediction": "True/False -Opt: False- The pipeline starts after the topology prediction in the case of the seek functionality. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-afr/--after-analysis-processes": "True/False -Opt: False- The pipeline starts after all the analysis processes and only writes the annotation files for the proteins, both cases of the seek and taxonomy functionalities. The files generated by the previous steps should be present in the output directory provided by the user.",
        "-uts/--up-to-sra": "True/False -Opt: False- The pipeline ends after downloading and processing the sample corresponding to the SRA code.",
        "-utd/--up-to-databases": "True/False -Opt: False- The pipeline ends after creating the seek and taxonomy profile databases and the seek and taxonomy filtered protein databases are created.",
        "-utpc/--up-to-preprocessing-com": "True/False -Opt: False- The pipeline ends after the preprocessing of the reads.",
        "-utpu/--up-to-preprocessing-uncom": "True/False -Opt: False- The pipeline ends after decompressing the compressed preprocessed reads.",
        "-uta/--up-to-assembly": "True/False -Opt: False- The pipeline ends after the assembly of the reads.",
        "-sen/--sra-env": "Str -Opt: ps_sra_tools- The conda environment for sra tools. 'None/none': To not use an environment at all.",
        "-fen/--fastqc-env": "Str -Opt: ps_fastqc- The conda environment for FastQC. 'None/none': To not use an environment at all.",
        "-uen/--bbtools-env": "Str -Opt: ps_bbtools- The conda environment for bbtools. 'None/none': To not use an environment at all.",
        "-men/--megahit-env": "Str -Opt: ps_megahit- The conda environment for megahit. 'None/none': To not use an environment at all.",
        "-ken/--kraken-env": "Str -Opt: ps_kraken- The conda environment for kraken2. 'None/none': To not use an environment at all.",
        "-ren/--bracken-env": "Str -Opt: ps_bracken- The conda environment for bracken. 'None/none': To not use an environment at all.",
        "-nen/--metabinner-env": "Str -Opt: ps_metabinner- The conda environment for MetaBinner. 'None/none': To not use an environment at all.",
        "-sen/--comebin-env": "Str -Opt: ps_comebin- The conda environment for COMEBin. 'None/none': To not use an environment at all.",
        "-ien/--cdhit-env": "Str -Opt: ps_cd_hit- The conda environment for CD-HIT. 'None/none': To not use an environment at all.",
        "-gen/--genepred-env": "Str -Opt- The conda environment for FragGeneScanRs. 'None/none': To not use an environment at all.",
        "-hen/--hmmer-env": "Str -Opt: ps_hmmer- The conda environment for HMMER. 'None/none': To not use an environment at all.",
        "-den/--dimaond-env": "Str -Opt: ps_diamond- The conda environment for DIMAOND BLASTP. 'None/none': To not use an environment at all.",
        "-ten/--taxonkit-env": "Str -Opt: ps_taxonkit- The conda environment for taxonkit. 'None/none': To not use an environment at all.",
        "-pen/--phobius-env": "Str -Opt: ps_phobius- The conda environment for Phobius. 'None/none': To not use an environment at all.",
        "-ben/--bowtie-env": "Str -Opt: ps_bowtie- The conda environment for Bowtie2. 'None/none': To not use an environment at all.",
        "-adp/--anaconda-dir-path": "Str -Opt- The path to the anaconda installation directory. This directory includes directories like \"bin\" and \"etc\".",
        "-asp/--anaconda-sh-path": "Str -Opt- The path to conda.sh. If provided the path to conda.sh will not be automatically determined by the path to the conda installation directory (-adp).",
        "-rfp/--prefetch-path": "Str -Opt- The path to the prefetch executable.",
        "-vvp/--vdb_validate-path": "Str -Opt- The path to the vdb-validate executable.",
        "-fdp/--fastq-dump-path": "Str -Opt- The path to the fastq-dump executable.",
        "-fp/--fastqc-path": "Str -Opt- The path to the fastqc executable.",
        "-gzp/--gzip-path": "Str -Opt- The path to the gzip executable.",
        "-ctp/--cat-path": "Str -Opt- The path to the cat executable.",
        "-bdp/--bbduk-path": "Str -Opt- The path to the bbduk executable.",
        "-mp/--megahit-path": "Str -Opt- The path to the megahit executable.",
        "-kp/--kraken-path": "Str -Opt- The path to the kraken executable.",
        "-bp/--bracken-path": "Str -Opt- The path to the bracken executable.",
        "-ap/--alpha-diversity-path": "Str -Opt- The path to the alpha diversity executable from KrakenTools.",
        "-bfp/--binner-folder-path": "Str -Opt- The path to the bin directory of MetaBiner.",
        "-cfp/--comebin-folder-path": "Str -Opt- The path to the parent directory of \"run_comebin.sh\" of COMEBin.",
        "-chp/--cd-hit-path": "Str -Opt- The path to the CD-HIT executable.",
        "-fgp/--fraggenescars-path": "Str -Opt- The path to the FragGeneScanRs executable.",
        "-hp/--hmmscan-path": "Str -Opt- The path to the hmmscan executable.",
        "-hp/--hmmscan-path": "Str -Opt- The path to the hmmscan executable.",
        "-hpp/--hmmpress-path": "Str -Opt- The path to the hmmpress executable.",
        "-hfp/--hmmfetch-path": "Str -Opt- The path to the hmmfetch executable.",
        "-dp/--diamond-path": "Str -Opt- The path to the diamond executable.",
        "-tkp/--taxonkit-path": "Str -Opt- The path to the taxonkit executable.",
        "-php/--phobius-folder-path": "Str -Opt- The path to the directory of phobius.",
        "-bbp/--bowtie-build-path": "Str -Opt- The path to the bowtie build executable.",
        "-hfp/--bowtie-path": "Str -Opt- The path to the bowtie executable."
    }
    split_step = 60
    for key_hn in help_notes_dict.keys():
        description = help_notes_dict[key_hn]
        des_length = len(description)
        if split_step >= des_length:
            print('{:2} {:30} {}\n'.format("", key_hn, description))
        else:
            pieces = help_par(description, split_step)
            for pmi in range(0, len(pieces)):
                piece_mod = pieces[pmi]
                if pmi == 0:
                    print('{:2} {:30} {}'.format("", key_hn, piece_mod))
                elif pmi + 1 == len(pieces):
                    print('{:33} {}\n'.format("", piece_mod))
                else:
                    print('{:33} {}'.format("", piece_mod))
        if key_hn  == "-o/--output":
            print("---------Modes and Routes---------")
        elif key_hn == "-tr/--taxonomy-route":
            print("---------Pipeline options---------")
        elif key_hn == "-h/--help":
            print("---------FastQC options---------")
        elif key_hn == "-sf/--skip-fastqc":
            print("---------BBDuk options---------")
        elif key_hn == "-cs/--clear-space":
            print("---------Megahit options---------")
        elif key_hn == "-kl/--k-list":
            print("---------Kraken options---------")
        elif key_hn == "-bh/--bracken-threshold":
            print("---------Binning options---------")
        elif key_hn == "-cbs/--comebin-batch-size":
            print("---------CD-HIT options---------")
        elif key_hn == "-cmr/--cd-hit-max-ram":
            print("---------Gene prediction options---------")
        elif key_hn == "-gc/--genetic-code":
            print("---------HMMER options---------")
        elif key_hn == "-ndt/--no-domains-thr":
            print("---------Annotation options---------")
        elif key_hn == "-ai/--add-info":
            print("---------Threads---------")
        elif key_hn == "-ft/--filtering-threads":
            print("---------Tool parameters---------")
        elif key_hn == "-wpar/--bowtie-parameters":
            print("---------Processes performed after---------")
        elif key_hn == "-afr/--after-analysis-processes":
            print("---------Processes performed up to---------")
        elif key_hn == "-uta/--up-to-assembly":
            print("---------Tool environments---------")
        elif key_hn == "-ben/--bowtie-env":
            print("---------Tool paths---------")
