import os
import shutil
# ProteoSeeker modules
import command_process
import supportive_functions


def translation(sequence, genetic_code):
    translation_table_1 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_2 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "*", "AGG": "*", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_3 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W", "CTT": "T", "CTC": "T",
        "CTA": "T", "CTG": "T", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_4 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_5 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "S", "AGG": "S", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_11 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_table_12 = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y",
        "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
        "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
        "GGG": "G"
    }
    translation_tables = {
        1: translation_table_1,
        2: translation_table_2,
        3: translation_table_3,
        4: translation_table_4,
        5: translation_table_5,
        11: translation_table_11,
        12: translation_table_12
    }

    # Form the protein sequence.
    protein = ""
    for i in range(0, len(sequence), 3):
        # If the whole sequence is less than 3 nucleotides, no tranlation occurs.
        if len(sequence) - 3 < i:
            continue
        else:
            # If the last nucleotides do not add up to three, then no tranlation occurs for the last 1 or 2 nucleotides.
            if i+1 >= len(sequence):
                continue
            translation_table = translation_tables[genetic_code]
            tripeptide = sequence[i] + sequence[i+1] + sequence[i+2]
            # Capitalize all nucleotide letters
            tripeptide = tripeptide.upper()
            # If the tripeptide contains at least one N then directly transalte to the "X" amino acid
            if "N" in tripeptide:
                amino_acid = "X"
            else:
                amino_acid = translation_table[tripeptide]
            protein += amino_acid
    return protein


def gene_prediction(genepred_env, output_path_genepred, output_path_genepred_base, output_path_genepred_faa, output_fgs_protein_formatted_1, fullpath_contigs_formated, fraggenescanrs_path, thread_num, input_log_file, output_log_file, gene_bash_script, fraggenescanrs_version_path, fraggenescanrs_stdoe_path, conda_sh_path, fraggenescanrs_par):
    print("\nPerofrming gene prediction...")
    if os.path.exists(output_path_genepred):
        shutil.rmtree(output_path_genepred)
    os.mkdir(output_path_genepred)
    
    # In the output of FragGeneScanRs in the headers, the first two last two numbers, which are divided by an underscore, consist the starting and ending points of the gene on the contig which it was found to be located in.
    # E.g. Header: >k141_180_2_133_- : This means that the name of the contig is "k141_180", the starting nucleotide on the contig is 2 and ending nucleotide on the contig is 133.
    # FragGeneScanRs will either be run to the contigs or the reads after they have been filtered by the trimming tool.
    input_log_file.write("Gene prediction:\n")
    capture_status = True
    shell_status = True

    # Tool parameters
    if fraggenescanrs_par is None:
        fraggenescanrs_par = " -w 1 -p {} --training-file complete ".format(thread_num)
    else:
        if fraggenescanrs_par[0] != " ":
            fraggenescanrs_par = " {}".format(fraggenescanrs_par)
        if fraggenescanrs_par[-1] != " ":
            fraggenescanrs_par = "{} ".format(fraggenescanrs_par)

    # Run FragGeneScanRs
    if fraggenescanrs_path:
        phrase_1 = "\"{}\" -o \"{}\" -s \"{}\"{}&> \"{}\"".format(fraggenescanrs_path, output_path_genepred_base, fullpath_contigs_formated, fraggenescanrs_par, fraggenescanrs_stdoe_path)
        phrase_2 = "\"{}\" --version > \"{}\"".format(fraggenescanrs_path, fraggenescanrs_version_path)
    else:
        phrase_1 = "FragGeneScanRs -o \"{}\" -s \"{}\"{}&> \"{}\"".format(output_path_genepred_base, fullpath_contigs_formated, fraggenescanrs_par, fraggenescanrs_stdoe_path)
        phrase_2 = "FragGeneScanRs --version > \"{}\"".format(fraggenescanrs_version_path)
    
    # Create the Bash script.
    # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
    new_file_bash = open(gene_bash_script, "w")
    phrase = "#!/bin/bash"
    new_file_bash.write("{}\n".format(phrase))
    if genepred_env:
        phrase_s = "source \"{}\"".format(conda_sh_path)
        phrase_a = "conda activate \"{}\"".format(genepred_env)
        new_file_bash.write("{}\n".format(phrase_s))
        new_file_bash.write("{}\n".format(phrase_a))
    new_file_bash.write("{}\n".format(phrase_1))
    new_file_bash.write("{}\n".format(phrase_2))
    if genepred_env:
        phrase = "conda deactivate"
        new_file_bash.write("{}\n".format(phrase))
    new_file_bash.close()

    # Making the Bash script executable.
    # Sending command to run.
    phrase_b1 = "chmod +x {}".format(gene_bash_script)
    phrase_b2 = "chmod --version"
    title_1 = "Making a Bash script executable:"
    title_2 = "Version of chmod:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
    # Sending command to run.
    phrase_b1 = "bash {}".format(gene_bash_script)
    phrase_b2 = "bash --version"
    title_1 = "Running bash:"
    title_2 = "Version of bash:"
    capture_status = True
    shell_status = True
    pr_status = False
    command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    # Check the output file.
    if os.path.getsize(output_path_genepred_faa) != 0:
        split_size = 70
        new_file_fgs_formatted = open(output_fgs_protein_formatted_1, "w")
        with open("{}".format(output_path_genepred_faa)) as fgs_output_lines:
            for line in fgs_output_lines:
                line = line.rstrip("\n")
                if line[0] == ">":
                    new_file_fgs_formatted.write("{}\n".format(line))
                else:
                    fasta_splited = [line[ssi:ssi + split_size] for ssi in range(0, len(line), split_size)]
                    for fasta_line in fasta_splited:
                        new_file_fgs_formatted.write("{}\n".format(fasta_line))
        new_file_fgs_formatted.close()


def gene_annotation(output_path_genepred, output_path_genepred_ffn, output_fgs_protein_formatted_2, output_final_contigs_formated, gene_contig_dist_path, gene_info_file, genetic_code):
    print("\nExtracting information from the predicted gene sequences...")
    # If there is a file with the gene sequence from FragGeneScanRs then:
    # For each protein ID of the putative proteins parse the file and find its gene sequence.
    # Then check for starting (including those of bacterial type) and end codons in the sequence and also check the
    # sequence for the number of rare codos and ideally run it across a program/software that scores the gene sequence
    # based on its expression capability.
    gene_sequences_dict = {}
    if os.path.exists(output_path_genepred) and os.path.exists(output_path_genepred_ffn):
        gene_lines = supportive_functions.read_file(output_path_genepred_ffn)
        for i in range(0, len(gene_lines)):
            line = gene_lines[i]
            if line[0] == ">":
                protein_id = line[1:]
                gene_sequence = gene_lines[i+1]
                gene_sequences_dict[protein_id] = [gene_sequence]

    # Check for start and end codons for any gene sequence found.
    # For the next item after the gene sequence, if a start codon is present then 1 is placed, otherwise 0.
    # Similarly for the end codon.
    dict_genes = {}
    if gene_sequences_dict:
        for key in gene_sequences_dict.keys():
            if gene_sequences_dict[key][0][0:3] == "ATG" or gene_sequences_dict[key][0][0:3] == "GTG" or gene_sequences_dict[key][0][0:3] == "GTC":
                gene_sequences_dict[key].append("YES")
            else:
                gene_sequences_dict[key].append("NO")
            if gene_sequences_dict[key][0][-3:] == "TAA" or gene_sequences_dict[key][0][-3:] == "TAG" or gene_sequences_dict[key][0][-3:] == "TGA":
                gene_sequences_dict[key].append("YES")
            else:
                gene_sequences_dict[key].append("NO")
        new_file = open(gene_info_file, "w")
        for key in gene_sequences_dict.keys():
            new_file.write("{}\t{}\t{}\t{}\n".format(key, gene_sequences_dict[key][0], gene_sequences_dict[key][1], gene_sequences_dict[key][2]))
            dict_genes[key] = [gene_sequences_dict[key][0], gene_sequences_dict[key][1], gene_sequences_dict[key][2]]
        new_file.close()

    # Translate gene sequences
    pr_dict = {}
    if gene_sequences_dict:
        for key in gene_sequences_dict.keys():
            gene_seq = gene_sequences_dict[key][0]
            put_protein = translation(gene_seq, genetic_code)
            pr_dict[key] = put_protein

    # If proteines have been encoded by any genes, they are written in a file.
    if pr_dict:
        split_size = 70
        new_file_gc_formatted = open(output_fgs_protein_formatted_2, "w")
        for key in pr_dict.keys():
            pr_sqs = pr_dict[key]
            new_file_gc_formatted.write(">{}\n".format(key))
            fasta_splited = [pr_sqs[ssi:ssi + split_size] for ssi in range(0, len(pr_sqs), split_size)]
            for fasta_line in fasta_splited:
                new_file_gc_formatted.write("{}\n".format(fasta_line))
        new_file_gc_formatted.close()

    # Compute the distance of the gene's edges from the edges of its contig(s).
    # For megahit the headers of the contigs contain their lengths.
    contig_len_dict = None
    contig_gene_dist_dict = {}
    if os.path.exists(output_final_contigs_formated):
        contig_lines = supportive_functions.read_file(output_final_contigs_formated)
        header = None
        header_contig = None
        contig_seq = ""
        first_header = True
        contig_len_dict = {}
        for line in contig_lines:
            if line[0] == ">":
                if not first_header:
                    contig_len_dict[header_contig] = len(contig_seq)
                    contig_seq = ""
                header = line[1:]
                first_header = False
                # Process the header: "k141_49635 flag=1 multi=3.0000 len=461"
                # Keep only the contig's code without k: "141_49635"
                if " " in header:
                    header_splited = header.split(" ")
                    header_contig = header_splited[0]
                else:
                    header_contig = header
            else:
                contig_seq = "{}{}".format(contig_seq, line)
        contig_len_dict[header_contig] = len(contig_seq)

    # Find the code of the contig in which the gene was found. If the code of the contig has been corresponded to its length, then compare the length of the gene with length of the
    # contig to find the ditance of the gene from the start and end of the contig.
    for key in pr_dict.keys():
        # Example: key = contig23830_303_It-6_meta_12_197_- or key = contig24087_1789_It-6_meta_1062_1211_+
        if "_" in key:
            key_splited = key.split("_")
            origin_contig = "_".join(key_splited[0:2])
            contig_length = None
            if origin_contig in contig_len_dict.keys():
                contig_length = contig_len_dict[origin_contig]
            gene_start = int(key_splited[-3])
            gene_end = int(key_splited[-2])
            start_dist = gene_start
            end_dist = contig_length - gene_end
            contig_gene_dist_dict[key] = [start_dist, end_dist]
        else:
            contig_gene_dist_dict[key] = ["-", "-"]

    # Write the information in a file.
    new_file = open(gene_contig_dist_path, "w")
    for key in contig_gene_dist_dict.keys():
        new_file.write("{}\t{}\t{}\n".format(key, contig_gene_dist_dict[key][0], contig_gene_dist_dict[key][1]))
    new_file.close()
    return gene_sequences_dict, dict_genes, contig_gene_dist_dict