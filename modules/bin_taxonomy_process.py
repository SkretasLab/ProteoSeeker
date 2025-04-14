import re
import os
import csv
import copy
# ProteoSeeker modules
import command_process
import supportive_functions


def bin_tax_analysis(fpr_db_fasta, blastp_doms_info_nr_file, dict_contigs_bins, binning_results_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prstax, binphylo_path_binstax, binphylo_path_binstax_max, binphylo_path_prsids, binphylo_path_binstax_names, binphylo_path_binstax_max_names, binphylo_freq_taxids_path, binphylo_maxfreq_taxids_path, cd_hit_results_fasta_path, taxonkit_env, taxonkit_path, tax_freq_bash_script, tax_maxfreq_bash_script, taxoknit_freq_version_path, taxoknit_freq_stdoe_path, taxoknit_maxfreq_version_path, taxoknit_maxfreq_stdoe_path, thread_num, conda_sh_path, freq_taxids_path, maxfreq_taxids_path, taxonkit_freq_line_bash_script, taxonkit_maxfreq_line_bash_script, taxonkit_freq_line_version_path, taxonkit_freq_line_stdoe_path, taxonkit_maxfreq_line_version_path, taxonkit_maxfreq_line_stdoe_path, freq_lineage_path, maxfreq_lineage_path, freq_lineage_form_path, maxfreq_lineage_form_path, csvtk_freq_version_path, csvtk_freq_stdoe_path, csvtk_maxfreq_version_path, csvtk_maxfreq_stdoe_path, family_to_profile_phylo_dict, hmmer_enz_domains_all_proteins_phylo, family_profile_path, input_log_file, output_log_file):
    print("\nPerforming a taxonomic analysis of the bins...")
    # The process is:
    # 1. Run HMMER spec for the proteins of CD-HIT for each database selected for phylogenetic analysis (these databases are already created).
    # 2. For each identified hit, run the corresponding nr database (these databases are already created).
    # 3. At this stage use the results of each of these proteins (the bins do not change) and sum them up to concluce the taxonomy(ies) of each bin and of each protein.
    # The pnr database must be parsed and collect the correspondance between nr ACs and their organisms.
    dict_pnr_tax = {}
    pattern_acc = re.compile(r'^>(\S+)')
    pattern_tax_1 = re.compile(r'\[(.*?)\]')
    pattern_tax_2 = re.compile(r'Tax=(.*?) TaxID=')
    if os.path.exists(fpr_db_fasta):
        with open(fpr_db_fasta) as pnr_lines:
            for line in pnr_lines:
                line = line.rstrip("\n")
                if line[0] == ">":
                    nr_acc = None
                    nr_tax = None
                    # Search for the accession number.
                    result_acc = pattern_acc.search(line)
                    if result_acc:
                        # The taxa name inside the the first brackets is retained.
                        nr_acc = result_acc.group(1)
                    # Pattern 1
                    result_tax_1 = re.findall(pattern_tax_1, line)
                    if result_tax_1:
                        nr_tax = []
                        # Count the species related to the protein sequence.
                        nr_tax_dict = {}
                        for taxum in result_tax_1:
                            if taxum not in nr_tax_dict.keys():
                                nr_tax_dict[taxum] = 1
                            else:
                                nr_tax_dict[taxum] += 1
                        # If only one species was found, keep that species, otherwise keep the one(s) with the highest count.
                        if len(list(nr_tax_dict.keys()))  == 1:
                            nr_tax = [result_tax_1[0]]
                        else:
                            max_count = 0
                            for cur_count in nr_tax_dict.values():
                                if cur_count > max_count:
                                    max_count = cur_count
                            for cut_tax in nr_tax_dict.keys():
                                if nr_tax_dict[cut_tax] == max_count:
                                    nr_tax.append(cut_tax)
                    # Pattern 2
                    result_tax_2 = re.search(pattern_tax_2, line)
                    if result_tax_2:
                        nr_tax = []
                        taxum = result_tax_2.group(1)
                        nr_tax.append(taxum)
                    # Result
                    if nr_acc is not None:
                        if nr_tax is not None:
                            if nr_tax:
                                dict_pnr_tax[nr_acc] = nr_tax
    else:
        return {}, {}, {}, {}, {}, {}

    # NR results
    dict_nr = {}
    if os.path.exists(blastp_doms_info_nr_file):
        protein_acc = None
        min_evalue = None
        header = True
        with open(blastp_doms_info_nr_file) as lines_nr:
            for line_nr in lines_nr:
                if header:
                    header = False
                    continue
                line_nr = line_nr.rstrip("\n")
                splited_nr = line_nr.split("\t")
                protein_acc = splited_nr[0]
                e_value = float(splited_nr[3])
                # Line: k141_39729_24847_26913_-    cbw76918.1      43.0    1.70e-159   473
                # Line: k141_39729_24847_26913_-    wp_014905841.1  44.2    1.20e-156   464
                # Line: k141_29822_4122_5963_-      cab1404101.1    58.3    4.40e-246   686
                if protein_acc not in dict_nr.keys():
                    dict_nr[protein_acc] = line_nr
                    min_evalue = e_value
                else:
                    if e_value < min_evalue:
                        dict_nr[protein_acc] = line_nr
                        min_evalue = e_value
    else:
        print("\nThe path to the file with the nr information, for the phylogenetic analysis of the bins, is wrong. The phylogenetic analysis based on the bins is not performed.")
    
    # Binning results
    dict_prs_bins = {}
    dict_bins_prs = {}
    if not dict_contigs_bins:
        if os.path.exists(binning_results_path):
            with open(binning_results_path) as und_lines:
                tsv_lines = csv.reader(und_lines, delimiter="\t", quotechar='"')
                for line in tsv_lines:
                    contig_id = line[0]
                    bin_id = line[1]
                    if bin_id[:5] == "group":
                        bin_id = bin_id[5:]
                    dict_contigs_bins[contig_id] = bin_id
        dict_prs_contigs = {}
        if os.path.exists(output_path_genepred) and os.path.exists(output_path_genepred_faa):
            protein_lines = supportive_functions.read_file(output_path_genepred_faa)
            for line in protein_lines:
                if line[0] == ">":
                    protein_id = line[1:]
                    if "_" in protein_id:
                        protein_id_splited = protein_id.split("_")
                        origin_contig = "_".join(protein_id_splited[0:2])
                        dict_prs_contigs[protein_id] = origin_contig
        if dict_contigs_bins and dict_prs_contigs:
            for pr_id in dict_prs_contigs.keys():
                contig_id = dict_prs_contigs[pr_id]
                if contig_id in dict_contigs_bins.keys():
                    bin_id = dict_contigs_bins[contig_id]
                else:
                    bin_id = "-"
                dict_prs_bins[pr_id] = bin_id
                if bin_id not in dict_bins_prs.keys():
                    dict_bins_prs[bin_id] = [pr_id]
                else:
                    dict_bins_prs[bin_id].append(pr_id)
    if (not dict_prs_bins) or (not dict_bins_prs):
        return {}, {}, {}, {}, {}, {}

    # dict_pnr_tax:
    # wp_088159341.1: achromobacter xylosoxidans
    # dict_nr:
    # k141_10192_17880_18695_-: k141_10192_17880_18695_-\tmbh0196785.1\t59.0\t1.58e-55\t187
    # dict_prs_bins:
    # k141_49593_60106_60408_+: 2
    # k141_49596_2_409_+: -
    # The blastp_info_comb_nr_file contains the combined information from the proteins wihtout any domain of interest that scored e-values below the threshold against pnr and from the proteins
    # each which had at least one hit against the phmm database for the selected protein family.
    # Parse each bin:
    # Parse each of its proteins:
    # Check if the protein in part of the proteins in the file blastp_info_comb_nr_file.
    # If no, conitnue on.
    # If yes, collect the taxonomy of the protein.
    # If no such protein was found, then no taxonomy is assigned to the bin.
    # If at least one such protein was found for the bin then:
    # Assign to the bin, the the taxonomy of the highest frequency based on the proteins parse above.
    # If two or more proteins were found for the bin and two or more taxonomies had the highest frequency, then assign to the bin all taxonomies.
    # All genes and proteins of the bin are assigned the taxonomy(ies) of the bin, if any.
    # ---
    # Sort the bins in ascending order of the bin IDs.
    # Create a list of the bin IDs.
    dict_bins_prs_keys = list(dict_bins_prs.keys())
    # The bin IDs, as keys of a dictionary, are already unique.
    # Remove the "-" (undentified bin label) from the list.
    if "-" in dict_bins_prs_keys:
        dict_bins_prs_keys.remove("-")
    # Sort the keys in ascending order.
    dict_bins_prs_keys_sorted = sorted(dict_bins_prs_keys, key=int)
    # The sorted dictionary.
    dict_bins_prs_sorted = {}
    for bin_id in dict_bins_prs_keys_sorted:
        for key_bin in dict_bins_prs.keys():
            if bin_id == key_bin:
                dict_bins_prs_sorted[key_bin] = copy.deepcopy(dict_bins_prs[key_bin])
    # Write the bin and their protein IDs in a file.
    new_file_bin_proteins = open(binphylo_path_prsids, "w")
    for binid in dict_bins_prs_sorted.keys():
        for item in dict_bins_prs_sorted[binid]:
            new_file_bin_proteins.write("{}\t{}\n".format(binid, item))
    new_file_bin_proteins.close()

    # Dictionary with the information of bins - proteins - taxonomies.
    bins_prs_tax_dict = {}
    # Dictinoary that holds the frequencies fo the taxonomies for each bin.
    bin_tax_dict = {}
    # For a bin...
    for key_bin in dict_bins_prs_sorted.keys():
        if key_bin == "-":
            continue
        # Add the bin the dictionary of bins and taxonomies.
        bins_prs_tax_dict[key_bin] = {}
        bin_tax_dict[key_bin] = {}
        # For a protein...
        for pr_acc in dict_bins_prs_sorted[key_bin]:
            # Check if protein in the proteins with information from nr.
            if pr_acc in dict_nr.keys():
                # Analyze the information from the pnr database for the protein.
                nr_info = dict_nr[pr_acc]
                # Add the protien in the dictionary.
                bins_prs_tax_dict[key_bin][pr_acc] = {}
                # Collect the ID of the hit in the pnr/nr database.
                nr_info_splited = nr_info.split("\t")
                pr_nr_acc = nr_info_splited[1]
                # Collect the taxonomy of the accession number of the pnr/nr database.
                if pr_nr_acc in dict_pnr_tax.keys():
                    # pr_nr_tax is a list of highest occuring species from each protein
                    pr_nr_tax = dict_pnr_tax[pr_nr_acc]
                    for tax_item in pr_nr_tax:
                        # Bins - Proteins - Taxonmies Dicitonary
                        if tax_item not in bins_prs_tax_dict[key_bin][pr_acc].keys():
                            bins_prs_tax_dict[key_bin][pr_acc][tax_item] = 1
                        else:
                            bins_prs_tax_dict[key_bin][pr_acc][tax_item] += 1
                        # Bins - Taxonmies Dicitonary
                        if tax_item not in bin_tax_dict[key_bin].keys():
                            bin_tax_dict[key_bin][tax_item] = 1
                        else:
                            bin_tax_dict[key_bin][tax_item] += 1

    # Protein IDs to taxa with a flag for direct or indirect determination.
    # Select the highest occuring species from the bin_tax_dict.
    bin_tax_max_dict = {}
    for binid in bin_tax_dict.keys():
        bin_tax_max_dict[binid] = {}
        # Find the maximum taxonomy frequency.
        tax_freq_max = 0
        for species_id in bin_tax_dict[binid].keys():
            tax_freq_cur = bin_tax_dict[binid][species_id]
            # If the frequency of the current taxnomy is higher than the currently max one, replace the existing taxonomies for the current bin.
            # If the frequency of the current taxonomy is equal to the currently max one, add the existing taxonmy to the current bin.
            # Otherwise, do nothing.
            if tax_freq_cur > tax_freq_max:
                bin_tax_max_dict[binid] = {}
                bin_tax_max_dict[binid][species_id] = tax_freq_cur
                tax_freq_max = tax_freq_cur
            elif tax_freq_cur == tax_freq_max:
                bin_tax_max_dict[binid][species_id] = tax_freq_cur

    # Dictionary to write the information in the annotation files.
    bin_group_tax_dict = {}
    if not bin_group_tax_dict:
        for bin_id in bin_tax_max_dict.keys():
            bin_group_tax_dict[bin_id] = {
                "max_freq": None,
                "protein_ids_d": [],
                "protein_ids_i": [],
                "species": []
            }
            for species_id in bin_tax_max_dict[bin_id].keys():
                if bin_group_tax_dict[bin_id]["max_freq"] is None:
                    max_freq = bin_tax_max_dict[bin_id][species_id]
                    bin_group_tax_dict[bin_id]["max_freq"] = max_freq
                bin_group_tax_dict[bin_id]["species"].append(species_id)
                # For each protein of the bin, if the protein is present in bins_prs_tax_dict[bin_id] then it has immediate taxonomy inference otherwise intemediate.
                bin_protens_with_tax = list(bins_prs_tax_dict[bin_id].keys())
                for protein_id in dict_prs_bins.keys():
                    if dict_prs_bins[protein_id] == bin_id:
                        if protein_id in bin_protens_with_tax:
                            bin_group_tax_dict[bin_id]["protein_ids_d"].append(protein_id)
                        else:
                            bin_group_tax_dict[bin_id]["protein_ids_i"].append(protein_id)

    # Binned proteins FASTA file.
    bin_pr_accs = set()
    for key_bin in dict_bins_prs_sorted.keys():
        for pr_ac in dict_bins_prs_sorted[key_bin]:
            bin_pr_accs.add(pr_ac)

    # Protein - Bin
    new_file_bin_phylo_prs = open(binphylo_path_prstax, "w")
    for key_1 in bins_prs_tax_dict.keys():
        for key_2 in bins_prs_tax_dict[key_1].keys():
            for key_3 in bins_prs_tax_dict[key_1][key_2].keys():
                new_file_bin_phylo_prs.write("{}\t{}\t{}\t{}\n".format(key_1, key_2, key_3, bins_prs_tax_dict[key_1][key_2][key_3]))
    new_file_bin_phylo_prs.close()
    # Bin - Taxon - Frequency
    new_file_bin_phylo_bins = open(binphylo_path_binstax, "w")
    for key_1 in bin_tax_dict.keys():
        for key_2 in bin_tax_dict[key_1].keys():
            new_file_bin_phylo_bins.write("{}\t{}\t{}\n".format(key_1, key_2, bin_tax_dict[key_1][key_2]))
    new_file_bin_phylo_bins.close()
    # Taxon name - Frequency
    new_file_name_freq = open(binphylo_path_binstax_names, "w")
    for key_1 in bin_tax_dict.keys():
        for key_2 in bin_tax_dict[key_1].keys():
            new_file_name_freq.write("{}\n".format(key_2))
    new_file_name_freq.close()
    # Bin - Taxon - Max frequency
    new_file_bin_phylo_bins_max = open(binphylo_path_binstax_max, "w")
    for key_1 in bin_tax_max_dict.keys():
        for key_2 in bin_tax_max_dict[key_1].keys():
            new_file_bin_phylo_bins_max.write("{}\t{}\t{}\n".format(key_1, key_2, bin_tax_max_dict[key_1][key_2]))
    new_file_bin_phylo_bins_max.close()
    # Taxon name - Max frequency
    new_file_name_maxfreq = open(binphylo_path_binstax_max_names, "w")
    for key_1 in bin_tax_max_dict.keys():
        for key_2 in bin_tax_max_dict[key_1].keys():
            new_file_name_maxfreq.write("{}\n".format(key_2))
    new_file_name_maxfreq.close()

    # Finding taxonomy IDs
    if os.path.exists(binphylo_path_binstax_names):
        # Command
        if os.path.exists(taxonkit_path):
            phrase_1 = "\"{}\" name2taxid \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(taxonkit_path, binphylo_path_binstax_names, binphylo_freq_taxids_path, thread_num, taxoknit_freq_stdoe_path)
            phrase_2 = "\"{}\" -h > \"{}\"".format(taxonkit_path, taxoknit_freq_version_path)
        else:
            phrase_1 = "taxonkit name2taxid \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(binphylo_path_binstax_names, binphylo_freq_taxids_path, thread_num, taxoknit_freq_stdoe_path)
            phrase_2 = "taxonkit -h > \"{}\"".format(taxoknit_freq_version_path)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(tax_freq_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if taxonkit_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(taxonkit_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        if taxonkit_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(tax_freq_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(tax_freq_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    if os.path.exists(binphylo_path_binstax_max_names):
        # Command
        if os.path.exists(taxonkit_path):
            phrase_1 = "\"{}\" name2taxid \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(taxonkit_path, binphylo_path_binstax_max_names, binphylo_maxfreq_taxids_path, thread_num, taxoknit_maxfreq_stdoe_path)
            phrase_2 = "\"{}\" -h > \"{}\"".format(taxonkit_path, taxoknit_maxfreq_version_path)
        else:
            phrase_1 = "taxonkit name2taxid \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(binphylo_path_binstax_max_names, binphylo_maxfreq_taxids_path, thread_num, taxoknit_maxfreq_stdoe_path)
            phrase_2 = "taxonkit -h > \"{}\"".format(taxoknit_maxfreq_version_path)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(tax_maxfreq_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if taxonkit_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(taxonkit_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        if taxonkit_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(tax_maxfreq_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(tax_maxfreq_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    # Collected the unique taxonomy IDs
    if os.path.exists(binphylo_freq_taxids_path):
        unique_taxids = []
        freq_name_lines = supportive_functions.read_file(binphylo_freq_taxids_path)
        for line in freq_name_lines:
            line_splited = line.split("\t")
            taxid = line_splited[1]
            if taxid not in unique_taxids:
                unique_taxids.append(taxid)
        freq_taxids_file = open(freq_taxids_path, "w")
        for item in unique_taxids:
            freq_taxids_file.write("{}\n".format(item))
        freq_taxids_file.close()
    if os.path.exists(binphylo_maxfreq_taxids_path):
        unique_taxids = []
        freq_name_lines = supportive_functions.read_file(binphylo_maxfreq_taxids_path)
        for line in freq_name_lines:
            line_splited = line.split("\t")
            taxid = line_splited[1]
            if taxid not in unique_taxids:
                unique_taxids.append(taxid)
        freq_taxids_file = open(maxfreq_taxids_path, "w")
        for item in unique_taxids:
            freq_taxids_file.write("{}\n".format(item))
        freq_taxids_file.close()

    # Find the lineage
    if os.path.exists(freq_taxids_path):
        # Command
        if os.path.exists(taxonkit_path):
            phrase_1 = "\"{}\" lineage \"{}\" -o \"{}\" -j {} &> \"{}\"".format(taxonkit_path, freq_taxids_path, freq_lineage_path, thread_num, taxonkit_freq_line_stdoe_path)
            phrase_2 = "\"{}\" -h > \"{}\"".format(taxonkit_path, taxonkit_freq_line_version_path)
            phrase_3 = "\"{}\" lineage \"{}\" -j {} | csvtk pretty -Ht -x ';' -W 70 -S bold -o \"{}.txt\" &> \"{}\"".format(taxonkit_path, freq_taxids_path, thread_num, freq_lineage_form_path, csvtk_freq_stdoe_path)
            phrase_4 = "csvtk -h > {}".format(csvtk_freq_version_path)
        else:
            phrase_1 = "taxonkit lineage \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(freq_taxids_path, freq_lineage_path, thread_num, taxonkit_freq_line_stdoe_path)
            phrase_2 = "taxonkit -h > \"{}\"".format(taxonkit_freq_line_version_path)
            phrase_3 = "taxonkit lineage \"{}\" -j {} | csvtk pretty -Ht -x ';' -W 70 -S bold -o \"{}\" &> \"{}\"".format(freq_taxids_path, thread_num, freq_lineage_form_path, csvtk_freq_stdoe_path)
            phrase_4 = "csvtk -h > \"{}\"".format(csvtk_freq_version_path)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(taxonkit_freq_line_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if taxonkit_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(taxonkit_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        new_file_bash.write("{}\n".format(phrase_3))
        new_file_bash.write("{}\n".format(phrase_4))
        if taxonkit_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(taxonkit_freq_line_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(taxonkit_freq_line_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    if os.path.exists(maxfreq_taxids_path):
        # Command
        if os.path.exists(taxonkit_path):
            phrase_1 = "\"{}\" lineage \"{}\" -o \"{}\" -j {} &> \"{}\"".format(taxonkit_path, maxfreq_taxids_path, maxfreq_lineage_path, thread_num, taxonkit_maxfreq_line_stdoe_path)
            phrase_2 = "\"{}\" -h > \"{}\"".format(taxonkit_path, taxonkit_maxfreq_line_version_path)
            phrase_3 = "\"{}\" lineage \"{}\" -j {} | csvtk pretty -Ht -x ';' -W 70 -S bold -o \"{}.txt\" &> \"{}\"".format(taxonkit_path, maxfreq_taxids_path, thread_num, maxfreq_lineage_form_path, csvtk_maxfreq_stdoe_path)
            phrase_4 = "csvtk -h > \"{}\"".format(csvtk_maxfreq_version_path)
        else:
            phrase_1 = "taxonkit lineage \"{}\" -r -o \"{}\" -j {} &> \"{}\"".format(maxfreq_taxids_path, maxfreq_lineage_path, thread_num, taxonkit_maxfreq_line_stdoe_path)
            phrase_2 = "taxonkit -h > \"{}\"".format(taxonkit_maxfreq_line_version_path)
            phrase_3 = "taxonkit lineage \"{}\" -j {} | csvtk pretty -Ht -x ';' -W 70 -S bold -o \"{}\" &> \"{}\"".format(maxfreq_taxids_path, thread_num, maxfreq_lineage_form_path, csvtk_maxfreq_stdoe_path)
            phrase_4 = "csvtk -h > \"{}\"".format(csvtk_maxfreq_version_path)
        # Create the Bash script.
        # Four cases: 1: Conda environment and path needed for the script. 2: Conda environment and no path needed for the script. 3: No conda environment and path needed for the script. 4: No conda environment and no path needed for the script.
        new_file_bash = open(taxonkit_maxfreq_line_bash_script, "w")
        phrase = "#!/bin/bash"
        new_file_bash.write("{}\n".format(phrase))
        if taxonkit_env:
            phrase_s = "source \"{}\"".format(conda_sh_path)
            phrase_a = "conda activate \"{}\"".format(taxonkit_env)
            new_file_bash.write("{}\n".format(phrase_s))
            new_file_bash.write("{}\n".format(phrase_a))
        new_file_bash.write("{}\n".format(phrase_1))
        new_file_bash.write("{}\n".format(phrase_2))
        new_file_bash.write("{}\n".format(phrase_3))
        new_file_bash.write("{}\n".format(phrase_4))
        if taxonkit_env:
            phrase = "conda deactivate"
            new_file_bash.write("{}\n".format(phrase))
        new_file_bash.close()
        # Making the Bash script executable.
        # Sending command to run.
        phrase_b1 = "chmod +x {}".format(taxonkit_maxfreq_line_bash_script)
        phrase_b2 = "chmod --version"
        title_1 = "Making a Bash script executable:"
        title_2 = "Version of chmod:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)
        # Sending command to run.
        phrase_b1 = "bash {}".format(taxonkit_maxfreq_line_bash_script)
        phrase_b2 = "bash --version"
        title_1 = "Running bash:"
        title_2 = "Version of bash:"
        capture_status = True
        shell_status = True
        pr_status = False
        command_process.command_run(phrase_b1, phrase_b2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file)

    # Determine the frequency of profile usage for each family in locating the respective domains at the binned proteins.
    # The binned protein IDs hace been collected in bin_pr_accs.
    # The results of hmmscan of phylo profiles against the protein from cd-hit are in the hmmer_enz_domains_all_proteins_phylo file.
    family_to_profile_phylo_sum_dict = {}
    if bin_pr_accs and os.path.exists(hmmer_enz_domains_all_proteins_phylo):
        with open(hmmer_enz_domains_all_proteins_phylo) as hmmer_phylo_lines:
            for line in hmmer_phylo_lines:
                line_splited = line.split("\t")
                source_pr_acc = line_splited[0]
                if source_pr_acc in bin_pr_accs:
                    profile_id = line_splited[4]
                    # Remove any suffix after the dot, if a dot is present.
                    if "." in profile_id:
                        profile_id_splited = profile_id.split(".")
                        profile_id = profile_id_splited[0]
                    for key_fam in family_to_profile_phylo_dict.keys():
                        if key_fam not in family_to_profile_phylo_sum_dict.keys():
                            family_to_profile_phylo_sum_dict[key_fam] = 0
                        if profile_id in family_to_profile_phylo_dict[key_fam].keys():
                            family_to_profile_phylo_dict[key_fam][profile_id] += 1
                            family_to_profile_phylo_sum_dict[key_fam] += 1
    family_profile_file = open(family_profile_path, "w")
    family_profile_file.write("Protein family\tProfile ID\tFrequency\n")
    for key_fam in family_to_profile_phylo_dict.keys():
        for key_prof in family_to_profile_phylo_dict[key_fam].keys():
            fam_prof_freq = family_to_profile_phylo_dict[key_fam][key_prof]
            family_profile_file.write("{}\t{}\t{}\n".format(key_fam, key_prof, fam_prof_freq))
    if family_to_profile_phylo_sum_dict:
        family_profile_file.write("\n")
        family_profile_file.write("Protein family\tSum of profile hits\n")
        for key_fam in family_to_profile_phylo_sum_dict.keys():
            fam_prof_sum = family_to_profile_phylo_sum_dict[key_fam]
            family_profile_file.write("{}\t{}\n".format(key_fam, fam_prof_sum))
    family_profile_file.close()

    return dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict