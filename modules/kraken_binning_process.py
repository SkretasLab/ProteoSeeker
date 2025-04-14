import os
import copy
# ProteoSeeker modules
import supportive_functions


def kraken_binning(taxonomy_route, contigs_to_reads_dict, read_to_species_dict, taxid_to_species_dict, kraken_taxname_path, kraken_reads_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prsids, binned_ctb_path, binned_btc_path, binned_taxa_path, kraken_bin_info_path):
    print("\nPerforming binning based on the results of kraken2...")
    dict_contigs_bins = {}
    dict_bins_contigs = {}
    bin_group_tax_dict = {}
    if not taxid_to_species_dict:
        with open(kraken_taxname_path) as taxname_lines:
            for line in taxname_lines:
                line = line.rstrip("\n")
                line_splited = line.split("\t")
                taxid = int(line_splited[0])
                taxname = line_splited[1]
                taxid_to_species_dict[taxid] = taxname
    if taxonomy_route == 1:
        # For each contig:
        # Collect all its aligned reads. Get the species assigned to each read. Count each species
        # for the contig based on its reads. Keep the species of the highest count for each contig.
        # Based on the contigs that were assigned one species, group them based on their species of the highest counts.
        # Each group is a bin. Each bin has one species. The species of the bin is assigned to all genes of the bin.
        if not read_to_species_dict:
            with open(kraken_reads_path) as readtax_lines:
                for line in readtax_lines:
                    line = line.rstrip("\n")
                    line_splited = line.split("\t")
                    read_id = line_splited[0]
                    taxonid = int(line_splited[1])
                    thr_state = line_splited[2]
                    read_to_species_dict[read_id] = [taxonid, thr_state]
        # Contig -> Read -> Species -> Thr state
        contig_to_species = {}
        for key_contig in contigs_to_reads_dict.keys():
            if key_contig not in contig_to_species.keys():
                contig_to_species[key_contig] = {}
            if key_contig in contigs_to_reads_dict.keys():
                for read_id in contigs_to_reads_dict[key_contig][1:]:
                    if read_id in read_to_species_dict.keys():
                        read_species = read_to_species_dict[read_id][0]
                        thr_state = int(read_to_species_dict[read_id][1])
                        if thr_state == 1:
                            if read_species not in contig_to_species[key_contig].keys():
                                contig_to_species[key_contig][read_species] = 1
                            else:
                                contig_to_species[key_contig][read_species] += 1
        # Filter the species of maximum count for each contig.
        contig_to_species_max = {}
        for key_contig in contig_to_species.keys():
            max_count = 0
            for key_species in contig_to_species[key_contig].keys():
                species_count = contig_to_species[key_contig][key_species]
                if species_count > max_count:
                    max_count = species_count
            for key_species in contig_to_species[key_contig].keys():
                species_count = contig_to_species[key_contig][key_species]
                if max_count == species_count:
                    if key_contig not in contig_to_species_max:
                        contig_to_species_max[key_contig] = {}
                    contig_to_species_max[key_contig][key_species] = max_count
        # Filter the contigs with only one species of maximum count.
        contig_to_species_max_single = {}
        for key_contig in contig_to_species_max.keys():
            species_list = list(contig_to_species_max[key_contig].keys())
            species_num = len(species_list)
            if species_num == 1:
                contig_to_species_max_single[key_contig] = {}
                for key_species in contig_to_species_max[key_contig].keys():
                    count = contig_to_species_max[key_contig][key_species]
                    contig_to_species_max_single[key_contig][key_species] = count
        # Bin the contigs.
        bin_id_init = 0
        species_to_binids_dict = {}
        for key_contig in contig_to_species_max_single.keys():
            for key_species in contig_to_species_max_single[key_contig].keys():
                if key_species not in species_to_binids_dict.keys():
                    bin_id_init += 1
                    bin_id = bin_id_init
                    species_to_binids_dict[key_species] = bin_id
                    bin_id_str = str(bin_id)
                    dict_bins_contigs[bin_id_str] = [key_contig]
                else:
                    bin_id = species_to_binids_dict[key_species]
                    bin_id_str = str(bin_id)
                    dict_bins_contigs[bin_id_str].append(key_contig)
                dict_contigs_bins[key_contig] = bin_id_str
        # Binning results
        dict_prs_bins = {}
        dict_bins_prs = {}
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
        # Create a list of the bin IDs.
        dict_bins_prs_keys = list(dict_bins_prs.keys())
        if "-" in dict_bins_prs_keys:
            dict_bins_prs_keys.remove("-")
        dict_bins_prs_keys_sorted = sorted(dict_bins_prs_keys, key=int)
        # The sorted dictionary.
        dict_bins_prs_sorted = {}
        for bin_id in dict_bins_prs_keys_sorted:
            for key_bin in dict_bins_prs.keys():
                if bin_id == key_bin:
                    dict_bins_prs_sorted[key_bin] = copy.deepcopy(dict_bins_prs[key_bin])
        # Store proteins based on their bin IDs.
        prid_to_bin_dict = {}
        for protein_id, bin_id in dict_prs_bins.items():
            if bin_id not in prid_to_bin_dict:
                prid_to_bin_dict[bin_id] = [protein_id]
            else:
                prid_to_bin_dict[bin_id].append(protein_id)
        # In contrast with MetaBinner and COMEBin the taxonomy of the proteins in the bins created based on kraken2 are inferred all by the same means.
        # Dictionary to write the information in the annotation files.
        # In this case max_freq is the highest number of reads towards the primary species (the species of the bin) of any of the contigs of the bin.
        bin_group_tax_dict = {}
        for key_contig in contig_to_species_max_single.keys():
            bin_id = dict_contigs_bins[key_contig]
            if bin_id not in bin_group_tax_dict.keys():
                bin_group_tax_dict[bin_id] = {
                    "max_freq": None,
                    "protein_ids_k": [],
                    "species": []
                }
            for species_id in contig_to_species_max_single[key_contig].keys():
                if bin_group_tax_dict[bin_id]["max_freq"] is None:
                    max_freq = contig_to_species_max_single[key_contig][species_id]
                    bin_group_tax_dict[bin_id]["max_freq"] = max_freq
                species_name = taxid_to_species_dict[species_id]
                species_comb = "{}||{}".format(species_name, species_id)
                if species_comb not in bin_group_tax_dict[bin_id]["species"]:
                    bin_group_tax_dict[bin_id]["species"].append(species_comb)
                if bin_id in prid_to_bin_dict.keys():
                    bin_group_tax_dict[bin_id]["protein_ids_k"] = prid_to_bin_dict[bin_id]
        print("\nWriting information related to kraken2 binning.")
        # Write the binned contigs according to conigs
        binned_ctb_file = open(binned_ctb_path, "w")
        for key_contig in dict_contigs_bins.keys():
            bin_id = dict_contigs_bins[key_contig]
            binned_ctb_file.write("{}\t{}\n".format(key_contig, bin_id))
        binned_ctb_file.close()
        # Write the binned contigs according to bins
        bin_ids_list_dupl = list(dict_contigs_bins.values())
        bin_ids_list = []
        for item in bin_ids_list_dupl:
            if item not in bin_ids_list:
                bin_ids_list.append(item)
        bin_ids_list = sorted(bin_ids_list, key=int)
        binned_btc_file = open(binned_btc_path, "w")
        for bin_id in bin_ids_list:
            for key_contig in dict_contigs_bins.keys():
                if bin_id == dict_contigs_bins[key_contig]:
                    binned_btc_file.write("{}\t{}\n".format(bin_id, key_contig))
        binned_btc_file.close()
        # Write the binned contigs
        # Write this file an use it in info collection properly.
        binned_taxa_file = open(binned_taxa_path, "w")
        for key_bin in bin_group_tax_dict.keys():
            max_freq = bin_group_tax_dict[key_bin]["max_freq"]
            protein_list = bin_group_tax_dict[key_bin]["protein_ids_k"]
            species_list = bin_group_tax_dict[key_bin]["species"]
            binned_taxa_file.write("bin:\n")
            binned_taxa_file.write("{}\n".format(key_bin))
            binned_taxa_file.write("max frequency:\n")
            binned_taxa_file.write("{}\n".format(max_freq))
            binned_taxa_file.write("binned proteins:\n")
            for item in protein_list:
                binned_taxa_file.write("{}\n".format(item))
            binned_taxa_file.write("species:\n")
            for item in species_list:
                binned_taxa_file.write("{}\n".format(item))
            binned_taxa_file.write("//\n")
        binned_taxa_file.close()
        # Write the bin and their protein IDs in a file.
        new_file_bin_proteins = open(binphylo_path_prsids, "w")
        for binid in dict_bins_prs_sorted.keys():
            for item in dict_bins_prs_sorted[binid]:
                new_file_bin_proteins.write("{}\t{}\n".format(binid, item))
        new_file_bin_proteins.close()
        # Write all the information realted to binning in one file.
        # Bin ID -> Protein IDs -> Contig IDs -> Read IDs
        # There is no need to check in the contig of contigs_to_reads_dict whether the read/contig belong to a bin or not. Each contig parsed from dict_bins_contigs is part of a bin, thus each
        # read of those contigs is also part of the bin.
        bin_sum_file = open(kraken_bin_info_path, "w")
        for key_bin in bin_group_tax_dict.keys():
            max_freq = bin_group_tax_dict[key_bin]["max_freq"]
            protein_list = bin_group_tax_dict[key_bin]["protein_ids_k"]
            species_list = bin_group_tax_dict[key_bin]["species"]
            bin_sum_file.write("Bin ID:\n")
            bin_sum_file.write("{}\n".format(key_bin))
            bin_sum_file.write("Max read frequency (to any of the binned contigs):\n")
            bin_sum_file.write("{}\n".format(max_freq))
            bin_sum_file.write("Species:\n")
            for sp_id in species_list:
                bin_sum_file.write("{}\n".format(sp_id))
            bin_sum_file.write("Binned protein IDs:\n")
            for pr_id in protein_list:
                bin_sum_file.write("{}\n".format(pr_id))
            bin_sum_file.write("Binned contig IDs:\n")
            for cg_id in dict_bins_contigs[key_bin]:
                bin_sum_file.write("{}\n".format(cg_id))
            bin_sum_file.write("Binned read IDs:\n")
            for cg_id in dict_bins_contigs[key_bin]:
                for rd_id in contigs_to_reads_dict[cg_id][1:]:
                    bin_sum_file.write("{}\n".format(rd_id))
            bin_sum_file.write("//\n")
        bin_sum_file.close()
    return bin_group_tax_dict