import os
import csv
# ProteoSeeker modules
import supportive_functions


def info_collection(comb_all_domains_proteins, proteins_combined_file, blastp_info_swissprot_file, blastp_info_comb_nr_file, topology_info_path, gene_info_file, gene_contig_dist_path, input_motifs_results, family_info_path, binning_results_path, output_path_genepred, output_path_genepred_faa, binphylo_path_prstax, binphylo_path_binstax, binphylo_path_binstax_max, binphylo_path_prsids, dict_hm, dict_sp, dict_nr, dict_top, swiss_fams_len_comp_dict, dict_input_motifs, dict_genes, contig_gene_dist_dict, dict_contigs_bins, dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict, taxonomy_route, binned_taxa_path, dict_seqs, seek_mode, taxonomy_mode, cd_hit_results_fasta_path, add_seek_info, add_taxonomy_info):
    print("\nCollecting the information associated with the annotation of the proteins...")
    # Binned proteins
    if not dict_bins_prs_sorted:
        if os.path.exists(binphylo_path_prsids):
            bin_phylo_bins_lines = supportive_functions.read_file(binphylo_path_prsids)
            for line in bin_phylo_bins_lines:
                line_splited = line.split("\t")
                bin_id = line_splited[0]
                protein_id = line_splited[1]
                if bin_id not in dict_bins_prs_sorted.keys():
                    dict_bins_prs_sorted[bin_id] = [protein_id]
                else:
                    dict_bins_prs_sorted[bin_id].append(protein_id)
    # Find the binned proteins.
    binned_protein_ids = set()
    if dict_bins_prs_sorted:
        for key_bin in dict_bins_prs_sorted.keys():
            binned_prs = dict_bins_prs_sorted[key_bin]
            for pr_item in binned_prs:
                binned_protein_ids.add(pr_item)
    # Sequences
    # 1. Seek mode: proteins_combined_file
    # 2. Taxonomy mode: Kraken2: Proteins from the bins.
    # 3. Taxonomy mode: MetaBinner / COMEBin: Proteins from the bins.
    dict_seqs = {}
    if (not taxonomy_mode) or seek_mode:
        if add_seek_info:
            if os.path.exists(proteins_combined_file):
                with open(proteins_combined_file) as lines_seqs_cf:
                    for line_pr_seq in lines_seqs_cf:
                        line_pr_seq = line_pr_seq.rstrip("\n")
                        if line_pr_seq[0] == ">":
                            protein_acc = line_pr_seq[1:]
                            dict_seqs[protein_acc] = ""
                        else:
                            dict_seqs[protein_acc] += line_pr_seq
    if (not seek_mode) or taxonomy_mode:
        if add_taxonomy_info:
            found_binned_pr = False
            if os.path.exists(cd_hit_results_fasta_path):
                with open(cd_hit_results_fasta_path) as lines_seqs_ch:
                    for line_pr_seq in lines_seqs_ch:
                        line_pr_seq = line_pr_seq.rstrip("\n")
                        if line_pr_seq[0] == ">":
                            found_binned_pr = False
                            protein_acc = line_pr_seq[1:]
                            if protein_acc in binned_protein_ids:
                                found_binned_pr = True
                            if found_binned_pr:
                                dict_seqs[protein_acc] = ""
                        else:
                            if found_binned_pr:
                                dict_seqs[protein_acc] += line_pr_seq
    annotated_proteins_num = len(list(dict_seqs.keys()))
    print("\nNumber of proteins in the results: {}".format(annotated_proteins_num))
    # Domains
    if not dict_hm:
        if os.path.exists(comb_all_domains_proteins):
            dict_hm = {}
            header = True
            with open(comb_all_domains_proteins) as dom_lines:
                for line in dom_lines:
                    line = line.rstrip("\n")
                    line_splited = line.split("\t")
                    pr_name = line_splited[0]
                    if not header:
                        if pr_name not in dict_hm.keys():
                            dict_hm[pr_name] = [line]
                        else:
                            dict_hm[pr_name].append(line)
                    header = False
    # Phobius
    if not dict_top:
        if os.path.exists(topology_info_path):
            dict_top = {}
            with open(topology_info_path) as phobius_lines:
                for line in phobius_lines:
                    line = line.rstrip("\n")
                    # k141_12880_94839_96755_+\t1\t41\tSIGNAL
                    # k141_12880_94839_96755_+\t906\t926\tDOMAIN||C-REGION
                    # k141_12880_94839_96755_+\t906\t926\tDOMAIN||NON CYTOPLASMIC
                    # k127_2204_1502_2944_-	250	268	TRANSMEM
                    if "\t" in line:
                        line_splited = line.split("\t")
                        pr_id = line_splited[0]
                        temp_phrase = "\t".join(line_splited[1:])
                        if pr_id not in dict_top.keys():
                            dict_top[pr_id] = []
                        dict_top[pr_id].append(temp_phrase)
    # Swiss-Prot. Keeping only the protein with the lowest e-value from the hits against the Swiss-Prot database.
    if not dict_sp:
        if os.path.exists(blastp_info_swissprot_file):
            dict_sp = {}
            protein_acc = None
            min_evalue = None
            header = True
            with open(blastp_info_swissprot_file) as lines_sp:
                for line_sp in lines_sp:
                    if header:
                        header = False
                        continue
                    line_sp = line_sp.rstrip("\n")
                    splited_sp = line_sp.split("\t")
                    protein_acc = splited_sp[0]
                    e_value = float(splited_sp[3])
                    if protein_acc not in dict_sp.keys():
                        dict_sp[protein_acc] = line_sp
                        min_evalue = e_value
                    else:
                        if e_value < min_evalue:
                            dict_sp[protein_acc] = line_sp
                            min_evalue = e_value
    # NR. Keeping only the protein with the lowest e-calue from the hits against the nr database.
    if not dict_nr:
        if os.path.exists(blastp_info_comb_nr_file):
            dict_nr = {}
            protein_acc = None
            min_evalue = None
            header = True
            with open(blastp_info_comb_nr_file) as lines_nr:
                for line_nr in lines_nr:
                    if header:
                        header = False
                        continue
                    line_nr = line_nr.rstrip("\n")
                    splited_nr = line_nr.split("\t")
                    protein_acc = splited_nr[0]
                    e_value = float(splited_nr[3])
                    if protein_acc not in dict_nr.keys():
                        dict_nr[protein_acc] = line_nr
                        min_evalue = e_value
                    else:
                        if e_value < min_evalue:
                            dict_nr[protein_acc] = line_nr
                            min_evalue = e_value
    # Genes
    if not dict_genes:
        if os.path.exists(gene_info_file):
            dict_genes = {}
            gene_lines = supportive_functions.read_file(gene_info_file)
            for line in gene_lines:
                line_splited = line.split("\t")
                dict_genes[line_splited[0]] = [line_splited[1], line_splited[2], line_splited[3]]
    # Gene distances from contig's edges.
    if not contig_gene_dist_dict:
        if os.path.exists(gene_contig_dist_path):
            contig_gene_dist_dict = {}
            gene_dist_lines = supportive_functions.read_file(gene_contig_dist_path)
            for line in gene_dist_lines:
                line_splited = line.split("\t")
                protein_name = line_splited[0]
                start_dist = line_splited[1]
                end_dist = line_splited[2]
                contig_gene_dist_dict[protein_name] = [start_dist, end_dist]
    # Input motifs
    if not dict_input_motifs:
        dict_input_motifs = {}
        if os.path.exists(input_motifs_results):
            if os.path.getsize(input_motifs_results) != 0:
                motifs_lines = supportive_functions.read_file(input_motifs_results)
                for line in motifs_lines:
                    line_splited = line.split("\t")
                    protein_acc = line_splited[0]
                    mof = line_splited[1]
                    result_start = line_splited[2]
                    result_end = line_splited[3]
                    temp_mof_list = [mof, result_start, result_end]
                    if protein_acc not in dict_input_motifs.keys():
                        dict_input_motifs[protein_acc] = [temp_mof_list]
                    else:
                        dict_input_motifs[protein_acc].append(temp_mof_list)
    # SwissProt predicted family and length difference.
    if not swiss_fams_len_comp_dict:
        if os.path.exists(family_info_path):
            swiss_fams_len_comp_dict = {}
            swiss_fam_lines = supportive_functions.read_file(family_info_path)
            for line in swiss_fam_lines:
                line_splited = line.split("\t")
                protein_name = line_splited[0]
                pf_name = line_splited[1]
                pf_mean_length = line_splited[2]
                len_diff = line_splited[3]
                len_rel_change = line_splited[4]
                fam_dif = line_splited[5]
                temp_list = [pf_name, pf_mean_length, len_diff, len_rel_change, fam_dif]
                if pf_name not in swiss_fams_len_comp_dict.keys():
                    swiss_fams_len_comp_dict[protein_name] = [temp_list]
                else:
                    swiss_fams_len_comp_dict[protein_name].append(temp_list)
    # Binning results
    if (not dict_contigs_bins) or (not dict_prs_bins):
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
    # Bin taxonomies
    if not bins_prs_tax_dict:
        if os.path.exists(binphylo_path_prstax):
            bin_phylo_prs_lines = supportive_functions.read_file(binphylo_path_prstax)
            for line in bin_phylo_prs_lines:
                line_splited = line.split("\t")
                bin_id = line_splited[0]
                protein_id = line_splited[1]
                species_id = line_splited[2]
                species_num = int(line_splited[3])
                if bin_id not in bins_prs_tax_dict.keys():
                    bins_prs_tax_dict[bin_id] = {}
                if protein_id not in bins_prs_tax_dict[bin_id].keys():
                    bins_prs_tax_dict[bin_id][protein_id] = {}
                if species_id not in bins_prs_tax_dict[bin_id][protein_id].keys():
                    bins_prs_tax_dict[bin_id][protein_id][species_id] = species_num
    if not bin_tax_dict:
        if os.path.exists(binphylo_path_binstax):
            bin_phylo_bins_lines = supportive_functions.read_file(binphylo_path_binstax)
            for line in bin_phylo_bins_lines:
                line_splited = line.split("\t")
                bin_id = line_splited[0]
                species_id = line_splited[1]
                species_num = int(line_splited[2])
                if bin_id not in bin_tax_dict.keys():
                    bin_tax_dict[bin_id] = {}
                if species_id not in bin_tax_dict[bin_id].keys():
                    bin_tax_dict[bin_id][species_id] = species_num
    if not bin_tax_max_dict:
        if os.path.exists(binphylo_path_binstax_max):
            bin_phylo_bins_max_lines = supportive_functions.read_file(binphylo_path_binstax_max)
            for line in bin_phylo_bins_max_lines:
                line_splited = line.split("\t")
                bin_id = line_splited[0]
                species_id = line_splited[1]
                species_num = int(line_splited[2])
                if bin_id not in bin_tax_max_dict.keys():
                    bin_tax_max_dict[bin_id] = {}
                if species_id not in bin_tax_max_dict[bin_id].keys():
                    bin_tax_max_dict[bin_id][species_id] = species_num
    if not bin_group_tax_dict:
        # As for the kraken analysis a similar file "binned_taxa_path" should be created for binnin with MetaBinner and COMEBin and ue it below.
        if taxonomy_route == 1:
            if os.path.exists(binned_taxa_path):
                pr_line = ""
                protein_status = False
                species_status = False
                with open(binned_taxa_path) as bintaxa_lines:
                    for line in bintaxa_lines:
                        line = line.rstrip("\n")
                        if pr_line == "bin:":
                            bin_id = line
                            bin_group_tax_dict[bin_id] = {
                                "max_freq": None,
                                "protein_ids_k": [],
                                "species": []
                            }
                        elif pr_line == "max frequency:":
                            bin_group_tax_dict[bin_id]["max_freq"] = line
                        elif pr_line == "binned proteins:":
                            protein_status = True
                            species_status = False
                        elif pr_line == "species:":
                            protein_status = False
                            species_status = True
                        elif line == "//":
                            protein_status = False
                            species_status = False
                        if protein_status:
                            bin_group_tax_dict[bin_id]["protein_ids_k"].append(line)
                        if species_status:
                            bin_group_tax_dict[bin_id]["species"].append(line)
                        pr_line = line
        else:
            if os.path.exists(binning_results_path):
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
    return dict_seqs, dict_hm, dict_top, dict_sp, dict_nr, dict_genes, contig_gene_dist_dict, dict_input_motifs, swiss_fams_len_comp_dict, dict_contigs_bins, dict_prs_bins, dict_bins_prs_sorted, bins_prs_tax_dict, bin_tax_dict, bin_tax_max_dict, bin_group_tax_dict