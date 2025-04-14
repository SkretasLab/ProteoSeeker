import csv


def tsv_results(annotation_file_tsv_name, dict_hm, dict_seqs, dict_top, dict_sp, dict_nr, dict_input_motifs, dict_genes, swiss_fams_len_comp_dict, contig_gene_dist_dict, pr_tax_info_dcit, add_type, add_info):
    print("\nWriting the results in a TSV File...")

    # If the file path to which the results are added is None then a new TSV file must be created. Otherwise, the results are saved to the already existing file.
    with open(annotation_file_tsv_name, 'w', newline='') as tsvf:
        tsv_writer = csv.writer(tsvf, delimiter='\t', lineterminator='\n')
        # First line headers.
        first_row_headers = ["Putative Protein ID", "Protein Sequence", "Sequence Length", "Bin ID", "Bin Taxonomies", "Domains", "Domains", "Domains", "Domains", "Domains", "Domains", "Domains", "Domains",
                    "Domains", "Domains", "Domains", "Input motifs", "Input motifs", "Input motifs", "Protein family", "Protein family", "Protein family", "Protein family", "Protein family",
                    "Signal Peptide or Transmembrane Region", "Signal Peptide or Transmembrane Region", "Signal Peptide or Transmembrane Region", "Swiss-Prot", "Swiss-Prot", "Swiss-Prot", "Swiss-Prot",
                    "FPD", "FPD", "FPD", "FPD", "Gene Information", "Gene Information", "Gene Information", "Gene Information", "Gene Information"]
        if add_type:
            print(len(add_type))
            temp_headers = ["Custom Information"] * len(add_type)
            first_row_headers = first_row_headers + temp_headers
        # Write row.
        tsv_writer.writerow(first_row_headers)

        # Second line headers.
        second_row_headers = ["Putative Protein ID", "Protein Sequence", "Sequence Length", "Bin ID", "Bin Taxonomies", "Description", "Name", "Pfam ID", "Sequence E-value", "Bitscore", 
                               "Domain c-E-value", "Domain i-E-value", "Domain from", "Domain to", "Sequence from", "Sequence to", "Input motif", "Region(s) from", "Region(s) to", "Name", "Mean length",
                               "Length difference", "Length relative change", "Family difference", "Region", "Region from", "Region to", "Swiss-Prot Accession Number", "Percentage Identity", "E-value",
                               "Bitscore", "Protein ID", "Percentage Identity", "E-value", "Bitscore", "Gene Sequence", "Start Codon" , "End Codon", "Start Distance", "End Distance"]
        if add_type:
            for add_i in range(0, len(add_type)):
                type_item = add_type[add_i]
                second_row_headers.append(type_item)
        # Write row.
        tsv_writer.writerow(second_row_headers)

        # Total protein number to be written in the TSV file..
        total_seqs = len(dict_seqs.keys())

        # If the results are going to be added in an TSV file then find the maximum row of that TSV file
        counter_sqs_1 = 0
        counter_sqs_2 = 0
        perc_sqs = total_seqs*0.1
        for key_pr_ac in dict_seqs:
            counter_sqs_1 += 1
            counter_sqs_2 += 1
            if counter_sqs_1 >= perc_sqs:
                print("{} / {}".format(counter_sqs_2, total_seqs))
                counter_sqs_1 = 0
            
            # List to hold the data for the current row.
            data_row = []

            # Protein sequence and length
            seq = dict_seqs[key_pr_ac]
            seq_length = len(seq)
            # Protein ID
            data_row.append(key_pr_ac)
            # Protein sequence
            data_row.append(seq)
            # Sequence length
            data_row.append(seq_length)

            # Bin information
            if pr_tax_info_dcit:
                if key_pr_ac in pr_tax_info_dcit.keys():
                    pr_bin_id = pr_tax_info_dcit[key_pr_ac][0]
                    pr_taxonomies = pr_tax_info_dcit[key_pr_ac][1]
                    pr_cat_type = pr_tax_info_dcit[key_pr_ac][2]
                    cell_phrase = "-"
                    data_row.append(pr_bin_id)
                    pr_taxonomies_wr = ", ".join(pr_taxonomies)
                    cell_phrase = "{}: {}".format(pr_cat_type, pr_taxonomies_wr)
                    data_row.append(cell_phrase)
                else:
                    cell_phrase = "-"
                    data_row.append(cell_phrase)
                    data_row.append(cell_phrase)
            else:
                cell_phrase = "-"
                data_row.append(cell_phrase)
                data_row.append(cell_phrase)

            # Protein domains
            if key_pr_ac in dict_hm.keys():
                domain_info_seq_evalue = ""
                domain_info_bitscore = ""
                domain_info_name = ""
                domain_info_pfam_id = ""
                domain_info_description = ""
                domain_info_domain_cevalue = ""
                domain_info_domain_ievalue = ""
                domain_info_domain_from = ""
                domain_info_domain_to = ""
                domain_info_seq_from = ""
                domain_info_seq_to = ""
                for pr_i in range(0, len(dict_hm[key_pr_ac])):
                    pr_domains = dict_hm[key_pr_ac][pr_i]
                    domains_splited = pr_domains.split("\t")
                    if pr_i == 0:
                        domain_info_seq_evalue = domains_splited[1]
                        domain_info_bitscore = domains_splited[2]
                        domain_info_name = domains_splited[3]
                        domain_info_pfam_id = domains_splited[4]
                        domain_info_description = domains_splited[5]
                        domain_info_domain_cevalue = domains_splited[6]
                        domain_info_domain_ievalue = domains_splited[7]
                        domain_info_domain_from = domains_splited[8]
                        domain_info_domain_to = domains_splited[9]
                        domain_info_seq_from = domains_splited[10]
                        domain_info_seq_to = domains_splited[11]
                    else:
                        domain_info_seq_evalue = "{}||{}".format(domain_info_seq_evalue, domains_splited[1])
                        domain_info_bitscore = "{}||{}".format(domain_info_bitscore, domains_splited[2])
                        domain_info_name = "{}||{}".format(domain_info_name, domains_splited[3])
                        domain_info_pfam_id = "{}||{}".format(domain_info_pfam_id, domains_splited[4])
                        domain_info_description = "{}||{}".format(domain_info_description, domains_splited[5])
                        domain_info_domain_cevalue = "{}||{}".format(domain_info_domain_cevalue, domains_splited[6])
                        domain_info_domain_ievalue = "{}||{}".format(domain_info_domain_ievalue, domains_splited[7])
                        domain_info_domain_from = "{}||{}".format(domain_info_domain_from, domains_splited[8])
                        domain_info_domain_to = "{}||{}".format(domain_info_domain_to, domains_splited[9])
                        domain_info_seq_from = "{}||{}".format(domain_info_seq_from, domains_splited[10])
                        domain_info_seq_to = "{}||{}".format(domain_info_seq_to, domains_splited[11])
                data_row.append(domain_info_description)
                data_row.append(domain_info_name)
                data_row.append(domain_info_pfam_id)
                data_row.append(domain_info_seq_evalue)
                data_row.append(domain_info_bitscore)
                data_row.append(domain_info_domain_cevalue)
                data_row.append(domain_info_domain_ievalue)
                data_row.append(domain_info_domain_from)
                data_row.append(domain_info_domain_to)
                data_row.append(domain_info_seq_from)
                data_row.append(domain_info_seq_to)
            else:
                temp_data = ["-"] * 11
                data_row = data_row + temp_data

            # Input motifs
            if key_pr_ac not in dict_input_motifs.keys():
                temp_data = ["-"] * 3
                data_row = data_row + temp_data
            else:
                motif_sqs = None
                motifs_start = None
                motifs_end = None
                for pr_i in range(0, len(dict_input_motifs[key_pr_ac])):
                    item = dict_input_motifs[key_pr_ac][pr_i]
                    if pr_i == 0:
                        motif_pre_sqs = item[0]
                        motif_pre_sqs_splited = motif_pre_sqs.split("\\S")
                        motif_sqs = "X".join(motif_pre_sqs_splited)
                        motifs_start = item[1]
                        motifs_end = item[2]
                    else:
                        motif_sqs_cur = item[0]
                        motif_sqs_cur_splited = motif_sqs_cur.split("\\S")
                        motif_sqs_cur = "X".join(motif_sqs_cur_splited)
                        motif_sqs = "{}||{}".format(motif_sqs, motif_sqs_cur)
                        motifs_start = "{}||{}".format(motifs_start, item[1])
                        motifs_end = "{}||{}".format(motifs_end, item[2])
                data_row.append(motif_sqs)
                data_row.append(motifs_start)
                data_row.append(motifs_end)
            
            # Protein family predictions.
            if (not swiss_fams_len_comp_dict) or (key_pr_ac not in swiss_fams_len_comp_dict.keys()):
                temp_data = ["-"] * 5
                data_row = data_row + temp_data
            else:
                pf_name = None
                pf_mean_length = None
                len_diff = None
                len_rel_change = None
                for pr_i in range(0, len(swiss_fams_len_comp_dict[key_pr_ac])):
                    item = swiss_fams_len_comp_dict[key_pr_ac][pr_i]
                    if pr_i == 0:
                        pf_name = item[0]
                        pf_mean_length = item[1]
                        len_diff = item[2]
                        len_rel_change = item[3]
                        fam_dif = item[4]
                    else:
                        pf_name_cur = item[0]
                        pf_mean_length_cur = item[1]
                        len_diff_cur = item[2]
                        len_rel_change_cur = item[3]
                        fam_dif_cur = item[4]
                        pf_name = "{}||{}".format(pf_name, pf_name_cur)
                        pf_mean_length = "{}||{}".format(pf_mean_length, pf_mean_length_cur)
                        len_diff = "{}||{}".format(len_diff, len_diff_cur)
                        len_rel_change = "{}||{}".format(len_rel_change, len_rel_change_cur)
                        fam_dif = "{}||{}".format(fam_dif, fam_dif_cur)
                data_row.append(pf_name)
                data_row.append(pf_mean_length)
                data_row.append(len_diff)
                data_row.append(len_rel_change)
                data_row.append(fam_dif)

            # Signal peptide and transmembrane regions
            if dict_top:
                if key_pr_ac not in dict_top.keys():
                    temp_data = ["-"] * 3
                    data_row = data_row + temp_data
                else:
                    region_type = None
                    region_from = None
                    region_to = None
                    for sp_i in range(0, len(dict_top[key_pr_ac])):
                        sp_tr_region = dict_top[key_pr_ac][sp_i]
                        sp_tr_region_splited = sp_tr_region.split("\t")
                        if sp_i == 0:
                            region_type = sp_tr_region_splited[2]
                            region_from = sp_tr_region_splited[0]
                            region_to = sp_tr_region_splited[1]
                        else:
                            region_type_cur = sp_tr_region_splited[2]
                            region_from_cur = sp_tr_region_splited[0]
                            region_to_cur = sp_tr_region_splited[1]
                            region_type = "{}||{}".format(region_type, region_type_cur)
                            region_from = "{}||{}".format(region_from, region_from_cur)
                            region_to = "{}||{}".format(region_to, region_to_cur)
                        data_row.append(region_type)
                        data_row.append(region_from)
                        data_row.append(region_to)
            else:
                temp_data = ["-"] * 3
                data_row = data_row + temp_data

            # Best hit in the SwissProt protein database
            if key_pr_ac in dict_sp.keys():
                sp_hit = dict_sp[key_pr_ac]
                sp_hit = sp_hit.split("\t")
                sp_name = sp_hit[1]
                sp_evalue = sp_hit[2]
                sp_identity = sp_hit[3]
                sp_bitscore = sp_hit[4]
                data_row.append(sp_name)
                data_row.append(sp_evalue)
                data_row.append(sp_identity)
                data_row.append(sp_bitscore)
            else:
                temp_data = ["-"] * 4
                data_row = data_row + temp_data

            # Best hit in the NR protein database
            if key_pr_ac in dict_nr.keys():
                nr_hit = dict_nr[key_pr_ac]
                nr_hit_splited = nr_hit.split("\t")
                nr_protein_id = nr_hit_splited[1]
                nr_evalue = nr_hit_splited[2]
                nr_identity = nr_hit_splited[3]
                nr_bitscore = nr_hit_splited[4]
                data_row.append(nr_protein_id)
                data_row.append(nr_evalue)
                data_row.append(nr_identity)
                data_row.append(nr_bitscore)
            else:
                temp_data = ["-"] * 4
                data_row = data_row + temp_data

            # Information for gene sequences.
            if key_pr_ac in dict_genes.keys():
                gene_sequence = dict_genes[key_pr_ac][0]
                start_codon = dict_genes[key_pr_ac][1]
                end_codon = dict_genes[key_pr_ac][2]
                data_row.append(gene_sequence)
                data_row.append(start_codon)
                data_row.append(end_codon)
            else:
                temp_data = ["-"] * 3
                data_row = data_row + temp_data

            # Information for gene distance from contig's edges.
            if (contig_gene_dist_dict is not None) and (key_pr_ac in contig_gene_dist_dict.keys()):
                start_dist = contig_gene_dist_dict[key_pr_ac][0]
                end_dist = contig_gene_dist_dict[key_pr_ac][1]
                data_row.append(start_dist)
                data_row.append(end_dist)
            else:
                temp_data = ["-"] * 2
                data_row = data_row + temp_data

            # Adding information set by the user
            if add_info:
                for add_i in range(0, len(add_info)):
                    type_info = add_info[add_i]
                    data_row.append(type_info)

            # Writting the row with the data in the TSV file.
            tsv_writer.writerow(data_row)

        # Final count of the proteins.
        print("{} / {}".format(total_seqs, total_seqs))
