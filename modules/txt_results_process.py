import os
import shutil


def txt_results(annotation_file_txt_name, output_path_annotation, dict_hm, dict_seqs, dict_top, dict_sp, dict_nr, dict_input_motifs, dict_genes, swiss_fams_len_comp_dict, contig_gene_dist_dict, pr_tax_info_dcit, add_type, add_info):
    print("\nWriting the results in a TXT File...")
    # Final folder with annotation files is created.
    if os.path.exists(output_path_annotation):
        shutil.rmtree(output_path_annotation)
    os.mkdir(output_path_annotation)
    # Creating the TXT file for the results
    annotation_file = open(annotation_file_txt_name, "a")
    counter_seq = 0
    split_size = 70
    percentage = 10
    total_seqs = len(dict_seqs.keys())
    ten_percent_total_seqs = total_seqs * 0.1
    for key_pr_ac in dict_seqs:
        counter_seq += 1
        if counter_seq >= ten_percent_total_seqs:
            print("{}%".format(percentage))
            percentage += 10
            counter_seq = 0
        # Protein sequence and length
        annotation_file.write("\\\\{}\n".format(100*"-"))
        annotation_file.write("Putative protein ID:\n{}\n\n".format(key_pr_ac))
        seq = dict_seqs[key_pr_ac]
        seq_length = len(seq)
        annotation_file.write("Protein sequence:\n")
        ann_fasta_splited = [seq[anni:anni + split_size] for anni in range(0, len(seq), split_size)]
        for ann_fasta_line in ann_fasta_splited:
            annotation_file.write("{}\n".format(ann_fasta_line))
        annotation_file.write("\nSequence length:")
        annotation_file.write("\n{}\n\n".format(seq_length))
        annotation_file.write("Domains (Description\tName\tPfam ID\tSequence E-value\tBitscore\tDomain c-E-value\tDomain i-E-value\tDomain from\tDomain to\tSequence from\tSequence to):")
        if key_pr_ac in dict_hm.keys():
            for val_list in dict_hm[key_pr_ac]:
                domains_splited = val_list.split("\t")
                domains_line = "\t".join([domains_splited[5], domains_splited[3], domains_splited[4], domains_splited[1], domains_splited[2], domains_splited[6], domains_splited[7], domains_splited[8], domains_splited[9], domains_splited[10], domains_splited[11]])
                annotation_file.write("\n{}".format(domains_line))
        else:
            annotation_file.write("\nNo domains found.")
        # Bin taxonomies
        # There is a chance that the proteins, which are encoded from genes whose taxonomy was directly inferred, are not present in the annotation files because they do not contain
        # seek domains.
        if pr_tax_info_dcit:
            if key_pr_ac in pr_tax_info_dcit.keys():
                pr_bin_id = pr_tax_info_dcit[key_pr_ac][0]
                pr_taxonomies = pr_tax_info_dcit[key_pr_ac][1]
                pr_cat_type = pr_tax_info_dcit[key_pr_ac][2]
                annotation_file.write("\n\nBin ID:\n{}".format(pr_bin_id))
                annotation_file.write("\n\nProtein taxonomy - {}:\n".format(pr_cat_type))
                for ti in range(0, len(pr_taxonomies)):
                    item = pr_taxonomies[ti]
                    if (ti + 1) == len(pr_taxonomies):
                        annotation_file.write("{}".format(item))
                    else:
                        annotation_file.write("{}\t".format(item))
            else:
                annotation_file.write("\n\nBin ID:\n-")
                annotation_file.write("\n\nProtein taxonomy:\nCould not be inferred.")
        # Signal peptide and transmembrane regions
        if dict_top:
            annotation_file.write("\n\nSignal peptide and Transmembrane regions (Region\tRegion from\tRegion to):")
            if key_pr_ac not in dict_top.keys():
                annotation_file.write("\nNo signal peptide or transmembrane regions found.")
            else:
                for sp_tr_region in dict_top[key_pr_ac]:
                    sp_tr_region_splited = sp_tr_region.split("\t")
                    annotation_file.write("\n{}\t{}\t{}".format(sp_tr_region_splited[2], sp_tr_region_splited[0],  sp_tr_region_splited[1]))
        # Best hit in the Swissprot protein database.
        if key_pr_ac in dict_sp.keys():
            annotation_file.write("\n\nSwissProt database hit (SwissProt Accession Number\tE-value\tPercentage Identity\tBitscore):")
            sp_hit = dict_sp[key_pr_ac]
            sp_hit = sp_hit.split("\t")
            sp_line = "\t".join([sp_hit[1], sp_hit[3], sp_hit[2], sp_hit[4]])
            annotation_file.write("\n{}".format(sp_line))
        # Best hit in the NR protein database
        if key_pr_ac in dict_nr.keys():
            annotation_file.write("\n\nNon-redudant database hit (NR protein ID\tE-value\tPercentage Identity\tBitscore):")
            nr_hit = dict_nr[key_pr_ac]
            nr_hit = nr_hit.split("\t")
            nr_line = "\t".join([nr_hit[1], nr_hit[3], nr_hit[2], nr_hit[4]])
            annotation_file.write("\n{}".format(nr_line))
        # Any input motif found for a protein is written
        if key_pr_ac in dict_input_motifs.keys():
            annotation_file.write("\n\nInput motifs identified (Motif\tSequence Start\tSequence End):")
            motif_hit = dict_input_motifs[key_pr_ac]
            for motif_list in motif_hit:
                motif_seq = motif_list[0]
                motif_seq_splited = motif_seq.split("\\S")
                motif_seq = "X".join(motif_seq_splited)
                annotation_file.write("\n{}\t{}\t{}".format(motif_seq, motif_list[1], motif_list[2]))
        # Gene information
        if key_pr_ac in dict_genes.keys():
            annotation_file.write("\n\nGene Sequence:")
            gene_sequence = dict_genes[key_pr_ac][0]
            annotation_file.write("\n>{}".format(key_pr_ac))
            ann_fasta_splited = [gene_sequence[anni:anni + split_size] for anni in range(0, len(gene_sequence), split_size)]
            for ann_fasta_line in ann_fasta_splited:
                annotation_file.write("\n{}".format(ann_fasta_line))
            star_codon_status = dict_genes[key_pr_ac][1]
            end_codon_status = dict_genes[key_pr_ac][2]
            annotation_file.write("\nStart Codon:\n{}".format(star_codon_status))
            annotation_file.write("\nEnd Codon:\n{}".format(end_codon_status))
        # Gene distance from contig's edges.
        if contig_gene_dist_dict and (key_pr_ac in contig_gene_dist_dict.keys()):
            annotation_file.write("\n\nGene distance from contigs (Distance from contig start\tDistance from contig end):")
            start_dist = contig_gene_dist_dict[key_pr_ac][0]
            end_dist = contig_gene_dist_dict[key_pr_ac][1]
            annotation_file.write("\n{}\t{}".format(start_dist, end_dist))
        # Family prediction and family mean length comparison.
        if swiss_fams_len_comp_dict and (key_pr_ac in swiss_fams_len_comp_dict.keys()):
            annotation_file.write("\n\nPredicted protein family (PPF) and length comparison (PPF name\tPPF mean length\tPPF median length\tLength difference\tLength relative change\tFamily difference:")
            for item in swiss_fams_len_comp_dict[key_pr_ac]:
                pf_name = item[0]
                pf_mean_length = item[1]
                len_diff = item[2]
                len_rel_change = item[3]
                fam_dif = item[4]
                annotation_file.write("\n{}\t{}\t{}\t{}\t{}".format(pf_name, pf_mean_length, len_diff, len_rel_change, fam_dif))
        # Added types and information
        if add_type:
            annotation_file.write("\n\nCustom information:")
            for add_i in range(0, len(add_type)):
                type_item = add_type[add_i]
                info_item = add_info[add_i]
                annotation_file.write("\n{}:\n{}".format(type_item, info_item))
        # Termination symbol(s) for current protein accession.
        annotation_file.write("\n//{}\n\n\n".format(100*"-"))
    print("100%")
    # The new Excel and TXT files are saved in memory.
    annotation_file.close()