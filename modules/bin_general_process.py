import re
import os
import csv
import sam_file_process
# ProteoSeeker modules
import supportive_functions


def binning_analysis(mapped_reads_path, binning_results_path, file_read_info_path, reads_to_contigs_dict, contigs_to_reads_dict, output_final_contigs_formated, dict_contigs_bins, contig_read_summary_path, binphylo_path_binstax_max, binphylo_maxfreq_taxids_path, maxfreq_lineage_path, bin_summary_info_path, bin_group_tax_dict, binner_bin_info_path):
    print("\nAnalyzing mapped reads and binned contigs...")
    if (not reads_to_contigs_dict) or (not contigs_to_reads_dict):
        # Find the aligned reads.
        reads_to_contigs_dict, contigs_to_reads_dict = sam_file_process.readmapal(mapped_reads_path)

    dict_bins_contigs = {}
    if os.path.exists(binning_results_path) and os.path.exists(file_read_info_path) and reads_to_contigs_dict and contigs_to_reads_dict:
        # Get the contig IDs of all the contigs.
        # Contig status: 1 for binned, 0 for not binned
        contig_status_dict = {}
        pattern_contig_id = re.compile(r'^>(\S+)')
        with open(output_final_contigs_formated) as conitg_lines:
            for line in conitg_lines:
                line = line.rstrip()
                result_contig_id = pattern_contig_id.search(line)
                if result_contig_id:
                    contig_id = result_contig_id.group(1)
                    contig_status_dict[contig_id] = 0
        # Contigs to bins
        if not dict_contigs_bins:
            with open(binning_results_path) as und_lines:
                tsv_lines = csv.reader(und_lines, delimiter="\t", quotechar='"')
                for line in tsv_lines:
                    contig_id = line[0]
                    bin_id = line[1]
                    if bin_id[:5] == "group":
                        bin_id = bin_id[5:]
                    contig_status_dict[contig_id] = 1
                    dict_contigs_bins[contig_id] = bin_id
                    bin_id_str = str(bin_id)
                    if bin_id_str not in dict_bins_contigs.keys():
                        dict_bins_contigs[bin_id_str] = []
                    dict_bins_contigs[bin_id_str].append(contig_id)
                    # The condition below exists because we have contigs that belong to bins and contigs to which reads were mapped. There is a chance, that a read was mapped to a contig
                    # that was not binned or that a contig was binned and to which contig no read was mapped.
                    # Therefore, there is a chance that contig_id from binning_results_path (results of binning) does not exist in the contigs mapped to reads.
                    if contig_id in contigs_to_reads_dict.keys():
                        contigs_to_reads_dict[contig_id][0] = 1
        else:
            for key_contig in dict_contigs_bins.keys():
                contig_status_dict[key_contig] = 1
                if key_contig in contigs_to_reads_dict.keys():
                    contigs_to_reads_dict[key_contig][0] = 1
                bin_id = dict_contigs_bins[key_contig]
                bin_id_str = str(bin_id)
                if bin_id_str not in dict_bins_contigs.keys():
                    dict_bins_contigs[bin_id_str] = []
                dict_bins_contigs[bin_id_str].append(key_contig)
        # Reads to bins
        # The read IDs in bins_to_reads are unique because they come from as keys from the reads_to_contigs_dict.
        bins_to_reads = {}
        for key_read in reads_to_contigs_dict.keys():
            contig_id = reads_to_contigs_dict[key_read][0]
            if contig_id in dict_contigs_bins.keys():
                bin_id = dict_contigs_bins[contig_id]
                if bin_id not in bins_to_reads.keys():
                    bins_to_reads[bin_id] = []
                bins_to_reads[bin_id].append(key_read)
        # Bins to read numbers
        bins_to_readnums = {}
        for key_bin in bins_to_reads.keys():
            read_number = len(bins_to_reads[key_bin])
            bins_to_readnums[key_bin] = read_number
        # Get the read numbers from the file.
        input_total_reads = None
        input_paired_reads = None
        trimmed_total_reads = None
        trimmed_paired_reads = None
        readinfo_lines = supportive_functions.read_file(file_read_info_path)
        line_index = 0
        for line in readinfo_lines:
            line_splited = line.split("\t")
            if line_index == 0:
                input_total_reads = int(line_splited[1])
            elif line_index == 1:
                input_paired_reads = int(line_splited[1])
            elif line_index == 2:
                trimmed_total_reads = int(line_splited[1])
            elif line_index == 3:
                trimmed_paired_reads = int(line_splited[1])
            line_index += 1
        # Comparing the bin numbers to the initial number of reads
        bin_to_readfreq_dict = {}
        for key_bin in bins_to_readnums.keys():
            bin_to_readfreq_dict[key_bin] = [None, None, None, None]
            bin_readnum = bins_to_readnums[key_bin]
            if input_paired_reads is not None:
                input_paired_readfreq = bin_readnum / input_paired_reads
                input_paired_readfreq = 100 * input_paired_readfreq
                input_paired_readfreq = round(input_paired_readfreq, 2)
                bin_to_readfreq_dict[key_bin][1] = input_paired_readfreq
            if trimmed_paired_reads is not None:
                trimmed_paired_readfreq = bin_readnum / trimmed_paired_reads
                trimmed_paired_readfreq = 100 * trimmed_paired_readfreq
                trimmed_paired_readfreq = round(trimmed_paired_readfreq, 2)
                bin_to_readfreq_dict[key_bin][3] = trimmed_paired_readfreq
        # Get the contigs that were not binned.
        binned_contigs_num = 0
        nonbinned_contigs_num = 0
        for contig_id in contig_status_dict.keys():
            if contig_status_dict[contig_id] == 1:
                binned_contigs_num += 1
            else:
                nonbinned_contigs_num += 1
        # Count the reads that belong to binned contigs.
        binned_reads_num = 0
        nonbinned_reads_num = 0
        for key_contig in contigs_to_reads_dict.keys():
            if contigs_to_reads_dict[key_contig][0] == 1:
                binned_reads_num += len(contigs_to_reads_dict[key_contig]) - 1
            else:
                nonbinned_reads_num += len(contigs_to_reads_dict[key_contig]) - 1
        # Input total reads:
        # 1) Reads filtered out by preprocessing.
        # 2) Reads not assembled to contigs.
        # 3) Reads assmbled to contigs.
        # Reads from contigs:
        # 1) Reads of binned contigs.
        # 2) Reads of non-binned contigs.
        # Reads from binned contigs:
        # 1) Reads from binned contigs with organism(s).
        # 2) Reads from binned contigs without organism(s).
        processed_paired_reads = input_paired_reads - trimmed_paired_reads
        assembled_paired_reads = len(list(reads_to_contigs_dict.keys()))
        non_assembled_paired_reads = trimmed_paired_reads - assembled_paired_reads
        print("\nTotal contigs: {}".format(len(list(contig_status_dict.keys()))))
        print("Binned contigs: {}".format(binned_contigs_num))
        print("Non-binned contigs: {}".format(nonbinned_contigs_num))
        print("\nTotal reads: {}".format(input_paired_reads))
        print("Filtered reads: {}".format(processed_paired_reads))
        print("Contig reads: {}".format(assembled_paired_reads))
        print("Non-contig reads: {}".format(non_assembled_paired_reads))
        print("Binned reads: {}".format(binned_reads_num))
        print("Non-binned reads, from the contig reads: {}".format(nonbinned_reads_num))
        # Store the general information for the contigs and the reads.
        contig_read_summary_file = open(contig_read_summary_path, "w")
        contig_read_summary_file.write("Type\tCount\n")
        contig_read_summary_file.write("Total contigs\t{}\n".format(len(list(contig_status_dict.keys()))))
        contig_read_summary_file.write("Binned contigs\t{}\n".format(binned_contigs_num))
        contig_read_summary_file.write("Non-binned contigs\t{}\n".format(nonbinned_contigs_num))
        contig_read_summary_file.write("Total reads\t{}\n".format(input_paired_reads))
        contig_read_summary_file.write("Filtered reads\t{}\n".format(processed_paired_reads))
        contig_read_summary_file.write("Contig reads\t{}\n".format(assembled_paired_reads))
        contig_read_summary_file.write("Non-contig reads\t{}\n".format(non_assembled_paired_reads))
        contig_read_summary_file.write("Binned reads\t{}\n".format(binned_reads_num))
        contig_read_summary_file.write("Non-binned reads, from the contig reads\t{}\n".format(nonbinned_reads_num))
        contig_read_summary_file.close()
        bins_to_summary_dict = {}
        if os.path.exists(binphylo_path_binstax_max) and os.path.exists(binphylo_maxfreq_taxids_path) and os.path.exists(maxfreq_lineage_path):
            # 1    Intestinimonas butyriciproducens    3
            lines_1 = supportive_functions.read_file(binphylo_path_binstax_max)
            for line in lines_1:
                line_splited = line.split("\t")
                bin = line_splited[0]
                taxname_1 = line_splited[1]
                if bin not in bins_to_summary_dict.keys():
                    bins_to_summary_dict[bin] = {}
                bins_to_summary_dict[bin][taxname_1] = []
            # Intestinimonas butyriciproducens    1297617    species
            max_taxids_lines = supportive_functions.read_file(binphylo_maxfreq_taxids_path)
            for line in max_taxids_lines:
                line_splited = line.split("\t")
                taxname_2 = line_splited[0]
                taxid_1 = line_splited[1]
                taxrank = line_splited[2]
                for key_bin in bins_to_summary_dict.keys():
                    if taxname_2 in bins_to_summary_dict[key_bin].keys():
                        bins_to_summary_dict[key_bin][taxname_2] = [taxid_1, taxrank]
            # 1297617    cellular organisms;Bacteria;Terrabacteria group;Bacillota;Clostridia;Eubacteriales;Eubacteriales incertae sedis;Intestinimonas;Intestinimonas butyriciproducens    species
            max_lineages_lines = supportive_functions.read_file(maxfreq_lineage_path)
            for line in max_lineages_lines:
                line_splited = line.split("\t")
                taxid_2 = line_splited[0]
                lineage = line_splited[1]
                for key_bin in bins_to_summary_dict.keys():
                    for name_bin in bins_to_summary_dict[key_bin].keys():
                        if bins_to_summary_dict[key_bin][name_bin][0] == taxid_2:
                            bins_to_summary_dict[key_bin][name_bin].append(lineage)
            for key_bin in bin_to_readfreq_dict.keys():
                # There is a chance that no species were found for a bin.
                if key_bin in bins_to_summary_dict.keys():
                    for taxname in bins_to_summary_dict[key_bin].keys():
                        readfreq_1 = bin_to_readfreq_dict[key_bin][1]
                        readfreq_3 = bin_to_readfreq_dict[key_bin][3]
                        bins_to_summary_dict[key_bin][taxname].append(readfreq_1)
                        bins_to_summary_dict[key_bin][taxname].append(readfreq_3)
        # Write information.
        bin_summary_info_file = open(bin_summary_info_path, "w")
        for key_bin in bins_to_summary_dict.keys():
            for key_name in bins_to_summary_dict[key_bin].keys():
                bin_summary_info_file.write("{}\t{}".format(key_bin, key_name))
                for item in bins_to_summary_dict[key_bin][key_name]:
                    bin_summary_info_file.write("\t{}".format(item))
                bin_summary_info_file.write("\n")
        bin_summary_info_file.close()
        # Write all the information realted to binning in one file.
        # Bin ID -> Protein IDs -> Contig IDs -> Read IDs
        # There is no need to check in the contig of contigs_to_reads_dict whether the read/contig belong to a bin or not. Each contig parsed from dict_bins_contigs is part of a bin, thus each
        # read of those contigs is also part of the bin.
        bin_sum_file = open(binner_bin_info_path, "w")
        for key_bin in bin_group_tax_dict.keys():
            max_freq = bin_group_tax_dict[key_bin]["max_freq"]
            if "protein_ids_d" in bin_group_tax_dict[key_bin].keys():
                protein_list = bin_group_tax_dict[key_bin]["protein_ids_d"]
            elif "protein_ids_i" in bin_group_tax_dict[key_bin].keys():
                protein_list = bin_group_tax_dict[key_bin]["protein_ids_i"]
            else:
                protein_list = []
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
                if cg_id in contigs_to_reads_dict.keys():
                    for rd_id in contigs_to_reads_dict[cg_id][1:]:
                        bin_sum_file.write("{}\n".format(rd_id))
            bin_sum_file.write("//\n")
        bin_sum_file.close()
    else:
        print("\nOne of the files needed for the binning analysis is missing.")