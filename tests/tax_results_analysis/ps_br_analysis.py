import os
import sys
import copy
import math
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sklearn.metrics import r2_score
from matplotlib.font_manager import FontProperties


def check_i_value(i_value, option_str):
    if i_value in ["False", "false", "f", "0"]:
        i_value = False
    elif i_value in ["True", "true", "t", "1"]:
        i_value = True
    else:
        print("Wrong value given for the '{}' option. Exiting.".format(option_str))
        exit()
    return i_value


def read_file(filename):
    file_handle = open(filename, "r")
    lines_raw = file_handle.readlines()
    lines = []
    for line in lines_raw:
        line = line.rstrip("\n")
        lines.append(line)
    file_handle.close()
    return lines


def br_analysis(benchmark_path, benchmark_mod_path, benchmark_mod_spabu_path):
    br_info_dict = {}
    write_status = False
    ready_to_stop = False
    benchmark_mod_file = open(benchmark_mod_path, "w")
    file_one_lines = read_file(benchmark_path)
    for line in file_one_lines:
        if line[:7] == "#sample":
            write_status = True
            ready_to_stop = False
            sample_id = line[-2:]
            sample_id = int(sample_id)
            br_info_dict[sample_id] = {}
        elif line == "@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE":
            ready_to_stop = True
        elif (not line) or (line == ""):
            if ready_to_stop:
                write_status = False
                benchmark_mod_file.write("{}\n".format(line))
        if write_status:
            if (line[:2] != "@@") and ("\t" in line):
                line_splited = line.split("\t")
                taxid = line_splited[0]
                abundance = line_splited[4]
                if abundance != "0":
                    benchmark_mod_file.write("{}\n".format(line))
                    # Add the information for the species in the dictionary.
                    # The information is added if the abudance is not 0.
                    if taxid not in br_info_dict[sample_id].keys():
                        br_info_dict[sample_id][taxid] = abundance
                    else:
                        print("Error. Duplicate taxid found for the same sample. Exiting.")
                        exit()
            else:
                benchmark_mod_file.write("{}\n".format(line))
    benchmark_mod_file.close()
    # Write information in the file specific for the species and their abundancies.
    benchmark_mod_spabu_file = open(benchmark_mod_spabu_path, "w")
    if br_info_dict:
        for key_sm in br_info_dict.keys():
            benchmark_mod_spabu_file.write("#{}\n".format(key_sm))
            for key_sp in br_info_dict[key_sm].keys():
                abu = br_info_dict[key_sm][key_sp]
                benchmark_mod_spabu_file.write("{} - {}\n".format(key_sp, abu))
            benchmark_mod_spabu_file.write("//\n")
    benchmark_mod_spabu_file.close()
    # Return the dictionary.
    return br_info_dict


def crfiles(ps_dir, ps_results, analysis_results_dict, time_dir, time_stats_dir_path):
    # Create the directory for the results.
    if os.path.exists(ps_dir):
        shutil.rmtree(ps_dir)
    os.mkdir(ps_dir)
    # Execution time information
    if os.path.exists(time_dir):
        shutil.rmtree(time_dir)
    os.mkdir(time_dir)
    # Create the sub-directories for each of the methods.
    for key in analysis_results_dict.keys():
        if os.path.exists(analysis_results_dict[key]):
            shutil.rmtree(analysis_results_dict[key])
        os.mkdir(analysis_results_dict[key])
    dir_sample_names = os.listdir(ps_results)
    for dsn in dir_sample_names:
        if dsn != "dbs":
            dsn_spliited = dsn.split("_")
            sample_id = dsn_spliited[1]
            dir_sample_path = "{}/{}".format(ps_results, dsn)
            dir_method_names = os.listdir(dir_sample_path)
            # Problem in copying the files. More files are copied in the directories than it should.
            for dmn in dir_method_names:
                dir_method_path = "{}/{}".format(dir_sample_path, dmn)
                if "kraken" in dmn:
                    kraken_dir = "{}/kraken_results".format(dir_method_path)
                    kraken_names = os.listdir(kraken_dir)
                    for kn in kraken_names:
                        if kn not in ["binned_taxa.txt", "binphylo_prsids.tsv",  "bins_to_contigs.tsv",  "contigs_to_bins.tsv",  "kraken.sh",  "kraken_bin_summary.txt",  "kraken_reads.tsv",  "kraken_stde.txt",  "kraken_taxa_names.tsv",  "kraken_version.txt",  "report_info.tsv",  "taxon_results.tsv"]:
                            kraken_source_path = "{}/{}".format(kraken_dir, kn)
                            if "kraken_8" in dmn:
                                target_dir = analysis_results_dict["k8"]
                                method_db = "k8"
                            elif "kraken_16" in dmn:
                                target_dir = analysis_results_dict["k16"]
                                method_db = "k16"
                            elif "kraken_72" in dmn:
                                target_dir = analysis_results_dict["k72"]
                                method_db = "k72"
                            if kn == "time_analysis.tsv":
                                kraken_dest_path = "{}/mdb_{}_sample_{}_{}".format(time_stats_dir_path, method_db, sample_id, kn)
                            else:
                                kraken_dest_path = "{}/sample_{}_{}".format(target_dir, sample_id, kn)
                            shutil.copyfile(kraken_source_path, kraken_dest_path)
                elif "metabinner" in dmn:
                    metabinner_dir = "{}/metabinner_results".format(dir_method_path)
                    metabinner_source_path_1 = "{}/binphylo_freq_taxids.tsv".format(metabinner_dir)
                    metabinner_source_path_2 = "{}/binphylo_maxfreq_taxids.tsv".format(metabinner_dir)
                    metabinner_source_path_3 = "{}/b_summary_info_metabinner.tsv".format(metabinner_dir)
                    metabinner_source_path_4 = "{}/time_analysis.tsv".format(metabinner_dir)
                    if "metabinner_nr" in dmn:
                        target_dir = analysis_results_dict["mnr"]
                        method_db = "mnr"
                    metabinner_dest_path_1 = "{}/sample_{}_binphylo_freq_taxids.tsv".format(target_dir, sample_id)
                    metabinner_dest_path_2 = "{}/sample_{}_binphylo_maxfreq_taxids.tsv".format(target_dir, sample_id)
                    metabinner_dest_path_3 = "{}/sample_{}_b_summary_info_metabinner.tsv".format(target_dir, sample_id)
                    metabinner_dest_path_4 = "{}/mdb_{}_sample_{}_time_analysis.tsv".format(time_stats_dir_path, method_db, sample_id)
                    if os.path.exists(metabinner_source_path_1):
                        shutil.copyfile(metabinner_source_path_1, metabinner_dest_path_1)
                    if os.path.exists(metabinner_source_path_2):
                        shutil.copyfile(metabinner_source_path_2, metabinner_dest_path_2)
                    if os.path.exists(metabinner_source_path_3):
                        shutil.copyfile(metabinner_source_path_3, metabinner_dest_path_3)
                    if os.path.exists(metabinner_source_path_4):
                        shutil.copyfile(metabinner_source_path_4, metabinner_dest_path_4)
                elif "comebin" in dmn:
                    comebin_dir = "{}/comebin_results".format(dir_method_path)
                    comebin_source_path_1 = "{}/binphylo_freq_taxids.tsv".format(comebin_dir)
                    comebin_source_path_2 = "{}/binphylo_maxfreq_taxids.tsv".format(comebin_dir)
                    comebin_source_path_3 = "{}/b_summary_info_comebin.tsv".format(comebin_dir)
                    comebin_source_path_4 = "{}/time_analysis.tsv".format(comebin_dir)
                    if "comebin_nr" in dmn:
                        target_dir = analysis_results_dict["cnr"]
                        method_db = "cnr"
                    comebin_dest_path_1 = "{}/sample_{}_binphylo_freq_taxids.tsv".format(target_dir, sample_id)
                    comebin_dest_path_2 = "{}/sample_{}_binphylo_maxfreq_taxids.tsv".format(target_dir, sample_id)
                    comebin_dest_path_3 = "{}/sample_{}_b_summary_info_comebin.tsv".format(target_dir, sample_id)
                    comebin_dest_path_4 = "{}/mdb_{}_sample_{}_time_analysis.tsv".format(time_stats_dir_path, method_db, sample_id)
                    if os.path.exists(comebin_source_path_1):
                        shutil.copyfile(comebin_source_path_1, comebin_dest_path_1)
                    if os.path.exists(comebin_source_path_2):
                        shutil.copyfile(comebin_source_path_2, comebin_dest_path_2)
                    if os.path.exists(comebin_source_path_3):
                        shutil.copyfile(comebin_source_path_3, comebin_dest_path_3)
                    if os.path.exists(comebin_source_path_4):
                        shutil.copyfile(comebin_source_path_4, comebin_dest_path_4)


def collect_kraken_filters(kraken_ps_dir, filter_method_label, kraken_filters_dict):
    for i in range(1, 20):
        i_str = str(i)
        kraken_filters_dict[filter_method_label][i_str] = ["-", "-"]
    kraken_filenames = os.listdir(kraken_ps_dir)
    for kfn in kraken_filenames:
        if "kraken_filters" in kfn:
            kfn_splited = kfn.split("_")
            kfn_sample_id = kfn_splited[1]
            kfn_path = "{}/{}".format(kraken_ps_dir, kfn)
            kfn_lines = read_file(kfn_path)
            kraken_filter_value = None
            shannon_index_value = None
            for line in kfn_lines:
                if "Filter: -2" in line:
                    print("Reached the non-gut filter for \"-2\" selection. Exiting.")
                    exit()
                if "Shannon Index: " in line:
                    line_splited_shannon = line.split("Shannon Index: ")
                    shannon_index_value = line_splited_shannon[1]
                    shannon_index_value = float(shannon_index_value)
                    shannon_index_value = round(shannon_index_value, 2)
                if "Kraken percentage filter: " in line:
                    line_splited_kraken = line.split("Kraken percentage filter: ")
                    kraken_filter_value = line_splited_kraken[1]
                    break
            if (kraken_filter_value is None) or (shannon_index_value is None):
                print("Kraken filtering value not found. Exiting.")
                exit()
            kraken_filters_dict[filter_method_label][kfn_sample_id] = [shannon_index_value, kraken_filter_value]
    return kraken_filters_dict


def process_kraken_filter_files(ps_dir, stats_dir_path):
    kraken_8_ps_dir = "{}/kraken_8".format(ps_dir)
    kraken_16_ps_dir = "{}/kraken_16".format(ps_dir)
    kraken_72_ps_dir = "{}/kraken_72".format(ps_dir)
    kraken_filters_dict = {
        "8": {},
        "16": {},
        "72": {}
    }
    if os.path.exists(kraken_8_ps_dir):
        filter_method_label = "8"
        kraken_filters_dict = collect_kraken_filters(kraken_8_ps_dir, filter_method_label, kraken_filters_dict)
    if os.path.exists(kraken_16_ps_dir):
        filter_method_label = "16"
        kraken_filters_dict = collect_kraken_filters(kraken_16_ps_dir, filter_method_label, kraken_filters_dict)
    if os.path.exists(kraken_72_ps_dir):
        filter_method_label = "72"
        kraken_filters_dict = collect_kraken_filters(kraken_72_ps_dir, filter_method_label, kraken_filters_dict)
    # Storing the information in a file.
    kraken_filters_stats_path = "{}/kraken_ng_filtering_values.tsv".format(stats_dir_path)
    kraken_filters_stats_file = open(kraken_filters_stats_path, "w")
    kraken_filters_stats_file.write("database\tsample\tshannon index\tnon-gut filtering value (%)\n")
    for key_method in kraken_filters_dict.keys():
        for key_sample in kraken_filters_dict[key_method].keys():
            si_value = kraken_filters_dict[key_method][key_sample][0]
            ft_value = kraken_filters_dict[key_method][key_sample][1]
            kraken_filters_stats_file.write("{}\t{}\t{}\t{}\n".format(key_method, key_sample, si_value, ft_value))
    kraken_filters_stats_file.close()
    

def ps_kraken_analyze(ps_kraken_dir, filter_name, sp_dir_path, spec_label):
    kraken_info_dict = {}
    unclassified_dict = {}
    kraken_files = os.listdir(ps_kraken_dir)
    # Kraken2
    # taxid read_number percentage
    # 2892996   729281  13.71
    # 1964449	193404	3.63
    # 1812480	56036	1.05
    for kf in kraken_files:
        if filter_name in kf:
            perc_sum = 0
            kf_nosuffix = kf[:-4]
            kf_splited = kf_nosuffix.split("_")
            sample_id = int(kf_splited[1])
            kraken_info_dict[sample_id] = {}
            unclassified_dict[sample_id] = {}
            relpath_kf = "{}/{}".format(ps_kraken_dir, kf)
            kraken_lines = read_file(relpath_kf)
            for line in kraken_lines:
                line_splited = line.split("\t")
                taxid = line_splited[0]
                if (taxid != "Total") and (taxid != "Toal"):
                    percentage = line_splited[2]
                    kraken_info_dict[sample_id][taxid] = percentage
                    percentage_float = float(percentage)
                    perc_sum += percentage_float
            perc_sum = round(perc_sum, 2)
            perc_sum_unclassified = 100 - perc_sum
            perc_sum_unclassified = round(perc_sum_unclassified, 2)
            unclassified_dict[sample_id] = perc_sum_unclassified
    # Write the information in a file.
    for key_sample in kraken_info_dict.keys():
        sp_sample_method_path = "{}/sample_{}_kraken_{}.tsv".format(sp_dir_path, key_sample, spec_label)
        sp_dir_file = open(sp_sample_method_path, "w")
        sp_dir_file.write("sample\ttaxid\trelative abundance (%)\n")
        for key_item in kraken_info_dict[key_sample].keys():
            key_percentage = kraken_info_dict[key_sample][key_item]
            sp_dir_file.write("{}\t{}\t{}\n".format(key_sample, key_item, key_percentage))
        perc_sum_unclassified = unclassified_dict[key_sample]
        sp_dir_file.write("{}\tunclassified\t{}\n".format(key_sample, perc_sum_unclassified))
    sp_dir_file.close()
    return kraken_info_dict


def ps_comebin_analyze(ps_cmbn, filter_name, sp_dir_path, spec_label):
    # bind      species_name                taxid   taxonomy_rank   lineage                                                                                                                                                                                     percentage_input_reads      percentage_preprocessed_reads
    # 0	        Cupriavidus	                106589	genus	        cellular organisms;Bacteria;Pseudomonadota;Betaproteobacteria;Burkholderiales;Burkholderiaceae;Cupriavidus                                                                                  1.48	                    1.48
    # 1	        Micromonospora krabiensis	307121	species	        cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora;Micromonospora krabiensis                              	1.45	                    1.46
    # 1	        Micromonospora	            1873	genus	        cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora	                                                        1.45	                    1.46
    # 1	        Nonomuraea sp. TT08I-71	    2733864	species	        cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Streptosporangiales;Streptosporangiaceae;Nonomuraea;unclassified Nonomuraea;Nonomuraea sp. TT08I-71	        1.45	                    1.46
    # 1	        Micromonospora sp. WMMA2032	2039870	species	        cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora;unclassified Micromonospora;Micromonospora sp. WMMA2032	1.45	                    1.46
    cmbn_info_dict = {}
    cmbn_sole_info_dict = {}
    cmbn_perc_sum_dict = {}
    unclassified_dict = {}
    cmbn_files = os.listdir(ps_cmbn)
    for cf in cmbn_files:
        if filter_name in cf:
            cf_nosuffix = cf[:-4]
            cf_splited = cf_nosuffix.split("_")
            sample_id = int(cf_splited[1])
            cmbn_info_dict[sample_id] = {}
            cmbn_sole_info_dict[sample_id] = {}
            cmbn_perc_sum_dict[sample_id] = {}
            relpath_cf = "{}/{}".format(ps_cmbn, cf)
            cmbn_lines = read_file(relpath_cf)
            header_line = True
            for line in cmbn_lines:
                if header_line:
                    header_line = False
                    continue
                line_splited = line.split("\t")
                line_splited_length = len(line_splited)
                # The scientific name of the organism associated with a bin is automatically extracted by the protein headers of the protein database,
                # hence sometimes it may not contain an actual species name.
                if line_splited_length != 7:
                    continue
                taxid = line_splited[2]
                rank = line_splited[3]
                if rank == "species":
                    percentage = float(line_splited[6])
                    # If the taxon ID has already been found then add its percentage to the current one.
                    # A dictionary will hold only the species and their percentages, which species are assigned only
                    # as the sole species to a bin ID. If the same species is found assigned to another (new) bin:
                    # A. If it is the sole species for the new bin ID its percentage is increased by that its percentage for the new bin ID.
                    # Other taxonomy ranks predicted for the bin are not taken into consideration.
                    # B. If it is not the sole species, none of the species are collected from the new bin and its percentage is not added to its existing entry in the dictionary.
                    # This makes sure that all the predicted species do not have a sum of percentages above 100 % and L1 norm can be computed. The species might not add up to 100 %.
                    # The remaining percentage for the "unknown" species is added as an immediate difference to the L1 norm.
                    if taxid not in cmbn_info_dict.keys():
                        cmbn_info_dict[sample_id][taxid] = percentage
                    else:
                        cmbn_info_dict[sample_id][taxid] += percentage
            # At first the bin IDs are filtered based on whether they correspond to one species or more.
            found_bin_ids = set()
            bin_ids_with_one_species = []
            bin_ids_with_mult_species = set()
            for line in cmbn_lines:
                line_splited = line.split("\t")
                line_splited_length = len(line_splited)
                # The scientific name of the organism associated with a bin is automatically extracted by the protein headers of the protein database,
                # hence sometimes it may not contain an actual species name. Any line with not the expected number of information is skipped.
                if line_splited_length != 7:
                    continue
                # Find the bin ID
                bin_id = line_splited[0]
                if bin_id not in found_bin_ids:
                    found_bin_ids.add(bin_id)
                else:
                    bin_ids_with_mult_species.add(bin_id)
            for fnid in found_bin_ids:
                if fnid not in bin_ids_with_mult_species:
                    bin_ids_with_one_species.append(fnid)
            # Collecting information for the bins asigned to one species.
            perc_sum = 0
            for line in cmbn_lines:
                line_splited = line.split("\t")
                line_splited_length = len(line_splited)
                # The scientific name of the organism associated with a bin is automatically extracted by the protein headers of the protein database,
                # hence sometimes it may not contain an actual species name. Any line with not the expected number of information is skipped.
                if line_splited_length != 7:
                    continue
                bin_id = line_splited[0]
                taxid = line_splited[2]
                rank = line_splited[3]
                # If the bin ID is one of the bin IDs for bins win one species then proceed.
                if bin_id in bin_ids_with_one_species:
                    if rank == "species":
                        sole_percentage = line_splited[6]
                        sole_percentage_fl = float(sole_percentage)
                        if taxid not in cmbn_sole_info_dict[sample_id].keys():
                            cmbn_sole_info_dict[sample_id][taxid] = sole_percentage_fl
                            # Compute the sum of the percentages of the identified species.
                            perc_sum += sole_percentage_fl
                        else:
                            cmbn_sole_info_dict[sample_id][taxid] += sole_percentage_fl
                            perc_sum += sole_percentage_fl
            perc_sum = round(perc_sum, 2)
            perc_sum_unclassified = 100 - perc_sum
            perc_sum_unclassified = round(perc_sum_unclassified, 2)
            cmbn_perc_sum_dict[sample_id] = [perc_sum, perc_sum_unclassified]
            unclassified_dict[sample_id] = perc_sum_unclassified
    # Write the information in a file.
    for key_sample in cmbn_sole_info_dict.keys():
        sp_sample_method_path = "{}/sample_{}_comebin_{}.tsv".format(sp_dir_path, key_sample, spec_label)
        sp_dir_file = open(sp_sample_method_path, "w")
        sp_dir_file.write("sample\ttaxid\trelative abundance (%)\n")
        for key_item in cmbn_sole_info_dict[key_sample].keys():
            key_percentage = cmbn_sole_info_dict[key_sample][key_item]
            sp_dir_file.write("{}\t{}\t{}\n".format(key_sample, key_item, key_percentage))
        perc_sum_unclassified = unclassified_dict[key_sample]
        sp_dir_file.write("{}\tunclassified\t{}\n".format(key_sample, perc_sum_unclassified))
    sp_dir_file.close()
    return cmbn_info_dict, cmbn_sole_info_dict, cmbn_perc_sum_dict


def ps_metabinner_analyze(ps_mtbr, filter_name, sp_dir_path, spec_label):
    # MetaBinner
    # bind      species_name    taxid   taxonomy_rank   lineage     percentage_input_reads      percentage_preprocessed_reads
    # 1	Listeria monocytogenes	1639	species	cellular organisms;Bacteria;Terrabacteria group;Bacillota;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria monocytogenes	9.6	10.25
    # 2	Pseudomonas aeruginosa	287	species	cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas aeruginosa group;Pseudomonas aeruginosa	14.76	15.76
    mtbr_info_dict = {}
    mtbr_sole_info_dict = {}
    mtbr_perc_sum_dict = {}
    unclassified_dict = {}
    mtbr_files = os.listdir(ps_mtbr)
    for mf in mtbr_files:
        if filter_name in mf:
            mf_nosuffix = mf[:-4]
            mf_splited = mf_nosuffix.split("_")
            sample_id = int(mf_splited[1])
            mtbr_info_dict[sample_id] = {}
            mtbr_sole_info_dict[sample_id] = {}
            mtbr_perc_sum_dict[sample_id] = {}
            relpath_mf = "{}/{}".format(ps_mtbr, mf)
            mtbr_lines = read_file(relpath_mf)
            for line in mtbr_lines:
                line_splited = line.split("\t")
                taxid = line_splited[2]
                rank = line_splited[3]
                if rank == "species":
                    percentage = line_splited[5]
                    # If the taxon ID has already been found then add its percentage to the current one.
                    # A dictionary will hold only the species and their percentages, which species are assigned only
                    # as the sole species to a bin ID. If the same species is found assigned to another (new) bin:
                    # A. If it is the sole species for the new bin ID its percentage is increased by that its percentage for the new bin ID.
                    # Other taxonomy ranks predicted for the bin are not taken into consideration.
                    # B. If it is not the sole species, none of the species are collected from the new bin and its percentage is not added to its existing entry in the dictionary.
                    # This makes sure that all the predicted species do not have a sum of percentages above 100 % and L1 norm can be computed. The species might not add up to 100 %.
                    # The remaining percentage for the "unknown" species is added as an immediate difference to the L1 norm.
                    if taxid not in mtbr_info_dict.keys():
                        mtbr_info_dict[sample_id][taxid] = percentage
                    else:
                        mtbr_info_dict[sample_id][taxid] += percentage
            # At first the bin IDs are filtered based on whether they correspond to one species or more.
            found_bin_ids = set()
            bin_ids_with_one_species = []
            bin_ids_with_mult_species = set()
            for line in mtbr_lines:
                line_splited = line.split("\t")
                line_splited_length = len(line_splited)
                # The scientific name of the organism associated with a bin is automatically extracted by the protein headers of the protein database,
                # hence sometimes it may not contain an actual species name. Any line with not the expected number of information is skipped.
                if line_splited_length != 7:
                    continue
                # Find the bin ID
                bin_id = line_splited[0]
                if bin_id not in found_bin_ids:
                    found_bin_ids.add(bin_id)
                else:
                    bin_ids_with_mult_species.add(bin_id)
            for fnid in found_bin_ids:
                if fnid not in bin_ids_with_mult_species:
                    bin_ids_with_one_species.append(fnid)
            # Collecting information for the bins asigned to one species.
            perc_sum = 0
            for line in mtbr_lines:
                line_splited = line.split("\t")
                line_splited_length = len(line_splited)
                # The scientific name of the organism associated with a bin is automatically extracted by the protein headers of the protein database,
                # hence sometimes it may not contain an actual species name. Any line with not the expected number of information is skipped.
                if line_splited_length != 7:
                    continue
                bin_id = line_splited[0]
                taxid = line_splited[2]
                rank = line_splited[3]
                # If the bin ID is one of the bin IDs for bins win one species then proceed.
                if bin_id in bin_ids_with_one_species:
                    if rank == "species":
                        sole_percentage = line_splited[6]
                        sole_percentage_fl = float(sole_percentage)
                        if taxid not in mtbr_sole_info_dict[sample_id].keys():
                            mtbr_sole_info_dict[sample_id][taxid] = sole_percentage_fl
                            # Compute the sum of the percentages of the identified species.
                            perc_sum += sole_percentage_fl
                        else:
                            mtbr_sole_info_dict[sample_id][taxid] += sole_percentage_fl
                            perc_sum += sole_percentage_fl
            perc_sum = round(perc_sum, 2)
            perc_sum_unclassified = 100 - perc_sum
            perc_sum_unclassified = round(perc_sum_unclassified, 2)
            mtbr_perc_sum_dict[sample_id] = [perc_sum, perc_sum_unclassified]
            unclassified_dict[sample_id] = perc_sum_unclassified
    # Write the information in a file.
    for key_sample in mtbr_sole_info_dict.keys():
        sp_sample_method_path = "{}/sample_{}_metabinner_{}.tsv".format(sp_dir_path, key_sample, spec_label)
        sp_dir_file = open(sp_sample_method_path, "w")
        sp_dir_file.write("sample\ttaxid\trelative abundance (%)\n")
        for key_item in mtbr_sole_info_dict[key_sample].keys():
            key_percentage = mtbr_sole_info_dict[key_sample][key_item]
            sp_dir_file.write("{}\t{}\t{}\n".format(key_sample, key_item, key_percentage))
        perc_sum_unclassified = unclassified_dict[key_sample]
        sp_dir_file.write("{}\tunclassified\t{}\n".format(key_sample, perc_sum_unclassified))
    sp_dir_file.close()
    return mtbr_info_dict, mtbr_sole_info_dict, mtbr_perc_sum_dict


def basic_stats(br_info_dict, pred_info_dict, pred_sole_info_dict, target_stats_dict, sample_id):
    # Benchmark | Sample | Result
    #    Yes    |   Yes  |   TP
    #    No     |   No   |   TN
    #    Yes    |   No   |   FN
    #    No     |   Yes  |   FP
    # True positive = Common species
    # False positive = Unique to prediction - Common species
    # False negative = Unique to validation - Common species
    # Sensitivity / Recall / True positive rate = TP / (TP + FN)
    # Specificity = TN / (TN + FP) = 0
    # Precision = TP / (TP + FP)
    # F1 Score = Sorensen-Dice Index = 2 * (Precision * Sensitivity) / (Precision + Sensitivity) = 2 * TP / (2 * TP + FP + FN)
    # L1 Norm = Sum of the absolute difference of the percentages of each species found in both of the groups, without parsing the same species more than once.
    # Accuracy = (TP + TN) / (TP + TN + FP + FN) = TP / (TP + FN)
    # Jaccard index = common_species / unique_species = TP / (TP + FP + FN)
    # Note: The results of MetaBinner and COMEBin may have more than one species for one bin.
    # All species for each bin are considered in the results. Also all species of the same bin
    tp = 0
    tn = 0
    fn = 0
    fp = 0
    sensitivity = 0
    precision = 0
    accuracy = 0
    f1_score = 0
    jaccard_index = 0
    l1_norm = 0

    # Number of decimals for each metric
    round_dec_num = 2

    # Gold standard sample information
    gold_si_info = br_info_dict[sample_id]

    # Determine whether the information from the predictions, which are to be compared with the standard samples,
    # originate from the predictions based on all species found or based on the species from the bins which bins
    # are each assigned one species. When the predictions come from the Kraken2 methods the "pred_sole_info_dict"
    # is None as there are no predictions to filter for bins with sole species because each bin is assigned always
    # one species.
    if pred_sole_info_dict is None:
        target_info_dict = copy.deepcopy(pred_info_dict)
    else:
        target_info_dict = copy.deepcopy(pred_sole_info_dict)
    target_si_info = target_info_dict[sample_id]
    
    # Remove the unclassified group and retain it for usage in computing the statistics.
    if "unclassified" in target_si_info.keys():
        del target_si_info['unclassified']

    # Create a group of all the items from gold and predicted items.
    # unique_items: Union of the species from the gold standard sample and the predicted group of species
    unique_items = set()
    for item in gold_si_info:
        unique_items.add(item)
    for item in target_si_info:
        unique_items.add(item)

    # Parse the list of unique items and check the existance of each item in each of the lists.
    # common_items: Common species
    # group_val_items: Gold species
    # group_val_unique_items: Unique gold species
    # group_pred_items: Predicted species
    # group_pred_unique_items: Unique predicted species
    # Four conditions:
    # In both groups.
    # In the predicted group and not in the validation group.
    # In the validation group and not in the predicted group.
    # In none of the groups.
    common_items = set()
    group_val_items = set()
    group_val_unique_items = set()
    group_pred_items = set()
    group_pred_unique_items = set()
    for item in unique_items:
        if (item in gold_si_info) and (item in target_si_info):
            common_items.add(item)
            group_val_items.add(item)
            group_pred_items.add(item)
        if (item not in gold_si_info) and (item in target_si_info):
            group_pred_items.add(item)
            group_pred_unique_items.add(item)
        if (item in gold_si_info) and (item not in target_si_info):
            group_val_items.add(item)
            group_val_unique_items.add(item)
        if (item not in gold_si_info) and (item not in target_si_info):
            tn += 1
    # Check for true negatives.
    if tn != 0:
        print("True negatives found. Exiting.")
        exit()

    # Compute TP, FP, FN
    tp = len(common_items)
    fp = len(group_pred_unique_items)
    fn = len(group_val_unique_items)

    # Sensitivity
    if (tp + fn) != 0:
        sensitivity =  tp / (tp + fn)
        sensitivity = sensitivity * 100
        sensitivity = round(sensitivity, round_dec_num)

    # Precision
    if (tp + fp) != 0:
        precision = tp / (tp + fp)
        precision = precision * 100
        precision = round(precision, round_dec_num)

    # Accuracy
    if (tp + tn + fp + fn) != 0:
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        accuracy = accuracy * 100
        accuracy = round(accuracy, round_dec_num)

    # F1 Score
    if (2 * tp + fp + fn) != 0:
        f1_score = (2 * tp) / (2 * tp + fp + fn)
        f1_score = round(f1_score, round_dec_num)

    # Jaccard index:
    # The true positive hits are the number of the items of the intersection of the two gruops.
    # The number of unique items of the two groups is the number of the items of their union.
    jaccard_index = tp / (tp + fp + fn)
    jaccard_index = round(jaccard_index, round_dec_num)

    # L1 norm
    for item in group_val_items:
        if item in br_info_dict[sample_id]:
            gold_abu_prec = float(br_info_dict[sample_id][item])
        else:
            print("Item expected to be in the validated group not found. Exiting.")
            exit()
        if item in target_info_dict[sample_id]:
            pred_abu_perc = float(target_info_dict[sample_id][item])
        else:
            pred_abu_perc = 0
        temp_term = abs(gold_abu_prec - pred_abu_perc)
        l1_norm += temp_term
    # Divide by 100 and round the L1 Norm.
    l1_norm = l1_norm / 100
    l1_norm = round(l1_norm, round_dec_num)

    # Storing information in a dictionary.
    val_species_num = len(gold_si_info)
    pred_species_num = len(target_si_info)
    unique_spec_num = len(unique_items)
    unique_val_num = len(group_val_unique_items)
    unique_pred_num = len(group_pred_unique_items)
    target_stats_dict[sample_id] = {
        "gold_species_number": val_species_num,
        "predicted_species_number": pred_species_num,
        "common_species": tp,
        "unique_species_total": unique_spec_num,
        "unique_gold": unique_val_num,
        "unique_predicted": unique_pred_num,
        "true_positive": tp,
        "false_positive": fp,
        "false_negative": fn,
        "sensitivity": sensitivity,
        "precision": precision,
        "accuracy": accuracy,
        "f1_score": f1_score,
        "jaccard_index": jaccard_index,
        "l1_norm": l1_norm
    }
    # Return value.
    return target_stats_dict, common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items


def write_list(path_tow, list_tow):
    file_tow = open(path_tow, "w")
    for item in list_tow:
        file_tow.write("{}\n".format(item))
    file_tow.close()


def write_species(common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items, stats_dir_path, label, si):
    # File path.
    sample_synopsis_path = "{}/synopsis_{}_sample_{}.tsv".format(stats_dir_path, label, si)
    # Converting to lists.
    common_items = list(common_items)
    group_val_items = list(group_val_items)
    group_val_unique_items = list(group_val_unique_items)
    group_pred_items = list(group_pred_items)
    group_pred_unique_items = list(group_pred_unique_items)
    # Lengths.
    common_items_len = len(common_items)
    group_val_items_len = len(group_val_items)
    group_val_unique_items_len = len(group_val_unique_items)
    group_pred_items_len = len(group_pred_items)
    group_pred_unique_items_len = len(group_pred_unique_items)
    groups_lengths_list = [common_items_len, group_val_items_len, group_val_unique_items_len, group_pred_items_len, group_pred_unique_items_len]
    # Max length.
    group_max_len = max(groups_lengths_list)
    # Fill values.
    fill_value = "-"
    if common_items_len < group_max_len:
        fill_list_ci = [fill_value] * (group_max_len - common_items_len)
        common_items += fill_list_ci
    if group_val_items_len < group_max_len:
        fill_list_gvi = [fill_value] * (group_max_len - group_val_items_len)
        group_val_items += fill_list_gvi
    if group_val_unique_items_len < group_max_len:
        fill_list_gvu = [fill_value] * (group_max_len - group_val_unique_items_len)
        group_val_unique_items += fill_list_gvu
    if group_pred_items_len < group_max_len:
        fill_list_gpi = [fill_value] * (group_max_len - group_pred_items_len)
        group_pred_items += fill_list_gpi
    if group_pred_unique_items_len < group_max_len:
        fill_list_gpu = [fill_value] * (group_max_len - group_pred_unique_items_len)
        group_pred_unique_items += fill_list_gpu
    # Dictionary with the species. It is needed to create the dataframe.
    species_dict = {
        "common_species": common_items,
        "validation_species": group_val_items,
        "validation_species_unique": group_val_unique_items,
        "predicted_species": group_pred_items,
        "predicted_species_unique": group_pred_unique_items
    }
    # Create the pandas dataframe for the synopsis data. Each empty cell is filled with the NaN value.
    syn_df = pd.DataFrame.from_dict(species_dict)
    # CSV file for one sample and one analysis with all the stats for the species
    syn_df.fillna("-").to_csv(sample_synopsis_path, sep="\t", index=False)


def comp_stats(br_info_dict, comp_dict, comp_sole_dict, stats_dir_path, label):
    print("Comparing results for {} with the validation data...".format(label))
    # have the same abundance.
    sample_ids = []
    sample_ids_dupl = list(br_info_dict.keys())
    for item in sample_ids_dupl:
        if item in sample_ids:
            print("Error. Duplicate sample ID. Exiting.")
            exit()
    sample_ids = copy.deepcopy(sample_ids_dupl)
    sample_ids = sorted(sample_ids, key=int)
    # Create a TSV for the comp_dict, which includes the percentages of the species.
    if comp_sole_dict is not None:
        tsv_dict = copy.deepcopy(comp_sole_dict)
    else:
        tsv_dict = copy.deepcopy(comp_dict)
    species_rel_abu_path = "{}/{}_species_rel_abu.tsv".format(stats_dir_path, label)
    species_rel_abu_file = open(species_rel_abu_path, "w")
    species_rel_abu_file.write("sample ID\ttaxid\trelative abundance (%)\n")
    for key_sample in tsv_dict.keys():
        for key_taxonid in tsv_dict[key_sample].keys():
            rel_abu = tsv_dict[key_sample][key_taxonid]
            species_rel_abu_file.write("{}\t{}\t{}\n".format(key_sample, key_taxonid, rel_abu))
    species_rel_abu_file.close()
    # Compute sensitivity
    comb_stats_dict = {}
    for si in sample_ids:
        if si in comp_dict.keys():
            comb_stats_dict, common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items = basic_stats(br_info_dict, comp_dict, comp_sole_dict, comb_stats_dict, si)
            write_species(common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items, stats_dir_path, label, si)
    # Return value
    return comb_stats_dict


def design_grouped_plots(comb_info_dict, metric_label_dict_1, plot_dir_parh, stats_dir_path):
    print("Plotting set 1...\n")
    if os.path.exists(plot_dir_parh):
        shutil.rmtree(plot_dir_parh)
    os.mkdir(plot_dir_parh)
    # Font sizes
    fs_num_1 = 25
    fs_num_2 = 20
    fs_num_3 = 20
    # Metrics to skip and groups of metrics based on the value range of the y axis.
    group_pass = ["gold_species_number", "common_species", "unique_species_total", "unique_gold", "unique_predicted", "True Negative (TN)"]
    group_0_100 = ["sensitivity", "precision", "accuracy"]
    group_0_1 = ["f1_score", "jaccard_index"]
    group_percentages = ["sensitivity", "precision", "accuracy"]
    group_logs = ["predicted_species_number", "true_positive", "false_positive", "false_negative"]
    # The values of each metric are collected from each type of analysis.
    # Dataframe:
    # metric: F1
    # sample    method   value
    # 1         kraken 8        0.02
    # 1         kraken 16       0.15
    # 1         kraken 72       0.2
    # 1         kraken 72 f     0.2
    # 1         metabinner 50   0.15
    # 1         metabinner 90   0.26
    # 1         metabinner nr   0.1
    # 1         comebin 50      0.14
    # 1         comebin 90      0.29
    # 1         comebin nr      0.06
    # 2         kraken 8        0.09
    # 2         kraken 16       0.19
    # 2         kraken 72       0.78
    # ...
    # Dictionary with information for the metrics.
    metric_dict = {}
    for key in comb_info_dict.keys():
        info_dict = comb_info_dict[key]
        if info_dict:
            for key_sample in info_dict.keys():
                for key_metric in info_dict[key_sample].keys():
                    metric_value = info_dict[key_sample][key_metric]
                    tmp_values = [key_sample, key, metric_value]
                    if key_metric not in metric_dict.keys():
                        metric_dict[key_metric] = [tmp_values]
                    else:
                        metric_dict[key_metric].append(tmp_values)
    # Pandas dataframe.
    pandas_dict = {}
    col_labels = ["sample", "method", "value"]
    for key_metric in metric_dict.keys():
        metric_df = pd.DataFrame(columns=col_labels)
        row_index = 0
        for item in metric_dict[key_metric]:
            metric_df.loc[row_index] = item
            row_index += 1
        pandas_dict[key_metric] = metric_df
    # Find all the different thresholds used in the analysis.
    combination_types = []
    for key_metric in metric_dict.keys():
        metric_category = metric_dict[key_metric]
        for sample_item_search in metric_category:
            cur_combination_seach = sample_item_search[1]
            if cur_combination_seach not in combination_types:
                combination_types.append(cur_combination_seach)
    # Sort the samples based on their species abundances and in turn based on their biases.
    sab_sort_order = [19, 18, 8, 11, 5, 6, 7, 13, 15, 10, 1, 4, 14, 12, 3, 16, 17, 2, 9]
    sample_id_label_dict = {
        "1": "1-A",
        "2": "2-N",
        "3": "3-G",
        "4": "4-G",
        "5": "5-A",
        "6": "6-G",
        "7": "7-N",
        "8": "8-G",
        "9": "9-N",
        "10": "10-N",
        "11": "11-N",
        "12": "12-A",
        "13": "13-A",
        "14": "14-N",
        "15": "15-G",
        "16": "16-N",
        "17": "17-N",
        "18": "18-A",
        "19": "19-N"
    }
    metric_dict_sab = {}
    for key_metric in metric_dict.keys():
        metric_category = metric_dict[key_metric]
        metric_dict_sab[key_metric] = []
        # For each combination type:
        for cmbt in combination_types:
            for sabid in sab_sort_order:
                # For each sample ID in the order of sample IDs:
                # Parse all the items of the current metric until finding
                # the item with the specific sample ID and combination type.
                for sample_item_check in metric_category:
                    cur_sample_id = int(sample_item_check[0])
                    cur_combination_check = sample_item_check[1]
                    if (sabid == cur_sample_id) and (cmbt == cur_combination_check):
                        # Modifying the sample ID fron integer to string in the item as when plotted not to be
                        # sorted automatically in ascending order.
                        cur_sample_id_str = str(cur_sample_id)
                        sample_label = sample_id_label_dict[cur_sample_id_str]
                        sample_item_check_mod = copy.deepcopy(sample_item_check)
                        sample_item_check_mod[0] = sample_label
                        metric_dict_sab[key_metric].append(sample_item_check_mod)
                        break
    # Pandas dataframe based on the order of the samples by species abundances and bisases.
    pandas_sab_dict = {}
    col_labels = ["sample", "method", "value"]
    for key_sab_metric in metric_dict_sab.keys():
        metric_sab_df = pd.DataFrame(columns=col_labels)
        sab_row_index = 0
        for sab_item in metric_dict_sab[key_sab_metric]:
            metric_sab_df.loc[sab_row_index] = sab_item
            sab_row_index += 1
        pandas_sab_dict[key_sab_metric] = metric_sab_df
    # Creating the plots
    for key_metric in pandas_dict.keys():
        if key_metric in group_pass:
            continue
        df_metric = pandas_dict[key_metric]
        # File path
        metric_label = metric_label_dict_1[key_metric]
        metric_label_sep = metric_label.split(" ")
        metric_label_con = "_".join(metric_label_sep)
        metric_tsv_path = "{}/{}.tsv".format(stats_dir_path, metric_label_con)
        # Storing the pandas dataframe as a csv.
        df_metric_sorted = df_metric.sort_values(by='sample')
        df_metric_sorted.to_csv(metric_tsv_path, sep="\t", index=False)
        # Figure
        figure, cur_axis = plt.subplots(figsize=(20, 12))
        # Plotting the pandas dataframe: Grouped barplot
        sns.barplot(data=df_metric, x="sample", y="value", hue="method", palette="colorblind", ax=cur_axis)
        # Place vertical lines to divie the sample barplots.
        for i in range(1, 19):
            cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
        # Legend title and legend
        legend_title = "Taxonomy Method\nand Database"
        legend_props = FontProperties(weight='bold', size=fs_num_2)
        cur_axis.legend(title=legend_title, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': fs_num_2}, title_fontproperties=legend_props)
        # Axis labels
        metric_label = metric_label_dict_1[key_metric]
        x_axis_label = "Sample ID"
        if key_metric in group_percentages:
            y_axis_label = "{} (%)".format(metric_label)
        else:
            y_axis_label = metric_label
        cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
        cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
        axis_title_label = "{} vs Sample ID".format(metric_label)
        cur_axis.set_title(axis_title_label, fontsize=fs_num_1, fontweight='bold')
        # Set the value range for the y axis.
        if key_metric in group_0_100:
            cur_axis.set_ylim(0, 100)
        elif key_metric in group_0_1:
            cur_axis.set_ylim(0, 1)
        elif key_metric in group_logs:
            cur_axis.set_yscale('log')
            cur_axis.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        # Changing the size of the tick labels
        cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
        # Layout: Padding from left, bottom, right, top
        figure.tight_layout()
        # File paths
        metric_label = metric_label_dict_1[key_metric]
        metric_label_sep = metric_label.split(" ")
        metric_label_con = "_".join(metric_label_sep)
        plot_png_path = "{}/{}.png".format(plot_dir_parh, metric_label_con)
        plot_jpg_path = "{}/{}.jpg".format(plot_dir_parh, metric_label_con)
        # Store the plots
        figure.savefig(plot_png_path)
        figure.savefig(plot_jpg_path)
        plt.close(figure)
    return pandas_dict, pandas_sab_dict
    

def design_metric_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict_1, plot_dir_parh, sab_status):
    if sab_status:
        print("Plotting set 2 for sorted samples based on their species abundances and biases....\n")
    else:
        print("Plotting set 2 for sorted samples based on their sample IDs....\n")
    # Font sizes
    fs_num_1 = 25
    fs_num_2 = 20
    fs_num_3 = 20
    # Letters for annotation
    annotation_letters = list("abcdefghijklmnopqrstuvwxyz")
    # Groups of metrics based on the value range of the y axis.
    group_0_100 = ["sensitivity", "precision", "accuracy"]
    group_0_1 = ["f1_score", "jaccard_index"]
    group_percentages = ["sensitivity", "precision", "accuracy"]
    group_logs = ["predicted_species_number", "true_positive", "false_positive", "false_negative"]
    # Creating the plots
    for key_group_metric in metric_group_dict.keys():
        cur_metric_group = metric_group_dict[key_group_metric]
        row_num = len(cur_metric_group)
        # Create the figure to store the subplots.
        figure, axis = plt.subplots(row_num, 1, figsize=(20, 5 * row_num)) 
        row_fig_index = 0
        for item in cur_metric_group:
            cur_axis = axis[row_fig_index]
            df_metric = pandas_dict[item]
            # Collecting the unique values of the sample column from the dataframe.
            sns.barplot(ax=cur_axis, data=df_metric, x="sample", y="value", hue="method", palette="colorblind", width=0.7)
            # Removes the right, top, left and bottom axis.
            sns.despine(ax=cur_axis, left=True)
            # Labels for x axis, y axis and title.
            metric_label = metric_label_dict_1[item]
            x_axis_label = "Sample ID"
            if item in group_percentages:
                y_axis_label = "{} (%)".format(metric_label)
            else:
                y_axis_label = metric_label
            axis_title_label = "{} vs Sample ID and Bias Combinations".format(metric_label)
            # Labels
            cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
            # Tick label size
            plt.xticks(fontsize=10)
            # Title
            cur_axis.set_title(axis_title_label, pad=20, loc='center', fontsize=fs_num_1, fontweight='bold')
            # Letter
            cur_axis.annotate(annotation_letters[row_fig_index], xy=(0.02, 1.10), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
            # Remove the legend, if for the last plot.
            if row_fig_index < row_num:
                cur_axis.legend().remove()
            # Set the value range for the y axis.
            if item in group_0_100:
                cur_axis.set_ylim(0, 100)
            elif item in group_0_1:
                cur_axis.set_ylim(0, 1)
            elif item in group_logs:
                cur_axis.set_yscale('log')
                # Compute the number of digits of the maximum value. Add one, this equals x. Find
                # the closes 10^x number.
                max_value = df_metric['value'].max()
                max_value_int = int(max_value)
                max_value_str = str(max_value_int)
                digit_num = len(max_value_str)
                y_axis_max_limit_mod = 10**digit_num
                # The value of 1 is 0 at log scale.
                cur_axis.set_ylim(bottom=1, top=y_axis_max_limit_mod)
                cur_axis.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
            # Replace the y axis ticks. Helps to rewrite ticks. For example the plot may continue on above the highest tick (tick with a value)
            # and thus has a bit more empty space without ticks. By rewriting the ticks a new tick will be provided to the top of the empty space.
            cur_y_ticks = cur_axis.get_yticks()
            cur_y_ticks_list = list(cur_y_ticks)
            cur_axis.set_yticks(cur_y_ticks_list)
            # Place vertical lines to divie the sample barplots.
            for i in range(1, 19):
                if sab_status:
                    if i in [3, 6, 9, 12, 15, 19]:
                        cur_axis.axvline(x=i-0.5, linewidth=3, color='black')
                    else:
                        cur_axis.axvline(x=i-0.5, linewidth=1, color='black')
                else:
                    cur_axis.axvline(x=i-0.5, linewidth=1, color='black')
            # Changing the size of the tick labels
            cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
            # Increasing the row index that shows which subplot will be added in the figure.
            row_fig_index += 1
            # Add space between the barplots.
            figure.subplots_adjust(wspace=1.0)
        # Add one legend at the bottom. Make the legend as much horizontal as possible.
        # Legend
        legend_title = "Taxonomy Method and Database"
        legend_props = FontProperties(weight='bold', size=fs_num_2)
        cur_axis.legend(title=legend_title, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=4, prop={'size': fs_num_2}, title_fontproperties=legend_props)
        # Layout: Padding from left, bottom, right, top
        figure.tight_layout()
        # Space between the plots of the figure.
        figure.subplots_adjust(hspace=0.6)
        # File paths
        if sab_status:
            group_metrics_path_png = "{}/group_metrics_sab_{}.png".format(plot_dir_parh, key_group_metric)
            group_metrics_path_jpg = "{}/group_metrics_sab_{}.jpg".format(plot_dir_parh, key_group_metric)
        else:
            group_metrics_path_png = "{}/group_metrics_{}.png".format(plot_dir_parh, key_group_metric)
            group_metrics_path_jpg = "{}/group_metrics_{}.jpg".format(plot_dir_parh, key_group_metric)
        # Saving the figure.
        figure.savefig(group_metrics_path_png)
        figure.savefig(group_metrics_path_jpg)
        # Closing the figure.
        plt.close(figure)


def rank_stats(pandas_dict, stats_dir_path, metric_label_dict_1, metric_label_dict_2):
    print("Ranking the methods...\n")
    # A group of methods is skipped, a group which contains each method that shows higher "similarity" of the two groups while its scores is increased
    # and a group which contains each method that shows higher "similarity" of the two groups while its scores is decreased.
    group_pass = ["gold_species_number", "predicted_species_number", "common_species", "unique_species_total", "unique_gold", "unique_predicted", "True Negative (TN)"]
    group_max = ["true_positive", "sensitivity", "precision", "accuracy", "f1_score", "jaccard_index"]
    group_min = ["false_positive", "false_negative", "l1_norm"]
    # Write all computed metrics for all samples in a TSV file.
    all_metrics_path = "{}/all_metrics.tsv".format(stats_dir_path)
    all_metrics_file = open(all_metrics_path, "w")
    for key_metric in pandas_dict.keys():
        if key_metric in group_pass:
            continue
        metric_label = metric_label_dict_2[key_metric]
        all_metrics_file.write("{}\n".format(metric_label))
        metric_df = pandas_dict[key_metric]
        metric_df_sorted = metric_df.sort_values(by='sample')
        metric_df_sorted.to_csv(all_metrics_file, sep="\t", index=False)
        all_metrics_file.write("\n")
    all_metrics_file.close()
    # The higher ranking methods across for each sample are identified.
    method_sample_rank_dict = {}
    for key_metric in pandas_dict.keys():
        if key_metric in group_pass:
            continue
        metric_df = pandas_dict[key_metric]
        method_sample_rank_dict[key_metric] = {}
        # Label and file path.
        filename_metric_label = metric_label_dict_1[key_metric]
        filename_metric_label_sep = filename_metric_label.split(" ")
        filename_metric_label_con = "_".join(filename_metric_label_sep)
        sorted_df_tsv_path ="{}/{}_sorted.tsv".format(stats_dir_path, filename_metric_label_con)
        sorted_df_tsv_file = open(sorted_df_tsv_path, "w")
        first_metric = True
        for sample_id in range(1, 20):
            sample_id_list = [sample_id]
            # Isolating the value from the sample column.
            df_metric_method = metric_df[metric_df['sample'].isin(sample_id_list)]
            if df_metric_method.empty:
                continue
            # Sorting the dataframe.
            df_sorted = df_metric_method.sort_values(by='value', ascending=False)
            # Identifying the methods with the highest or lowest value for the current metric across all samples.
            if key_metric in group_max:
                val_best = df_sorted['value'].max()
            elif key_metric in group_min:
                val_best = df_sorted['value'].min()
            else:
                print(key_metric)
                print("Metric was not found in the groups of determined lowest or highest value metrics. Exiting.\n")
                exit()
            max_val_methods = df_sorted[df_sorted['value'] == val_best]['method'].tolist()
            if first_metric:
                df_sorted.to_csv(sorted_df_tsv_file, sep="\t", index=False)
                first_metric = False
            else:
                df_sorted.to_csv(sorted_df_tsv_file, sep="\t", header=False, index=False)
            sorted_df_tsv_file.write("\n")
            # Add each method to the current metric.
            method_sample_rank_dict[key_metric][sample_id] = max_val_methods
    sorted_df_tsv_file.close()
    # Convert the dictionary to a dataframe.
    col_labels = ["metric", "sample", "top_methods"]
    method_sample_rank_df = pd.DataFrame(columns=col_labels)
    row_index = 0
    for key_metric in method_sample_rank_dict.keys():
        for key_sample in method_sample_rank_dict[key_metric].keys():
            top_methods = method_sample_rank_dict[key_metric][key_sample]
            top_methods_str = ",".join(top_methods)
            metric_label = metric_label_dict_2[key_metric]
            temp_list = [metric_label, key_sample, top_methods_str]
            method_sample_rank_df.loc[row_index] = temp_list
            row_index += 1
    # Store the dataframe.
    top_methods_tsv_path = "{}/sample_top_methods.tsv".format(stats_dir_path)
    method_sample_rank_df.to_csv(top_methods_tsv_path, sep="\t", index=False)
    # Top method counting across all samples for each metric.
    method_allsamples_count_dict = {}
    for key_metric in method_sample_rank_dict.keys():
        method_allsamples_count_dict[key_metric] = {}
        for key_sample in method_sample_rank_dict[key_metric].keys():
            temp_list = method_sample_rank_dict[key_metric][key_sample]
            for top_method in temp_list:
                if top_method not in method_allsamples_count_dict[key_metric].keys():
                    method_allsamples_count_dict[key_metric][top_method] = 1
                else:
                    method_allsamples_count_dict[key_metric][top_method] += 1
    # Convert the dictionary to a dataframe.
    col_labels = ["metric", "top_method", "frequency"]
    method_allsamples_rank_df = pd.DataFrame(columns=col_labels)
    row_index = 0
    for key_metric in method_allsamples_count_dict.keys():
        for key_top_method in method_allsamples_count_dict[key_metric].keys():
            top_method_freq = method_allsamples_count_dict[key_metric][key_top_method]
            metric_label = metric_label_dict_2[key_metric]
            temp_list = [metric_label, key_top_method, top_method_freq]
            method_allsamples_rank_df.loc[row_index] = temp_list
            row_index += 1
    # Sort the dataframe.
    method_allsamples_rank_sorted_df = method_allsamples_rank_df.sort_values(by=['metric', 'frequency'], ascending=[True, False])
    # Store the dataframe
    count_top_methods_tsv_path = "{}/total_top_methods.tsv".format(stats_dir_path)
    method_allsamples_rank_sorted_df.to_csv(count_top_methods_tsv_path, sep="\t", index=False)


def check_add(add_item, target_num):
    if add_item != "None":
        add_item = float(add_item)
        target_num += add_item
    return target_num


def collect_times(time_stats_dir_path):
    pre_time_dict = {}
    time_file_names = os.listdir(time_stats_dir_path)
    for tfn in time_file_names:
        if "time_analysis" in tfn:
            tfn_spliited = tfn.split("_")
            method_dp = tfn_spliited[1]
            sample_id = tfn_spliited[3]
            # Dictionary information
            if sample_id not in pre_time_dict.keys():
                pre_time_dict[sample_id] = {}
            pre_time_dict[sample_id][method_dp] = {}
            # File path
            tfn_path = "{}/{}".format(time_stats_dir_path, tfn)
            # File lines
            tfn_lines = read_file(tfn_path)
            for line in tfn_lines:
                line_splitted = line.split("\t")
                time_label = line_splitted[0]
                if line_splitted[1] != "None":
                    time_period = float(line_splitted[1])
                else:
                    time_period = "None"
                pre_time_dict[sample_id][method_dp][time_label] = time_period
    # Order the time dictionary.
    orderded_samples = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    orderded_methods = ["k8", "k16", "k72", "cnr", "mnr"]
    time_dict = {}
    for key_sample in orderded_samples:
        for key_method_source in orderded_methods:
            if key_sample in pre_time_dict.keys():
                for key_method_dest in pre_time_dict[key_sample].keys():
                    if key_method_source == key_method_dest:
                        if key_sample not in time_dict.keys():
                            time_dict[key_sample] = {}
                        time_dict[key_sample][key_method_source] = pre_time_dict[key_sample][key_method_dest]
    return time_dict


def time_stages(time_dict):
    # For the added time periods.
    time_full_dict = copy.deepcopy(time_dict)
    # K8: All stages
    # K16: After gene prediction: Add = FastQC_Initial + BBDuk + FastQc_Final + Assembly + Gene_Prediction
    # K72: After gene prediction: Add = FastQC_Initial + BBDuk + FastQc_Final + Assembly + Gene_Prediction
    # MNR: After gene prediction: Add = FastQC_Initial + BBDuk + FastQc_Final + Assembly + Gene_Prediction
    # CNR: After gene prediction: Add = FastQC_Initial + BBDuk + FastQc_Final + Assembly + Gene_Prediction
    # The earlier stages to be added to the other methods are taken from K8. The gene annotation time (not gene prediction time)
    # is only run by ProteoSeeker for the seek mode, thus in the taxonomy mode it does not run.
    # The total time is computed based on the previous total plus any common stage with the Kraken8 method.
    # For the Kraken8 method is the same as the previous total time.
    # This total time however will not be exactly equal to the sum of the stage-specific times of the full time table for a method.
    # Even with the stage-specific times filled by the Kraken8 method, the total time will be slightly higher because of the intermediate
    # parts of the pipeline which are not computed by these main processes whose times are computed in the different of the full time table.
    # ---------sample_{}_total_time.tsv---------
    # time_tool_dict: Contains the total time of each method based on its previous total time plus the stages from Kraken8
    # time_tool_dict -> total time converted from seconds to mins -> df_total_time_dict -> sample_{}_total_time.tsv
    # ---------sample_{}_full_time.tsv---------
    # time_full_dict: Contains the total time of each method based on its previous total time plus the stages from Kraken8 and each stage-specific time.
    # time_full_dict: Then, it gets the new stage-specific times from Kraken8 for all common stages.
    # time_full_dict: Then, its times are converted to mins and are combined (summed up) to main stages (e.g., first and second FastQC analysis to one) -> df_full_time_dict
    # df_full_time_dict -> sample_{}_full_time
    # Therefore, the full times are the combined (or not) times of stage-specific times. Their sum is expected not to match the total time of the method due to the
    # presence of the intermediate processes running between the main stages in the pipeline of ProteoSeeker.
    time_tool_dict = {}
    for key_sample in time_dict.keys():
        # Information to dictionary
        time_tool_dict[key_sample] = {}
        sum_common_times = None
        # Sum 1
        if "k8" in time_dict[key_sample].keys():
            k8_fastqc_init_time = time_dict[key_sample]["k8"]["fastqc_initial_time"]
            k8_preproc_time = time_dict[key_sample]["k8"]["preprocessing_time"]
            k8_fastqc_final_time = time_dict[key_sample]["k8"]["fastqc_final_time"]
            k8_assembly_time = time_dict[key_sample]["k8"]["assembly_time"]
            k8_gene_pred_time = time_dict[key_sample]["k8"]["gene_prediction_time"]
            if (k8_fastqc_init_time == "None") or (k8_preproc_time == "None") or (k8_fastqc_final_time == "None") or (k8_assembly_time == "None") or (k8_gene_pred_time == "None"):
                sum_common_times = None
            else:
                k8_fastqc_init_time = float(k8_fastqc_init_time)
                k8_preproc_time = float(k8_preproc_time)
                k8_fastqc_final_time = float(k8_fastqc_final_time)
                k8_assembly_time = float(k8_assembly_time)
                k8_gene_pred_time = float(k8_gene_pred_time)
                sum_common_times = k8_fastqc_init_time + k8_preproc_time + k8_fastqc_final_time + k8_assembly_time + k8_gene_pred_time
        # Perform the additions
        if sum_common_times is None:
            print("Error in computing the additional times. Exiting.\n")
            exit()
        else:
            # Add the new tool time in the dictionary which stored only tool times and also
            # change the tool time in the dicionary with all the time periods of the tool.
            for key_method in time_dict[key_sample].keys():
                tool_time = float(time_dict[key_sample][key_method]["tool_time"])
                if key_method == "k8":
                    new_tool_time = round(tool_time, 2)
                    time_full_dict[key_sample][key_method]["tool_time"] = new_tool_time
                elif key_method in ["k16", "k72", "cnr", "mnr"]:
                    new_tool_time = tool_time + sum_common_times
                    new_tool_time = round(new_tool_time, 2)
                    time_full_dict[key_sample][key_method]["tool_time"] = new_tool_time
                else:
                    print("Error. Unknown method: {}. Exiting.".format(key_method))
                    exit()
                time_tool_dict[key_sample][key_method] = new_tool_time
    return time_full_dict, time_tool_dict


def fill_times(time_dict, time_full_dict):
    # Each stage-specific execution time from Kraken 8 untill the stage of gene prediction is provided to the other methods.
    # The new sum for the total time of the method has already been computed and provided to the method in the "time_stages" function.
    for key_sample in time_dict.keys():
        if "k8" in time_dict[key_sample].keys():
            k8_fastqc_init_time = time_dict[key_sample]["k8"]["fastqc_initial_time"]
            k8_preproc_time = time_dict[key_sample]["k8"]["preprocessing_time"]
            k8_fastqc_final_time = time_dict[key_sample]["k8"]["fastqc_final_time"]
            k8_assembly_time = time_dict[key_sample]["k8"]["assembly_time"]
            k8_gene_pred_time = time_dict[key_sample]["k8"]["gene_prediction_time"]
        for key_method in time_full_dict[key_sample].keys():
            if key_method == "k8":
                continue
            time_full_dict[key_sample][key_method]["fastqc_initial_time"] = k8_fastqc_init_time
            time_full_dict[key_sample][key_method]["preprocessing_time"] = k8_preproc_time
            time_full_dict[key_sample][key_method]["fastqc_final_time"] = k8_fastqc_final_time
            time_full_dict[key_sample][key_method]["assembly_time"] = k8_assembly_time
            time_full_dict[key_sample][key_method]["gene_prediction_time"] = k8_gene_pred_time
    return time_full_dict


def convert_to_minutes(time_seconds):
    time_minutes = time_seconds / 60
    time_minutes = round(time_minutes, 2)
    time_hours = time_minutes / 60
    time_hours = round(time_hours, 2)
    return time_minutes, time_hours


def create_df_stacked(time_full_dict):
    # Dataframe:
    # Periods of time for processes that no place to certain methods are computed as having run 0 seconds for those methods. Processes, which take up little time are grouped together as "other processes" or something else.
    # Sample ID
    # method    tool_time    preprocessing_time    assembly_time    gene_prediction_time    cd_hit_time    kraken_specific_time/binning_time    bowtie_time    bin_taxonomy_cm_time    kraken_binning_time
    # fastqc = fastqc_initial + fastqc_final
    # hmmer = hmmer_spec_time + hmmer_broad
    # blastp = blastp_fpd_no_doms_time + blastp_fpd_swiss_doms_time
    # cm_taxonomy = hmmer_spec_taxonomy_time + blastp_fpd_swiss_taxonomy_time + bin_analysis_time + bin_taxonomy_time
    # results = info_collection_time + results_time
    # The only taxonomy mode does not include: gene_annotation_time, topology_time    motifs_time    family_prediction_time
    all_col_labels = ["sample", "method", "quality_control", "preprocessing", "assembly", "gene_prediction", "clustering", "binning", "taxonomy", "read_alignment", "results"]
    df_full_time_all = pd.DataFrame(columns=all_col_labels)
    df_full_time_dict = {}
    row_index_all = 0
    for key_sample in time_full_dict.keys():
        # Initialize an empty Pandas dataframe.
        col_labels = ["method", "quality_control", "preprocessing", "assembly", "gene_prediction", "clustering", "binning", "taxonomy", "read_alignment", "results"]
        time_metric_df = pd.DataFrame(columns=col_labels)
        row_index = 0
        for key_method in time_full_dict[key_sample].keys():
            df_total_time_str = time_full_dict[key_sample][key_method]["tool_time"]
            df_fastqc_init_time_str = time_full_dict[key_sample][key_method]["fastqc_initial_time"]
            df_preproc_time_str = time_full_dict[key_sample][key_method]["preprocessing_time"]
            df_fastqc_final_time_str = time_full_dict[key_sample][key_method]["fastqc_final_time"]
            df_assembly_time_str = time_full_dict[key_sample][key_method]["assembly_time"]
            df_gene_pred_time_str = time_full_dict[key_sample][key_method]["gene_prediction_time"]
            df_cd_hit_time_str = time_full_dict[key_sample][key_method]["cd_hit_time"]
            df_kraken_tax_time_str = time_full_dict[key_sample][key_method]["kraken_time"]
            df_cm_binning_time_str = time_full_dict[key_sample][key_method]["binning_time"]
            df_bowtie_time_str = time_full_dict[key_sample][key_method]["bowtie_time"]
            df_hmmer_tax_spec_time_str = time_full_dict[key_sample][key_method]["hmmer_spec_taxonomy_time"]
            df_blastp_3_time_str = time_full_dict[key_sample][key_method]["blastp_fpd_swiss_taxonomy_time"]
            df_cm_analysis_time_str = time_full_dict[key_sample][key_method]["bin_analysis_cm_time"]
            df_cm_tax_time_str = time_full_dict[key_sample][key_method]["bin_taxonomy_cm_time"]
            df_kraken_binning_time_str = time_full_dict[key_sample][key_method]["kraken_binning_time"]
            df_info_time_str = time_full_dict[key_sample][key_method]["info_collection_time"]
            df_results_time_str = time_full_dict[key_sample][key_method]["results_time"]
            # Initialize sums
            df_total_time = 0
            df_fastqc_time = 0
            df_preproc_time = 0
            df_assembly_time = 0
            df_gene_pred_time = 0
            df_cd_hit_time = 0
            df_kraken_taxonomy_time = 0
            df_kraken_binning_time = 0
            df_cm_taxonomy_time = 0
            df_cm_binning_time = 0
            df_bowtie_time = 0
            df_results_time = 0
            df_total_time = check_add(df_total_time_str, df_total_time)
            df_fastqc_time = check_add(df_fastqc_init_time_str, df_fastqc_time)
            df_fastqc_time = check_add(df_fastqc_final_time_str, df_fastqc_time)
            df_preproc_time = check_add(df_preproc_time_str, df_preproc_time)
            df_assembly_time = check_add(df_assembly_time_str, df_assembly_time)
            df_gene_pred_time = check_add(df_gene_pred_time_str, df_gene_pred_time)
            df_cd_hit_time = check_add(df_cd_hit_time_str, df_cd_hit_time)
            df_kraken_taxonomy_time = check_add(df_kraken_tax_time_str, df_kraken_taxonomy_time)
            df_kraken_binning_time = check_add(df_kraken_binning_time_str, df_kraken_binning_time)
            df_cm_binning_time = check_add(df_cm_binning_time_str, df_cm_binning_time)
            df_cm_taxonomy_time = check_add(df_hmmer_tax_spec_time_str, df_cm_taxonomy_time)
            df_cm_taxonomy_time = check_add(df_blastp_3_time_str, df_cm_taxonomy_time)
            df_cm_taxonomy_time = check_add(df_cm_analysis_time_str, df_cm_taxonomy_time)
            df_cm_taxonomy_time = check_add(df_cm_tax_time_str, df_cm_taxonomy_time)
            df_bowtie_time = check_add(df_bowtie_time_str, df_bowtie_time)
            df_results_time = check_add(df_info_time_str, df_results_time)
            df_results_time = check_add(df_results_time_str, df_results_time)
            # Converting seconds to minutes.
            # df_total_time_m, df_total_time_h = convert_to_minutes(df_total_time)
            df_fastqc_time_m, df_fastqc_time_h = convert_to_minutes(df_fastqc_time)
            df_preproc_time_m, df_preproc_time_h = convert_to_minutes(df_preproc_time)
            df_assembly_time_m, df_assembly_time_h = convert_to_minutes(df_assembly_time)
            df_gene_pred_time_m, df_gene_pred_time_h = convert_to_minutes(df_gene_pred_time)
            df_cd_hit_time_m, df_cd_hit_time_h = convert_to_minutes(df_cd_hit_time)
            df_kraken_binning_time_m, df_kraken_binning_time_h = convert_to_minutes(df_kraken_binning_time)
            df_kraken_taxonomy_time_m, df_kraken_taxonomy_time_h = convert_to_minutes(df_kraken_taxonomy_time)
            df_cm_binning_time_m, df_cm_binning_time_h = convert_to_minutes(df_cm_binning_time)
            df_cm_taxonomy_time_m, df_cm_taxonomy_time_h = convert_to_minutes(df_cm_taxonomy_time)
            df_bowtie_time_m, df_bowtie_time_h = convert_to_minutes(df_bowtie_time)
            df_results_time_m, df_results_time_h = convert_to_minutes(df_results_time)
            # Passing the information to a list.
            temp_list = [key_method, df_fastqc_time_m, df_preproc_time_m, df_assembly_time_m, df_gene_pred_time_m, df_cd_hit_time_m]
            if key_method in ["k8", "k16", "k72"]:
                temp_list.append(df_kraken_binning_time_m)
                temp_list.append(df_kraken_taxonomy_time_m)
            else:
                temp_list.append(df_cm_binning_time_m)
                temp_list.append(df_cm_taxonomy_time_m)
            temp_list += [df_bowtie_time_m, df_results_time_m]
            all_temp_list = copy.deepcopy(temp_list)
            key_sample_int = int(key_sample)
            all_temp_list = [key_sample_int] + all_temp_list
            # Set the list as the row to the current dataframe index (row index).
            time_metric_df.loc[row_index] = temp_list
            df_full_time_all.loc[row_index_all] = all_temp_list
            row_index += 1
            row_index_all += 1
        # Set the first column of the dataframe to be the indexes (labels) of the rows.
        time_metric_df.set_index('method', inplace=True)
        # Store the dataframe to a dictionary.
        df_full_time_dict[key_sample] = time_metric_df
    return df_full_time_dict, df_full_time_all


def create_df_total_time(time_tool_dict):
    # A dataframe for each sample.
    df_total_time_dict = {}
    for key_sample in time_tool_dict.keys():
        # Initialize an empty Pandas dataframe.
        col_labels = ["method", "total_time"]
        time_metric_df_temp = pd.DataFrame(columns=col_labels)
        row_index = 0
        for key_method in time_tool_dict[key_sample].keys():
            total_time = time_tool_dict[key_sample][key_method]
            total_time_m = total_time / 60
            total_time_m = round(total_time_m, 2)
            temp_list = [key_method, total_time_m]
            time_metric_df_temp.loc[row_index] = temp_list
            row_index += 1
        df_total_time_dict[key_sample] = time_metric_df_temp
    # A dataframe for all samples.
    # Initialize an empty Pandas dataframe.
    col_labels = ["sample", "method", "total_time"]
    total_time_all_df = pd.DataFrame(columns=col_labels)
    row_index = 0
    for key_sample in time_tool_dict.keys():
        for key_method in time_tool_dict[key_sample].keys():
            total_time = time_tool_dict[key_sample][key_method]
            total_time_m = total_time / 60
            total_time_m = round(total_time_m, 2)
            key_sample_int = int(key_sample)
            temp_list = [key_sample_int, key_method, total_time_m]
            total_time_all_df.loc[row_index] = temp_list
            row_index += 1
    return df_total_time_dict, total_time_all_df


def count_common_time(df_full_time_dict):
    # Counter to sum the time up untill the stage of megahit.
    common_time_dict = {}
    for key_sample in df_full_time_dict.keys():
        full_time_metric_df = df_full_time_dict[key_sample]
        # Retaining the columns up to the "gene_prediction" column.
        up_to_genpred_df = full_time_metric_df.loc[:, :"gene_prediction"]
        # Summing the values of the columns across each row and storing the results in a dataframe compomsed of one column with the headers of the
        # columns of the "up_to_genpred_df" and one column with the sums.
        row_sums_df = up_to_genpred_df.sum(axis=1)
        # Converting the column with the sums to a list.
        common_stages_sum_list = row_sums_df.tolist()
        # Check if all sums are the same. If not exit with error.
        item_0 = common_stages_sum_list[0]
        for item in common_stages_sum_list:
            if item_0 != item:
                print("Error. Different sum of time up untill the megahit stage. Exiting.\n")
                exit()
        item_0 = float(item_0)
        # Rounding the sum. 
        common_stages_sum = round(item_0, 2)
        # This sum is stored as the common time period for a given sample.
        common_time_dict[key_sample] = common_stages_sum
    return common_time_dict


def plot_full_sample_time(df_full_time_dict, common_time_dict, time_dir, time_stats_dir_path):
    print("Full time plots for each sample...\n")
    # Font sizes
    fs_num_1 = 12
    fs_num_2 = 10
    fs_num_3 = 10
    # Grouped barplots
    for key_sample in df_full_time_dict.keys():
        # Get the dataframe from the dictionary.
        full_time_metric_df = df_full_time_dict[key_sample]
        # Get the time for megahit.
        cur_megahit_time = common_time_dict[key_sample]
        # Create the filename.
        stacked_tsv_path = "{}/sample_{}_full_time.tsv".format(time_stats_dir_path, key_sample)
        # Save to CSV
        full_time_metric_df.to_csv(stacked_tsv_path, sep="\t")
        # Color palette
        color_num = len(full_time_metric_df.columns)
        cld_palette = sns.color_palette("colorblind", color_num + 1)
        # Removing the 9nth color wich is a very dim yellow.
        del cld_palette[8]
        # Create the plot
        ax = full_time_metric_df.plot(kind='bar', stacked=True, color=cld_palette)
        # Place a horizontal line at the sum of time up to megahit.
        ax.axhline(y=cur_megahit_time, linewidth=0.5, color='black')
        # Labels
        ax.set_xlabel('Taxonomy Route and Database', fontsize=fs_num_2, fontweight='bold')
        ax.set_ylabel('Execution Time (min)', fontsize=fs_num_2, fontweight='bold')
        ax.set_title('Execution Time vs Taxonomy Route and Database', fontsize=fs_num_1, fontweight='bold')
        # Place the legend outside the plot
        # Custom legend labels
        custom_labels = ["quality control", "preprcessing", "assembly", "gene prediction", "protein clustering", "binning", "taxonomy", "read alignment", "results"]
        # Get the current legend handles and labels
        legend_mods, labels = ax.get_legend_handles_labels()
        # Update the labels with custom labels
        for handle, new_label in zip(legend_mods, custom_labels):
            handle.set_label(new_label)        
        legend_props = FontProperties(weight='bold', size=fs_num_3)
        legend_title = "Stage"
        ax.legend(handles=legend_mods, title=legend_title, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': fs_num_3}, title_fontproperties=legend_props)
        # Changing the size of the tick labels
        ax.tick_params(axis='both', which='major', labelsize=fs_num_3)
        # Layout
        plt.tight_layout()
        # Image paths.
        full_time_df_png_path = "{}/sample_{}_full_time.png".format(time_dir, key_sample)
        full_time_df_jpg_path = "{}/sample_{}_full_time.jpg".format(time_dir, key_sample)
        # Save the plot
        plt.savefig(full_time_df_png_path)
        plt.savefig(full_time_df_jpg_path)
        # Close the plot
        plt.close()


def plot_total_sample_time(df_total_time_dict, common_time_dict, time_dir, time_stats_dir_path):
    print("Total time plot for each sample...\n")
    # Font sizes
    fs_num_1 = 12
    fs_num_2 = 10
    fs_num_3 = 10
    # Create a barplot for each sample.
    for key_sample in df_total_time_dict.keys():
        # Get the dataframe from the dictionary.
        total_time_sample_df = df_total_time_dict[key_sample]
        # Get the time for megahit.
        cur_megahit_time = common_time_dict[key_sample]
        # Create the filename.
        stacked_tsv_path = "{}/sample_{}_total_time.tsv".format(time_stats_dir_path, key_sample)
        # Save to CSV
        total_time_sample_df.to_csv(stacked_tsv_path, sep="\t", index=False)
        # Color palette
        color_num = len(total_time_sample_df)
        cld_palette = sns.color_palette("colorblind", color_num)
        ax = total_time_sample_df.plot(kind='bar', x='method', y='total_time', color=cld_palette, legend=False)
        # Create the plot
        # Place a horizontal line at the sum of time up to megahit.
        ax.axhline(y=cur_megahit_time, linewidth=0.5, color='black')
        # Labels
        ax.set_xlabel('Taxonomy Route and Database', fontsize=fs_num_2, fontweight='bold')
        ax.set_ylabel('Execution Time (min)', fontsize=fs_num_2, fontweight='bold')
        ax.set_title('Execution Time vs Taxonomy Route and Database', fontsize=fs_num_1, fontweight='bold')
        # Changing the size of the tick labels        
        ax.tick_params(axis='both', which='major', labelsize=fs_num_3)
        # Layout
        plt.tight_layout()
        # Image paths
        total_time_sample_df_png_path = "{}/sample_{}_total_time.png".format(time_dir, key_sample)
        total_time_sample_df_jpg_path = "{}/sample_{}_total_time.jpg".format(time_dir, key_sample)
        # Save the plot
        plt.savefig(total_time_sample_df_png_path)
        plt.savefig(total_time_sample_df_jpg_path)
        # Close the plot
        plt.close()


def plot_total_all_time(total_time_all_df, time_dir, time_stats_dir_path):
    print("Total time plot for all samples...\n")
    # Font sizes
    fs_num_1 = 17
    fs_num_2 = 15
    fs_num_3 = 10
    # Create a grouped barplot for the total times
    total_time_all_tsv_path = "{}/total_time_all.tsv".format(time_stats_dir_path)
    # Storing the pandas dataframe as a csv.
    total_time_all_df.to_csv(total_time_all_tsv_path, sep="\t", index=False)
    # Plotting the pandas dataframe: Barplot
    total_time_all_barplot = sns.catplot(data=total_time_all_df, kind="bar", x="sample", y="total_time", hue="method", palette="colorblind", height=6, aspect=2)
    # Remove legend.
    total_time_all_barplot._legend.remove()
    # Labels
    x_axis_label = "Sample ID"
    y_axis_label = "Execution Time (min)"
    title_label = "Execution Time vs Sample ID"
    legend_label = "Taxonomy Route and Database"
    total_time_all_barplot.set_axis_labels(x_axis_label, y_axis_label, fontsize=fs_num_2, fontweight='bold')
    plt.suptitle(title_label, fontsize=fs_num_1, fontweight='bold')
    # Legend properties
    legend_specs = FontProperties(weight='bold', size=fs_num_3)
    plt.legend(title=legend_label, title_fontproperties=legend_specs, loc='upper right', prop={'size': fs_num_3})
    # Changing the size of the tick labels        
    plt.tick_params(axis='both', which='major', labelsize=fs_num_2)
    # Layout
    plt.tight_layout()
    # File paths
    plot_png_path = "{}/total_time_all.png".format(time_dir)
    plot_jpg_path = "{}/total_time_all.jpg".format(time_dir)
    # Saving plots
    plt.savefig(plot_png_path)
    plt.savefig(plot_jpg_path)
    plt.close()


def plot_full_group_sample_time(df_full_time_all, time_dir, time_stats_dir_path, time_syn_sample_group_dict, time_sample_group_label_dict, methods_group):
    print("Full time plot for all samples grouped...\n")
    # Font sizes
    fs_num_1 = 25
    fs_num_2 = 20
    fs_num_3 = 20
    # Letters for annotation
    annotation_letters = list("abcdefghijklmnopqrstuvwxyz")
    # The number of the methods to plot.
    methods_num = len(methods_group)
    # Save the dataframe to a CSV file.
    full_time_sample_grouped_time_tsv_path = "{}/full_time_sample_grouped_time.tsv".format(time_stats_dir_path)
    df_full_time_all.to_csv(full_time_sample_grouped_time_tsv_path, sep="\t", index=False)
    # Create the figure to store the subplots.
    row_num = 6
    figure, axis = plt.subplots(row_num, 1, figsize=(20, 5 * row_num))
    row_fig_index = 0
    # Filter specific samples.
    for key_group in time_syn_sample_group_dict.keys():
        sample_id_group = time_syn_sample_group_dict[key_group]
        # Group label
        group_label = time_sample_group_label_dict[key_group]
        # Filtering the dataframe.
        df_time_filtered = df_full_time_all[df_full_time_all['sample'].isin(sample_id_group)]
        if df_time_filtered.empty:
            break
        # Create the plot
        cur_axis = axis[row_fig_index]
        # If the inplace option was set to True then the dataframe would be modified and would not return a new one.
        df_time_indexed = df_time_filtered.set_index(['sample', 'method'])
        # Color palette
        cld_palette = sns.color_palette("colorblind")
        # Removing the 9nth color wich is a very dim yellow.
        del cld_palette[8]
        # Plot
        df_time_indexed.plot(kind='bar', stacked=True, ax=cur_axis, color=cld_palette)
        # Loop through the index of df_time_indexed.
        x_labels = []
        # The index has already been set above based on the two columns of sampel and method.
        for cur_sample, cur_method in df_time_indexed.index:
            cur_label = "{}-{}".format(cur_sample, cur_method)
            x_labels.append(cur_label)
        cur_axis.set_xticklabels(x_labels, rotation=45, ha="right")
        # Labels for x axis, y axis and title.
        x_axis_label = "Sample ID and Database"
        y_axis_label = "Execution Time (min)"
        axis_title_label = "Execution Time vs Sample ID, Taxonomy Approach and Database Combinations: {}".format(group_label)
        # Labels
        cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
        cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
        # Title
        cur_axis.set_title(axis_title_label, pad=20, loc='center', fontsize=fs_num_1, fontweight='bold')
        # Letter
        cur_axis.annotate(annotation_letters[row_fig_index], xy=(0.02, 1.10), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
        # Changing the size of the tick labels        
        cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
        # Remove the legend, if for the last plot.
        if row_fig_index < row_num:
            cur_axis.legend().remove()
        # Replace the y axis ticks. Helps to rewrite ticks. For example the plot may continue on above the highest tick (tick with a value)
        # and thus has a bit more empty space without ticks. By rewriting the ticks a new tick will be provided to the top of the empty space.
        cur_y_ticks = cur_axis.get_yticks()
        cur_y_ticks_list = list(cur_y_ticks)
        cur_axis.set_yticks(cur_y_ticks_list)
        # Place vertical lines to divie the sample barplots.
        max_ver_line_col = len(sample_id_group) * methods_num
        for i in range(methods_num, max_ver_line_col, methods_num):
            cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=1.5)
        # Increasing the row index that shows which subplot will be added in the figure.
        row_fig_index += 1
    # Legend
    legend_title = "Stage"
    legend_props = FontProperties(weight='bold', size=fs_num_3)
    # Custom legend labels
    custom_labels = ["quality control", "preprcessing", "assembly", "gene prediction", "protein clustering", "binning", "taxonomy", "read alignment", "results"]
    # Get the legend handles and labels from the last axis with a plot.
    legend_mods, labels = cur_axis.get_legend_handles_labels()
    # Get the last axis.
    tsgd_keys = list(time_syn_sample_group_dict.keys())
    tsgd_keys_len = len(tsgd_keys) - 1
    last_axis = axis[tsgd_keys_len]
    # Update the labels with custom labels and use the legend from the last axis.
    for handle, new_label in zip(legend_mods, custom_labels):
        handle.set_label(new_label)
    last_axis.legend(handles=legend_mods, title=legend_title, loc='upper center', bbox_to_anchor=(0.5, -0.6), ncol=6, prop={'size': fs_num_3}, title_fontproperties=legend_props)
    # Padding from left, bottom, right, top.
    figure.tight_layout()
    # Space between the plots of the figure.
    figure.subplots_adjust(hspace=0.9)
    # File paths.
    full_time_grouped_path_png = "{}/full_time_sample_grouped.png".format(time_dir)
    full_time_grouped_path_jpg = "{}/full_time_sample_grouped.jpg".format(time_dir)
    # Saving the figure.
    figure.savefig(full_time_grouped_path_png)
    figure.savefig(full_time_grouped_path_jpg)
    # Closing the figure.
    plt.close(figure)


def plot_full_group_method_time(total_time_all_df, time_dir, time_syn_sample_group_dict, time_sample_group_label_dict, time_stats_dir_path):
    print("Full time plot for all methods grouped...\n")
    # Font sizes
    fs_num_1 = 35
    fs_num_2 = 25
    fs_num_3 = 25
    # Letters for annotation
    annotation_letters = list("abcdefghijklmnopqrstuvwxyz")
    # Create the figure to store the subplots.
    col_num = 6
    figure, axis = plt.subplots(1, col_num, figsize=(40, 20))
    row_fig_index = 0
    # Filter specific samples.
    for key_group in time_syn_sample_group_dict.keys():
        sample_id_group = time_syn_sample_group_dict[key_group]
        # Filtering the dataframe.
        total_time_all_filtered_df = total_time_all_df[total_time_all_df['sample'].isin(sample_id_group)]
        cur_axis = axis[row_fig_index]
        if total_time_all_filtered_df.empty:
            break
        # Plot
        sns.barplot(ax=cur_axis, data=total_time_all_filtered_df, x="sample", y="total_time", hue="method", palette="colorblind")
        # Replace the y axis ticks. Helps to rewrite ticks. For example the plot may continue on above the highest tick (tick with a value)
        # and thus has a bit more empty space without ticks. By rewriting the ticks a new tick will be provided to the top of the empty space.
        cur_y_ticks = cur_axis.get_yticks()
        cur_y_ticks_list = list(cur_y_ticks)
        cur_axis.set_yticks(cur_y_ticks_list)
        # Remove the legend, if for the last plot.
        if row_fig_index < col_num:
            cur_axis.legend().remove()
        # Labels
        x_axis_label = "Sample ID"
        y_axis_label = "Execution Time (min)"
        cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
        cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
        # Group label
        group_label = time_sample_group_label_dict[key_group]
        cur_axis.set_title(group_label, fontsize=fs_num_2, fontweight='bold')
        # Changing the size of the tick labels        
        cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
        # Letter
        cur_axis.annotate(annotation_letters[row_fig_index], xy=(0.95, 0.98), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
        # Place vertical lines to divie the sample barplots.
        sanme_num = len(sample_id_group)
        for i in range(1, sanme_num):
            cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
        # Row index
        row_fig_index += 1
    # Get the axis for the middle plot
    middle_axis = axis[3]
    # Title
    fig_title_label = "Execution Time vs Species Number"
    title_pos = middle_axis.title.get_position()
    figure.text(title_pos[0], title_pos[1] - 0.025, fig_title_label, ha='center', fontsize=fs_num_1, fontweight='bold')
    # Get the labels from the last axis with a plot.
    cur_axis = axis[row_fig_index-1]
    legend_mods, labels = cur_axis.get_legend_handles_labels()
    # Legend
    legend_title = "Taxonomy Method and Database"
    legend_props = FontProperties(weight='bold', size=fs_num_3)
    middle_axis.legend(handles=legend_mods, title=legend_title, loc='upper center', bbox_to_anchor=(0.0, -0.05), ncol=11, prop={'size': fs_num_3}, title_fontproperties=legend_props)
    # Padding from left, bottom, right, top
    figure.tight_layout(rect=[0, 0, 0.99, 0.97])
    # Space between the subplots.
    figure.subplots_adjust(wspace=0.3)
    # File paths
    full_time_grouped_path_png = "{}/full_time_method_grouped.png".format(time_dir)
    full_time_grouped_path_jpg = "{}/full_time_method_grouped.jpg".format(time_dir)
    # Saving the figure.
    figure.savefig(full_time_grouped_path_png)
    figure.savefig(full_time_grouped_path_jpg)
    # Closing the figure.
    plt.close(figure)


def plot_size_species(df_total_time_dict, sample_size_dict, methods_group, sample_speciesnum_dict, time_stats_dir_path, time_dir):
    # Font sizes
    fs_num_1 = 25
    fs_num_2 = 20
    fs_num_3 = 20
    # Letters for annotation
    annotation_letters = list("abcdefghijklmnopqrstuvwxyz")
    # Create a dictionary where the first key is the method and the second key the sample ID and its value the total time.
    method_sample_time_dict = {}
    for key_sample in df_total_time_dict.keys():
        temp_df = df_total_time_dict[key_sample]
        if not temp_df.empty:
            for mg in methods_group:
                if mg not in method_sample_time_dict.keys():
                    method_sample_time_dict[mg] = {}
                mg_df = temp_df.loc[temp_df['method'] == mg]
                if not mg_df.empty:
                    mg_total_time = mg_df['total_time'].iloc[0]
                    mg_total_time = float(mg_total_time)
                    mg_total_time = round(mg_total_time, 2)
                    method_sample_time_dict[mg][key_sample] = mg_total_time
    # Create a dataframe for each method.
    # Dataframe columns:
    # sample size total_time
    method_df_dict = {}
    for key_method in method_sample_time_dict.keys():
        mgtt_col_labels = ["sample", "size", "species_number", "total_time"]
        mgtt_col_df = pd.DataFrame(columns=mgtt_col_labels)
        row_index = 0
        for key_sample in method_sample_time_dict[key_method].keys():
            key_sample_int = int(key_sample)
            sample_size = sample_size_dict[key_sample_int]
            species_num = sample_speciesnum_dict[key_sample_int]
            total_time = method_sample_time_dict[key_method][key_sample]
            temp_list = [key_sample, sample_size, species_num, total_time]
            mgtt_col_df.loc[row_index] = temp_list
            row_index += 1
        method_df_dict[key_method] = mgtt_col_df
    
    # Create the figure to store the subplots.
    col_num = 5
    figure, axis = plt.subplots(1, col_num, figsize=(30, 10))
    row_fig_index = 0
    # Create the scatter plots for size vs time.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        if not plot_df.empty:
            # Save the dataframe in a file.
            size_time_tsv_path = "{}/size_time_{}.tsv".format(time_stats_dir_path, key_method)
            plot_df.to_csv(size_time_tsv_path, sep="\t", index=False)
            # The current axis.
            cur_axis = axis[row_fig_index]
            # Color
            cld_palette = sns.color_palette("colorblind")
            dot_color = cld_palette[0] 
            # Scatter plot
            sns.scatterplot(ax=cur_axis, data=plot_df, x="size", y="total_time", color=dot_color, s=100)
            # Regression plot
            # Convert the dataframe columns to lists of floats.
            x_points_pre = plot_df['size'].tolist()
            x_points = []
            for xpr in x_points_pre:
                xpr = float(xpr)
                x_points.append(xpr)
            y_points_pre = plot_df['total_time'].tolist()
            y_points = []
            for ypr in y_points_pre:
                ypr = float(ypr)
                y_points.append(ypr)
            # Convert the lists to numpy arrays.
            x_points_np = np.array(x_points)
            y_points_np = np.array(y_points)
            # Fit a line. It is a line because it set to fit a line of the first degree which is determiend by the "1" in the polyfit function.
            # This in turn returns two coefficients. To fit a second degree line "2" would be used. Same on for lines of a higher degree.
            a_coef, b_coef = np.polyfit(x_points_np, y_points_np, 1)
            cur_axis.plot(x_points_np, a_coef * x_points_np + b_coef)
            # Calculate the predicted y values
            y_points_pred = a_coef * x_points_np + b_coef
            # Calculate the R-squared value
            r_squared_value = r2_score(y_points_np, y_points_pred)
            # Rounding to two decimals.
            a_coef_r = round(a_coef, 2)
            b_coef_r = round(b_coef, 2)
            r_squared_value_r = round(r_squared_value, 2)
            # Write the coefficients and the R^2 values in a file.
            method_fit_path = "{}/method_{}_size_time_fit.txt".format(time_stats_dir_path, key_method)
            method_fit_file = open(method_fit_path, "w")
            method_fit_file.write("Method: {}\n".format(key_method))
            method_fit_file.write("y = {} * x + {}\n".format(a_coef, b_coef))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value))
            method_fit_file.write("Rounded:\n")
            method_fit_file.write("y = {} * x + {}\n".format(a_coef_r, b_coef_r))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value_r))
            method_fit_file.close()
            # Set the minimum value for the x and y axis.
            cur_axis.set_xlim(left=0)
            cur_axis.set_ylim(bottom=0)
            # Remove the 0.0 value from the x axis.
            x_axis_ticks = cur_axis.get_xticks()
            if x_axis_ticks[0] == 0.0:
                cur_axis.set_xticks(x_axis_ticks[1:])
            # Labels
            x_axis_label = "Sample Size (GB)"
            y_axis_label = "Execution Time (min)"
            if key_method == "k8":
                title_part = "Kraken2 db:8"
            elif key_method == "k16":
                title_part = "Kraken2 db:16"
            elif key_method == "k72":
                title_part = "Kraken2 db:72"
            elif key_method == "cnr":
                title_part = "COMEBin db:nr"
            elif key_method == "mnr":
                title_part = "MetaBinner db:nr"
            # Place labels
            cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_title(title_part, pad=20, loc='center', fontsize=fs_num_1, fontweight='bold')
            # Changing the size of the tick labels.
            cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
            # Letter
            cur_axis.annotate(annotation_letters[row_fig_index], xy=(-0.05, 1.05), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
            # Increament the index of the plot position.
            row_fig_index += 1
    # Get the axis for the middle plot.
    middle_axis = axis[2]
    # Title
    fig_title_label = "Execution Time vs Sample Size"
    title_pos = middle_axis.title.get_position()
    figure.text(title_pos[0], title_pos[1] - 0.04, fig_title_label, ha='center', fontsize=fs_num_1, fontweight='bold')
    # Layout
    figure.tight_layout(rect=[0, 0, 0.99, 0.95])
    # Space between the plots of the figure.
    figure.subplots_adjust(wspace=0.4)
    # File paths.
    size_time_png_path = "{}/size_time.png".format(time_dir)
    size_time_jpg_path = "{}/size_time.jpg".format(time_dir)
    # Saving the plots.
    figure.savefig(size_time_png_path)
    figure.savefig(size_time_jpg_path)
    # Closing the plot.
    plt.close(figure)
    
    # Create the figure to store the subplots.
    col_num = 5
    figure, axis = plt.subplots(1, col_num, figsize=(30, 10))
    row_fig_index = 0
    # Create the scatter plots for species number vs time.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        if not plot_df.empty:
            # Save the dataframe in a file.
            species_time_tsv_path = "{}/species_time_{}.tsv".format(time_stats_dir_path, key_method)
            plot_df.to_csv(species_time_tsv_path, sep="\t", index=False)
            # The current axis.
            cur_axis = axis[row_fig_index]
            # Returns a list of colors or continuous colormap defining a palette.
            species_num = len(plot_df['species_number'].unique())
            cld_palette = sns.color_palette("colorblind", n_colors=species_num)
            # Plot
            sns.scatterplot(ax=cur_axis, data=plot_df, x="species_number", y="total_time", hue="species_number", palette=cld_palette, s=100)
            # Regression plot
            # Convert the dataframe columns to lists of floats.
            x_points_pre = plot_df['species_number'].tolist()
            x_points = []
            for xpr in x_points_pre:
                xpr = float(xpr)
                x_points.append(xpr)
            y_points_pre = plot_df['total_time'].tolist()
            y_points = []
            for ypr in y_points_pre:
                ypr = float(ypr)
                y_points.append(ypr)
            # Convert the lists to numpy arrays.
            x_points_np = np.array(x_points)
            y_points_np = np.array(y_points)
            # Fit a line. It is a line because it set to fit a line of the first degree which is determiend by the "1" in the polyfit function.
            # This in turn returns two coefficients. To fit a second degree line "2" would be used. Same on for lines of a higher degree.
            a_coef, b_coef = np.polyfit(x_points_np, y_points_np, 1)
            cur_axis.plot(x_points_np, a_coef * x_points_np + b_coef)
            # Calculate the predicted y values
            y_points_pred = a_coef * x_points_np + b_coef
            # Calculate the R-squared value
            r_squared_value = r2_score(y_points_np, y_points_pred)
            # Rounding to two decimals.
            a_coef_r = round(a_coef, 2)
            b_coef_r = round(b_coef, 2)
            r_squared_value_r = round(r_squared_value, 2)
            # Write the coefficients and the R^2 values in a file.
            method_fit_path = "{}/method_{}_species_time_fit.txt".format(time_stats_dir_path, key_method)
            method_fit_file = open(method_fit_path, "w")
            method_fit_file.write("Method: {}\n".format(key_method))
            method_fit_file.write("y = {} * x + {}\n".format(a_coef, b_coef))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value))
            method_fit_file.write("Rounded:\n")
            method_fit_file.write("y = {} * x + {}\n".format(a_coef_r, b_coef_r))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value_r))
            method_fit_file.close()
            # Set the minimum value for the x and y axis.
            cur_axis.set_xlim(left=0)
            cur_axis.set_ylim(bottom=0)
            # Remove the 0.0 value from the x axis.
            x_axis_ticks = cur_axis.get_xticks()
            if x_axis_ticks[0] == 0.0:
                cur_axis.set_xticks(x_axis_ticks[1:])
            # Labels
            x_axis_label = "Species Number"
            y_axis_label = "Execution Time (min)"
            if key_method == "k8":
                title_part = "Kraken2 db:8"
            elif key_method == "k16":
                title_part = "Kraken2 db:16"
            elif key_method == "k72":
                title_part = "Kraken2 db:72"
            elif key_method == "cnr":
                title_part = "COMEBin db:nr"
            elif key_method == "mnr":
                title_part = "MetaBinner db:nr"
            # Place labels
            cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_title(title_part, pad=20, loc='center', fontsize=fs_num_1, fontweight='bold')
            # Remove the legend, if for the last plot.
            if row_fig_index < col_num:
                cur_axis.legend().remove()
            # Changing the size of the tick labels.
            cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
            # Letter
            cur_axis.annotate(annotation_letters[row_fig_index], xy=(-0.05, 1.05), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
            # Increament the index of the plot position.
            row_fig_index += 1
    # Get the axis for the middle plot.
    middle_axis = axis[2]
    # Title
    fig_title_label = "Execution Time vs Species Number"
    title_pos = middle_axis.title.get_position()
    figure.text(title_pos[0], title_pos[1] - 0.04, fig_title_label, ha='center', fontsize=fs_num_1, fontweight='bold')
    # Legend
    legend_title = "Species Number"
    legend_props = FontProperties(weight='bold', size=fs_num_3)
    middle_axis.legend(title=legend_title, loc='upper center', markerscale=2, bbox_to_anchor=(0.5, -0.2), ncol=11, prop={'size': fs_num_3}, title_fontproperties=legend_props)
    # Layout
    figure.tight_layout(rect=[0, 0, 0.99, 0.95])
    # Space between the plots of the figure.
    figure.subplots_adjust(wspace=0.4)
    # File paths.
    size_time_png_path = "{}/species_time.png".format(time_dir)
    size_time_jpg_path = "{}/species_time.jpg".format(time_dir)
    # Saving the plots.
    figure.savefig(size_time_png_path)
    figure.savefig(size_time_jpg_path)
    # Closing the plot.
    plt.close(figure)

    # Create the figure to store the subplots.
    col_num = 5
    figure, axis = plt.subplots(1, col_num, figsize=(30, 10))
    row_fig_index = 0
    # Species that are 10 in number are divided in the following groups:
    group_10_synthetic = ["8", "18", "19"]
    group_10_cultured = ["2", "9", "16", "17"]
    # Labels
    label_dict = {
        "10": '10',
        "40": '40',
        "120": '120',
        "500": '500',
        "1000": '1000'
    }
    # Label keys
    label_keys_list = list(label_dict.keys())
    # Create the scatter plots for species number vs mean time for each group of species number.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        if not plot_df.empty:
            grouped_plot_df = plot_df.groupby('species_number')
            ur_mean_plot_df = grouped_plot_df['total_time'].mean()
            mean_plot_df = ur_mean_plot_df.reset_index()
            # Add the labels to the DataFrame.
            label_list = []
            for i in range(len(mean_plot_df)):
                species_num = mean_plot_df.iloc[i]['species_number']
                species_num_int = int(species_num)
                species_num_str = str(species_num_int)
                if species_num_str in label_keys_list:
                    cur_label = label_dict[species_num_str]
                    label_list.append(cur_label)
            mean_plot_df['label'] = label_list
            # Filter DataFrame for each subgroup.
            # Compute the mean total_time for each species_number in each subgroup.
            # 10 Species - Cultured
            cul_10_df = plot_df[plot_df['sample'].isin(group_10_cultured)]
            if not cul_10_df.empty:
                grouped_10_cul_df = cul_10_df.groupby('species_number')
                ur_mean_10_cul_df = grouped_10_cul_df['total_time'].mean()
                mean_10_cul_df = ur_mean_10_cul_df.reset_index()
                mean_10_cul_df['label'] = "10 - Cultured"
                # Add the dataframe to the existing one for the mean time values.
                mean_plot_df = pd.concat([mean_plot_df, mean_10_cul_df], ignore_index=True)
                # Move the row to the second place.
                last_row = mean_plot_df.iloc[-1:]
                part_mean_plot_df = mean_plot_df.iloc[:-1]
                ra_mean_plot_df = pd.concat([part_mean_plot_df.iloc[:1], last_row, part_mean_plot_df.iloc[1:]])
                mean_plot_df = ra_mean_plot_df.reset_index(drop=True)
            # 10 Species - Simulated
            syn_10_df = plot_df[plot_df['sample'].isin(group_10_synthetic)]
            if not syn_10_df.empty:
                grouped_syn_10_df = syn_10_df.groupby('species_number')
                ur_mean_10_syn_df = grouped_syn_10_df['total_time'].mean()
                mean_10_syn_df = ur_mean_10_syn_df.reset_index()
                mean_10_syn_df['label'] = "10 - Simulated"
                # Add the dataframe to the existing one for the mean time values.
                mean_plot_df = pd.concat([mean_plot_df, mean_10_syn_df], ignore_index=True)
                # Move the row to the second place.
                last_row = mean_plot_df.iloc[-1:]
                part_mean_plot_df = mean_plot_df.iloc[:-1]
                ra_mean_plot_df = pd.concat([part_mean_plot_df.iloc[:1], last_row, part_mean_plot_df.iloc[1:]])
                mean_plot_df = ra_mean_plot_df.reset_index(drop=True)
            # Rename the column from "total_time" to "mean_time".
            mean_plot_df = mean_plot_df.rename(columns={'total_time': 'mean_time'})
            # Rounding the values for the mean time.
            mean_plot_df['mean_time'] = mean_plot_df['mean_time'].round(2)
            # Save the dataframe in a file.
            species_meantime_tsv_path = "{}/species_meantime_{}.tsv".format(time_stats_dir_path, key_method)
            mean_plot_df.to_csv(species_meantime_tsv_path, sep="\t", index=False)
            # The current axis.
            cur_axis = axis[row_fig_index]
            # Returns a list of colors or continuous colormap defining a palette.
            labels_num = len(mean_plot_df['label'].unique())
            cld_palette = sns.color_palette("colorblind", n_colors=labels_num)
            # Scatter plot
            sns.scatterplot(ax=cur_axis, data=mean_plot_df, x="species_number", y="mean_time", hue="label", palette=cld_palette, s=100)
            # Regression plot
            # Convert the dataframe columns to lists of floats.
            # The data are extracted from the dataframe which was generated by computing the mean values of the execution times for each category of samples
            # of species number.
            x_points_pre = mean_plot_df['species_number'].tolist()
            x_points = []
            for xpr in x_points_pre:
                xpr = float(xpr)
                x_points.append(xpr)
            y_points_pre = mean_plot_df['mean_time'].tolist()
            y_points = []
            for ypr in y_points_pre:
                ypr = float(ypr)
                y_points.append(ypr)
            # Convert the lists to numpy arrays.
            x_points_np = np.array(x_points)
            y_points_np = np.array(y_points)
            # Fit a line. It is a line because it set to fit a line of the first degree which is determiend by the "1" in the polyfit function.
            # This in turn returns two coefficients. To fit a second degree line "2" would be used. Same on for lines of a higher degree.
            a_coef, b_coef = np.polyfit(x_points_np, y_points_np, 1)
            cur_axis.plot(x_points_np, a_coef * x_points_np + b_coef)
            # Calculate the predicted y values
            y_points_pred = a_coef * x_points_np + b_coef
            # Calculate the R-squared value
            r_squared_value = r2_score(y_points_np, y_points_pred)
            # Rounding to two decimals.
            a_coef_r = round(a_coef, 2)
            b_coef_r = round(b_coef, 2)
            r_squared_value_r = round(r_squared_value, 2)
            # Write the coefficients and the R^2 values in a file.
            method_fit_path = "{}/method_{}_species_mean_time_fit.txt".format(time_stats_dir_path, key_method)
            method_fit_file = open(method_fit_path, "w")
            method_fit_file.write("Method: {}\n".format(key_method))
            method_fit_file.write("y = {} * x + {}\n".format(a_coef, b_coef))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value))
            method_fit_file.write("Rounded:\n")
            method_fit_file.write("y = {} * x + {}\n".format(a_coef_r, b_coef_r))
            method_fit_file.write("R^2 = {}\n".format(r_squared_value_r))
            method_fit_file.close()
            # Set the minimum value for the x and y axis.
            cur_axis.set_xlim(left=0)
            cur_axis.set_ylim(bottom=0)
            # Remove the 0.0 value from the x axis.
            x_axis_ticks = cur_axis.get_xticks()
            if x_axis_ticks[0] == 0.0:
                cur_axis.set_xticks(x_axis_ticks[1:])
            # Labels
            x_axis_label = "Species Number"
            y_axis_label = "Mean Execution Time (min)"
            if key_method == "k8":
                title_part = "Kraken2 db:8"
            elif key_method == "k16":
                title_part = "Kraken2 db:16"
            elif key_method == "k72":
                title_part = "Kraken2 db:72"
            elif key_method == "cnr":
                title_part = "COMEBin db:nr"
            elif key_method == "mnr":
                title_part = "MetaBinner db:nr"
            # Place labels
            cur_axis.set_xlabel(x_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_ylabel(y_axis_label, fontsize=fs_num_2, fontweight='bold')
            cur_axis.set_title(title_part, pad=20, loc='center', fontsize=fs_num_1, fontweight='bold')
            # Remove the legend, if for the last plot.
            if row_fig_index < col_num:
                cur_axis.legend().remove()
            # Changing the size of the tick labels.
            cur_axis.tick_params(axis='both', which='major', labelsize=fs_num_3)
            # Letter
            cur_axis.annotate(annotation_letters[row_fig_index], xy=(-0.05, 1.05), xycoords='axes fraction', fontsize=fs_num_1, fontweight='bold', ha='center', va='center')
            # Increament the index of the plot position.
            row_fig_index += 1
    # Get the axis for the middle plot
    middle_axis = axis[2]
    # Title
    fig_title_label = "Mean Execution Time vs Species Number"
    title_pos = middle_axis.title.get_position()
    figure.text(title_pos[0], title_pos[1] - 0.04, fig_title_label, ha='center', fontsize=fs_num_1, fontweight='bold')
    # Legend
    legend_title = "Species Number"
    legend_props = FontProperties(weight='bold', size=fs_num_3)
    middle_axis.legend(title=legend_title, loc='upper center', markerscale=2, bbox_to_anchor=(0.5, -0.2), ncol=11, prop={'size': fs_num_3}, title_fontproperties=legend_props)
    # Layout
    figure.tight_layout(rect=[0, 0, 0.99, 0.95])
    # Space between the plots of the figure.
    figure.subplots_adjust(wspace=0.4)
    # File paths
    size_time_png_path = "{}/species_meantime.png".format(time_dir)
    size_time_jpg_path = "{}/species_meantime.jpg".format(time_dir)
    # Saving the plots.
    figure.savefig(size_time_png_path)
    figure.savefig(size_time_jpg_path)
    # Closing the plot.
    plt.close(figure)


def time_analysis(time_dir, time_stats_dir_path, time_syn_sample_group_dict, time_sample_group_label_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict):
    print("Performing analysis of the execution time...\n")
    # Collect the time periods
    time_dict = collect_times(time_stats_dir_path)

    # Seperate the time stages
    time_full_dict, time_tool_dict = time_stages(time_dict)

    # Fill the full dict with the missing values.
    time_full_dict = fill_times(time_dict, time_full_dict)

    # Create the dataframe for the complete time slots.
    df_full_time_dict, df_full_time_all = create_df_stacked(time_full_dict)
    
    # Create the dataframe for the total time.
    df_total_time_dict, total_time_all_df = create_df_total_time(time_tool_dict)

    # Count the time up to the last common stage for each sample.
    common_time_dict = count_common_time(df_full_time_dict)

    # Plot the dataframe
    plot_full_sample_time(df_full_time_dict, common_time_dict, time_dir, time_stats_dir_path)
    plot_total_sample_time(df_total_time_dict, common_time_dict, time_dir, time_stats_dir_path)
    plot_total_all_time(total_time_all_df, time_dir, time_stats_dir_path)
    plot_full_group_sample_time(df_full_time_all, time_dir, time_stats_dir_path, time_syn_sample_group_dict, time_sample_group_label_dict, methods_time_group)
    plot_full_group_method_time(total_time_all_df, time_dir, time_syn_sample_group_dict, time_sample_group_label_dict, time_stats_dir_path)

    # Plot total time vs size
    plot_size_species(df_total_time_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict, time_stats_dir_path, time_dir)


def benchstats(benchmark_path="12864_2022_8803_MOESM1_ESM.txt", ps_results="", ps_output_dir="analysis_results", methods_group="", time_status=True):
    # Check for the input files.
    if (not benchmark_path) or (not ps_results):
        print("Input file not set. Exiting.")
        exit()
    if ((not os.path.exists(benchmark_path)) or (not os.path.exists(ps_results))):
        print("Input file not found. Exiting.")
        exit()

    # Create the output directory, overwriting it if already present.
    if os.path.exists(ps_output_dir):
        shutil.rmtree(ps_output_dir)
    os.mkdir(ps_output_dir)

    # File paths
    benchmark_mod_path = "{}/12864_2022_8803_MOESM1_ESM_mod.txt".format(ps_output_dir)
    benchmark_mod_spabu_path = "{}/12864_2022_8803_MOESM1_ESM_mod_spabu.txt".format(ps_output_dir)

    # Directory paths
    ps_dir = "{}/ps_analysis".format(ps_output_dir)
    time_dir = "{}/time_info".format(ps_output_dir)
    plot_dir_parh = "{}/plots".format(ps_output_dir)
    sp_dir_path = "{}/sp_info".format(ps_output_dir)
    stats_dir_path = "{}/stats".format(ps_output_dir)
    time_stats_dir_path = "{}/time_stats".format(stats_dir_path)
    # Directories for the results
    ps_kraken_8 = "{}/kraken_8".format(ps_dir)
    ps_kraken_16 = "{}/kraken_16".format(ps_dir)
    ps_kraken_72 = "{}/kraken_72".format(ps_dir)
    # Add the filtered results for Kraken.
    ps_mtbr_nr = "{}/metabinner_nr".format(ps_dir)
    ps_cmbn_nr = "{}/comebin_nr".format(ps_dir)
    # Path dictinoary
    analysis_results_dict = {
        "k8": "{}/kraken_8".format(ps_dir),
        "k16": "{}/kraken_16".format(ps_dir),
        "k72": "{}/kraken_72".format(ps_dir),
        "cnr": "{}/comebin_nr".format(ps_dir),
        "mnr": "{}/metabinner_nr".format(ps_dir),
    }

    # Creating the directory which will contain the species and percentagies for each method and sample..
    if os.path.exists(sp_dir_path):
        shutil.rmtree(sp_dir_path)
    os.mkdir(sp_dir_path)

    # Craeting the statistics directory.
    if os.path.exists(stats_dir_path):
        shutil.rmtree(stats_dir_path)
    os.mkdir(stats_dir_path)

    # Creating the time-statistics directory.
    if os.path.exists(time_stats_dir_path):
        shutil.rmtree(time_stats_dir_path)
    os.mkdir(time_stats_dir_path)

    # Analyze the bencharking file.
    if benchmark_path == "":
        print("Error. No input file was found. Exiting.\n")
        exit()

    # Determine the methods for which to get statistics and plots
    methods_group = methods_group.split(",")

    # Setting the grid lines.
    sns.set_theme(style="whitegrid")
    
    # Analyzing the benchmark dataset.
    br_info_dict = br_analysis(benchmark_path, benchmark_mod_path, benchmark_mod_spabu_path)
    
    if not br_info_dict:
        print("Error. No information from file 1. Exiting.")
        exit()

    # Copy the files with the results in the proper directories and rename them.
    crfiles(ps_dir, ps_results, analysis_results_dict, time_dir, time_stats_dir_path)

    # Collect the kraken filters.
    process_kraken_filter_files(ps_dir, stats_dir_path)

    # Collect the information from the results of ProteoSeeker
    # COMEBin/MetaBinner results
    # binid     species_name                taxid   rank        lineage                                                                                                                                                                                    percentage_input_reads      percentage_preprocessed_reads
    # 0	        Cupriavidus	                106589	genus	    cellular organisms;Bacteria;Pseudomonadota;Betaproteobacteria;Burkholderiales;Burkholderiaceae;Cupriavidus	                                                                               1.48	                       1.48
    # 1     	Micromonospora krabiensis	307121	species 	cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora;Micromonospora krabiensis	                               1.45	                       1.46
    # 1	        Micromonospora	            1873	genus	    cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora	                                                       1.45	                       1.46
    # 1	        Nonomuraea sp. TT08I-71	    2733864	species	    cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Streptosporangiales;Streptosporangiaceae;Nonomuraea;unclassified Nonomuraea;Nonomuraea sp. TT08I-71	       1.45	                       1.46
    # 1	        Micromonospora sp. WMMA2032	2039870	species	    cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora;unclassified Micromonospora;Micromonospora sp. WMMA2032   1.45	                       1.46
    # 1	        Micromonospora sp. AMSO31t	2650566	species	    cellular organisms;Bacteria;Terrabacteria group;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Micromonospora;unclassified Micromonospora;Micromonospora sp. AMSO31t	   1.45	                       1.46
    # Notes:
    # 1. Only the taxonomy ranks of "species" are taken into account.
    # 2. More than one taxa may be associated with one bin.
    # 3. The same taxa may be associated with more than one bins.
    # 4. The taxa identified for each bin have the same percentage. The percentage is computed for the bin.
    # 5. The percentages of all unique bins do not necessarily add up to 100%. The might be bins which were not associated with a species. 
    # Each bin corresponds to specific contigs. Each contig in turn corresponds to a set of reads. Hence, each bin is associated with a set
    # of read. If the bin is not associated with a species, the corresponding set of reads will not be associated with a species.
    # 6. The end result of the process is meant to be the associated of each protein with a species. Whether one protein will be associated
    # with more than one species and whether the sum of the percentages of all the species will be 100% is not of significance for this process
    # and neither does it affect negatively the association process of each protein with a species.
    # Questions:
    # 1. How many taxa will be taken into consideration from each bin?
    # 2. Which percentage will be used for a taxa associated with more than one bin?
    # Handling:
    # 1. All taxa are taken into consideration from each bin.
    # 2. The percentage of a taxon is computed as the sum of the percentages of each bin associated with that taxon.
    # Observations:
    # Allowing more than one taxa to be taken into account from each bin means that the total percentage of the species will be higher than 100%.
    # It also means that metrics which take into account the percentages of the predicted species, such as the L1 norm should be based on the species which are sole predictions for one bin ID each.
    # Kraken2
    filter_name_k_base = "_kraken_species.tsv"
    spec_label_k8 = "k8"
    spec_label_k16 = "k16"
    spec_label_k72 = "k72"
    spec_label_k8_ng = "k8_ng"
    spec_label_k16_ng = "k16_ng"
    spec_label_k72_ng = "k72_ng"
    spec_label_k8_0c1 = "k8_0c1"
    spec_label_k16_0c1 = "k16_0c1"
    spec_label_k72_0c1 = "k72_0c1"
    spec_label_k8_1c0 = "k8_1c0"
    spec_label_k16_1c0 = "k16_1c0"
    spec_label_k72_1c0 = "k72_1c0"
    spec_label_k8_10 = "k8_10"
    spec_label_k16_10 = "k16_10"
    spec_label_k72_10 = "k72_10"
    spec_label_k8_100 = "k8_100"
    spec_label_k16_100 = "k16_100"
    spec_label_k72_100 = "k72_100"
    spec_label_cnr = "cnr"
    spec_label_mnr = "mnr"
    if spec_label_k8 in methods_group:
        kraken_8_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_base, sp_dir_path, spec_label_k8)
    if spec_label_k16 in methods_group:
        kraken_16_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_base, sp_dir_path, spec_label_k16)
    if spec_label_k72 in methods_group:
        kraken_72_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_base, sp_dir_path, spec_label_k72)
    # Kraken2: non-gut
    filter_name_k_f_base = "thr_-1_"
    if spec_label_k8_ng in methods_group:
        kraken_8_ng_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base, sp_dir_path, spec_label_k8_ng)
    if spec_label_k16_ng in methods_group:
        kraken_16_ng_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base, sp_dir_path, spec_label_k16_ng)
    if spec_label_k72_ng in methods_group:
        kraken_72_ng_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base, sp_dir_path, spec_label_k72_ng)
    # Kraken2: 0.1%
    filter_name_k_f_base = "thr_0.1_"
    if spec_label_k8_0c1 in methods_group:
        kraken_8_0c1_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base, sp_dir_path, spec_label_k8_0c1)
    if spec_label_k16_0c1 in methods_group:
        kraken_16_0c1_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base, sp_dir_path, spec_label_k16_0c1)
    if spec_label_k72_0c1 in methods_group:
        kraken_72_0c1_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base, sp_dir_path, spec_label_k72_0c1)
    # Kraken2: 1%
    filter_name_k_f_base = "thr_1.0_"
    if spec_label_k8_1c0 in methods_group:
        kraken_8_1c0_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base, sp_dir_path, spec_label_k8_1c0)
    if spec_label_k16_1c0 in methods_group:
        kraken_16_1c0_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base, sp_dir_path, spec_label_k16_1c0)
    if spec_label_k72_1c0 in methods_group:
        kraken_72_1c0_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base, sp_dir_path, spec_label_k72_1c0)
    # Kraken2: 10
    filter_name_k_f_base = "thr_10_"
    if spec_label_k8_10 in methods_group:
        kraken_8_10_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base, sp_dir_path, spec_label_k8_10)
    if spec_label_k16_10 in methods_group:
        kraken_16_10_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base, sp_dir_path, spec_label_k16_10)
    if spec_label_k72_10 in methods_group:
        kraken_72_10_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base, sp_dir_path, spec_label_k72_10)
    # Kraken2: 100
    filter_name_k_f_base = "thr_100_"
    if spec_label_k8_100 in methods_group:
        kraken_8_100_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base, sp_dir_path, spec_label_k8_100)
    if spec_label_k16_100 in methods_group:
        kraken_16_100_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base, sp_dir_path, spec_label_k16_100)
    if spec_label_k72_100 in methods_group:
        kraken_72_100_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base, sp_dir_path, spec_label_k72_100)
    # COMEBin nr
    filter_name_c = "_b_summary_info_comebin.tsv"
    if spec_label_cnr in methods_group:
        cmbn_nr_info_dict, cmbn_nr_sole_info_dict, cmbn_nr_perc_sum_dict = ps_comebin_analyze(ps_cmbn_nr, filter_name_c, sp_dir_path, spec_label_cnr)
    # MetaBinner nr
    filter_name_m =  "_b_summary_info_metabinner.tsv"
    if spec_label_mnr in methods_group:
        mtbr_nr_info_dict, mtbr_nr_sole_info_dict, mtbr_nr_perc_sum_dict = ps_metabinner_analyze(ps_mtbr_nr, filter_name_m, sp_dir_path, spec_label_mnr)

    kraken_8_stats_dict = {}
    kraken_16_stats_dict = {}
    kraken_72_stats_dict = {}
    kraken_8_0c1_stats_dict = {}
    kraken_16_0c1_stats_dict = {}
    kraken_72_0c1_stats_dict = {}
    kraken_8_1c0_stats_dict = {}
    kraken_16_1c0_stats_dict = {}
    kraken_72_1c0_stats_dict = {}
    kraken_8_10_stats_dict = {}
    kraken_16_10_stats_dict = {}
    kraken_72_10_stats_dict = {}
    kraken_8_100_stats_dict = {}
    kraken_16_100_stats_dict = {}
    kraken_72_100_stats_dict = {}
    kraken_8_ng_stats_dict = {}
    kraken_16_ng_stats_dict = {}
    kraken_72_ng_stats_dict = {}
    mtbr_nr_stats_dict = {}
    cmbn_nr_stats_dict = {}
    # Compare the results of ProteoSeeker with the benchmark dataset.
    # Kraken2
    if "k8" in methods_group:
        label = "kraken_8"
        kraken_8_stats_dict = comp_stats(br_info_dict, kraken_8_info_dict, None, stats_dir_path, label)
    if "k16" in methods_group:
        label = "kraken_16"
        kraken_16_stats_dict = comp_stats(br_info_dict, kraken_16_info_dict, None, stats_dir_path, label)
    if "k72" in methods_group:
        label = "kraken_72"
        kraken_72_stats_dict = comp_stats(br_info_dict, kraken_72_info_dict, None, stats_dir_path, label)
    # Kraken2 0.1%
    if "k8_0c1" in methods_group:
        label = "kraken_8_0.1%"
        kraken_8_0c1_stats_dict = comp_stats(br_info_dict, kraken_8_0c1_info_dict, None, stats_dir_path, label)
    if "k16_0c1" in methods_group:
        label = "kraken_16_0.1%"
        kraken_16_0c1_stats_dict = comp_stats(br_info_dict, kraken_16_0c1_info_dict, None, stats_dir_path, label)
    if "k72_0c1" in methods_group:
        label = "kraken_72_0.1%"
        kraken_72_0c1_stats_dict = comp_stats(br_info_dict, kraken_72_0c1_info_dict, None, stats_dir_path, label)
    # Kraken2 1.0%
    if "k8_1c0" in methods_group:
        label = "kraken_8_1.0%"
        kraken_8_1c0_stats_dict = comp_stats(br_info_dict, kraken_8_1c0_info_dict, None, stats_dir_path, label)
    if "k16_1c0" in methods_group:
        label = "kraken_16_1.0%"
        kraken_16_1c0_stats_dict = comp_stats(br_info_dict, kraken_16_1c0_info_dict, None, stats_dir_path, label)
    if "k72_1c0" in methods_group:
        label = "kraken_72_1.0%"
        kraken_72_1c0_stats_dict = comp_stats(br_info_dict, kraken_72_1c0_info_dict, None, stats_dir_path, label)
    # Kraken2 10
    if "k8_10" in methods_group:
        label = "kraken_8_10"
        kraken_8_10_stats_dict = comp_stats(br_info_dict, kraken_8_10_info_dict, None, stats_dir_path, label)
    if "k16_10" in methods_group:
        label = "kraken_16_10"
        kraken_16_10_stats_dict = comp_stats(br_info_dict, kraken_16_10_info_dict, None, stats_dir_path, label)
    if "k72_10" in methods_group:
        label = "kraken_72_10"
        kraken_72_10_stats_dict = comp_stats(br_info_dict, kraken_72_10_info_dict, None, stats_dir_path, label)
    # Kraken2 100
    if "k8_100" in methods_group:
        label = "kraken_8_100"
        kraken_8_100_stats_dict = comp_stats(br_info_dict, kraken_8_100_info_dict, None, stats_dir_path, label)
    if "k16_100" in methods_group:
        label = "kraken_16_100"
        kraken_16_100_stats_dict = comp_stats(br_info_dict, kraken_16_100_info_dict, None, stats_dir_path, label)
    if "k72_100" in methods_group:
        label = "kraken_72_100"
        kraken_72_100_stats_dict = comp_stats(br_info_dict, kraken_72_100_info_dict, None, stats_dir_path, label)
    # Kraken2 non-gut
    if "k8_ng" in methods_group:
        label = "kraken_8_ng"
        kraken_8_ng_stats_dict = comp_stats(br_info_dict, kraken_8_ng_info_dict, None, stats_dir_path, label)
    if "k16_ng" in methods_group:
        label = "kraken_16_ng"
        kraken_16_ng_stats_dict = comp_stats(br_info_dict, kraken_16_ng_info_dict, None, stats_dir_path, label)
    if "k72_ng" in methods_group:
        label = "kraken_72_ng"
        kraken_72_ng_stats_dict = comp_stats(br_info_dict, kraken_72_ng_info_dict, None, stats_dir_path, label)
    # Kraken2 1%
    # COMEBin
    if "cnr" in methods_group:
        label = "comebin_nr"
        cmbn_nr_stats_dict = comp_stats(br_info_dict, cmbn_nr_info_dict, cmbn_nr_sole_info_dict, stats_dir_path, label)
    # MetaBinner
    if "mnr" in methods_group:
        label = "metabinner_nr"
        mtbr_nr_stats_dict = comp_stats(br_info_dict, mtbr_nr_info_dict, mtbr_nr_sole_info_dict, stats_dir_path, label)
    print("The computation of the statistics was completed.\n")

    # Combine the information from all the dictionaries.
    comb_info_dict = {}
    # Kraken2 db:8 - Filters: No filter, 0.1%, 1.0%, 10, 100, non-gut
    if kraken_8_stats_dict:
        comb_info_dict["Kraken2 db:8"] = kraken_8_stats_dict
    if kraken_8_0c1_stats_dict:
        comb_info_dict["Kraken2 db:8 0.1%"] = kraken_8_0c1_stats_dict
    if kraken_8_1c0_stats_dict:
        comb_info_dict["Kraken2 db:8 1.0%"] = kraken_8_1c0_stats_dict
    if kraken_8_10_stats_dict:
        comb_info_dict["Kraken2 db:8 10"] = kraken_8_10_stats_dict
    if kraken_8_100_stats_dict:
        comb_info_dict["Kraken2 db:8 100"] = kraken_8_100_stats_dict
    if kraken_8_ng_stats_dict:
        comb_info_dict["Kraken2 db:8 non-gut"] = kraken_8_ng_stats_dict
    # Kraken2 db:16 - Filters: No filter, 0.1%, 1.0%, 10, 100, non-gut
    if kraken_16_stats_dict:
        comb_info_dict["Kraken2 db:16"] = kraken_16_stats_dict
    if kraken_16_0c1_stats_dict:
        comb_info_dict["Kraken2 db:16 0.1%"] = kraken_16_0c1_stats_dict
    if kraken_16_1c0_stats_dict:
        comb_info_dict["Kraken2 db:16 1.0%"] = kraken_16_1c0_stats_dict
    if kraken_16_10_stats_dict:
        comb_info_dict["Kraken2 db:16 10"] = kraken_16_10_stats_dict
    if kraken_16_100_stats_dict:
        comb_info_dict["Kraken2 db:16 100"] = kraken_16_100_stats_dict
    if kraken_16_ng_stats_dict:
        comb_info_dict["Kraken2 db:16 non-gut"] = kraken_16_ng_stats_dict
    # Kraken db:72 - Filters: No filter, 0.1%, 1.0%, 10, 100, non-gut
    if kraken_72_stats_dict:
        comb_info_dict["Kraken2 db:72"] = kraken_72_stats_dict
    if kraken_72_0c1_stats_dict:
        comb_info_dict["Kraken2 db:72 0.1%"] = kraken_72_0c1_stats_dict
    if kraken_72_1c0_stats_dict:
        comb_info_dict["Kraken2 db:72 1.0%"] = kraken_72_1c0_stats_dict
    if kraken_72_10_stats_dict:
        comb_info_dict["Kraken2 db:72 10"] = kraken_72_10_stats_dict
    if kraken_72_100_stats_dict:
        comb_info_dict["Kraken2 db:72 100"] = kraken_72_100_stats_dict
    if kraken_72_ng_stats_dict:
        comb_info_dict["Kraken2 db:72 non-gut"] = kraken_72_ng_stats_dict
    # COMEBin nr
    if cmbn_nr_stats_dict:
        comb_info_dict["COMEBin db:nr"] = cmbn_nr_stats_dict
    # MetaBinner nr
    if mtbr_nr_stats_dict:
        comb_info_dict["MetaBinner db:nr"] = mtbr_nr_stats_dict

    # Create plots.
    # Metric groups.
    metric_group_1 = ["true_positive", "false_positive", "false_negative"]
    metric_group_2 = ["sensitivity", "precision", "accuracy"]
    metric_group_3 = ["f1_score", "jaccard_index", "l1_norm"]
    metric_group_dict = {
        1: metric_group_1,
        2: metric_group_2,
        3: metric_group_3
    }
    metric_label_dict_1 = {
        "gold_species_number": "Abundance of Gold Standard Species",
        "predicted_species_number": "Abundance of Predicted Species",
        "common_species": "Abundance of Common Species",
        "unique_species_total": "Abundance of Gold Standard and Predicted Species",
        "unique_gold": "Abundance of Species Unique to the Gold Standard Group",
        "unique_predicted": "Abundance of Species Unique to the Predicted Group",
        "true_positive": "True Positive",
        "false_positive": "False Positive",
        "false_negative": "False Negative",
        "sensitivity": "Sensitivity",
        "precision": "Precision",
        "accuracy": "Accuracy",
        "f1_score": "F1 Score",
        "jaccard_index": "Jaccard Index",
        "l1_norm": "L1 Norm"
    }
    metric_label_dict_2 = {
        "gold_species_number": "Abundance of Gold Standard Species",
        "predicted_species_number": "Abundance of Predicted Species",
        "common_species": "Abundance of Common Species",
        "unique_species_total": "Abundance of Gold Standard and Predicted Species",
        "unique_gold": "Abundance of Species Unique to the Gold Standard Group",
        "unique_predicted": "Abundance of Species Unique to the Predicted Group",
        "true_positive": "True Positive (Abundance of Species Common to Both the Gold Standard and Predicted Groups)",
        "false_positive": "False Positive (Abundance of Species Unique to the Predicted Group)",
        "false_negative": "False Negative (Abundance of Species Unique to the Gold Standard Group)",
        "sensitivity": "Sensitivity",
        "precision": "Precision",
        "accuracy": "Accuracy",
        "f1_score": "F1 Score",
        "jaccard_index": "Jaccard Index",
        "l1_norm": "L1 Norm"
    }

    # Plots.
    pandas_dict, pandas_sab_dict = design_grouped_plots(comb_info_dict, metric_label_dict_1, plot_dir_parh, stats_dir_path)
    sab_status = False
    design_metric_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict_1, plot_dir_parh, sab_status)
    sab_status = True
    design_metric_grouped_plots(pandas_sab_dict, metric_group_dict, metric_label_dict_1, plot_dir_parh, sab_status)

    # Analyze the statistics.
    rank_stats(pandas_dict, stats_dir_path, metric_label_dict_1, metric_label_dict_2)

    if time_status:
        time_syn_sample_groups_1 = [8, 18, 19]
        time_syn_sample_groups_2 = [2, 9, 16, 17]
        time_syn_sample_groups_3 = [5, 6, 11]
        time_syn_sample_groups_4 = [7, 13, 15]
        time_syn_sample_groups_5 = [1, 4, 10]
        time_syn_sample_groups_6 = [3, 12, 14]
        time_syn_sample_group_dict = {
            1: time_syn_sample_groups_1,
            2: time_syn_sample_groups_2,
            3: time_syn_sample_groups_3,
            4: time_syn_sample_groups_4,
            5: time_syn_sample_groups_5,
            6: time_syn_sample_groups_6
        }
        time_sample_group_label_dict = {
            1: "10 Species and Simulated",
            2: "10 Species and Cultured",
            3: "40 Species",
            4: "120 Species",
            5: "500 Species",
            6: "1000 Species"
        }
        # Sample sizes
        sample_size_dict = {
            1: 11.64,
            2: 5.2,
            3: 2.2,
            4: 2.0,
            5: 0.6966,
            6: 0.4431,
            7: 1.3,
            8: 0.0911,
            9: 5.3,
            10: 1.6,
            11: 0.7589,
            12: 3.5,
            13: 1.3,
            14: 2.7,
            15: 1.3,
            16: 5.2,
            17: 5.2,
            18: 0.2454,
            19: 0.144
        }
        sample_speciesnum_dict = {
            1: 500,
            2: 10,
            3: 1000,
            4: 500,
            5: 40,
            6: 40,
            7: 120,
            8: 10,
            9: 10,
            10: 500,
            11: 40,
            12: 1000,
            13: 120,
            14: 1000,
            15: 120,
            16: 10,
            17: 10,
            18: 10,
            19: 10
        }
        # The time needed to filter the report/results of Kraken2 is negligent compared to the time needed for the analysis to take place. Therefore, the time of analysis for each filtering threshold
        # of the Kraken2 report is based primarily (almost completely) to the execution time of ProteoSeeker for the kraken database that the filter is based on.
        methods_time_group = ["k8", "k16", "k72", "cnr", "mnr"]
        time_analysis(time_dir, time_stats_dir_path, time_syn_sample_group_dict, time_sample_group_label_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict)


if __name__ == "__main__":
    arg_benchmark_path = "12864_2022_8803_MOESM1_ESM.txt"
    arg_ps_results = ""
    arg_ps_output_dir = "analysis_plots_ps"
    arg_methods_group = "k8,k16,k72,k8_ng,k16_ng,k72_ng,cnr,mnr"
    arg_time_status = True
    if len(sys.argv) > 1:
        arg_input_command = "{}".format(sys.argv[0])
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input-db":
                arg_benchmark_path = sys.argv[i+1]
            if sys.argv[i] == "-p" or sys.argv[i] == "--proteoseeker-results":
                arg_ps_results = sys.argv[i+1]
            if sys.argv[i] == "-o" or sys.argv[i] == "--output":
                arg_ps_output_dir = sys.argv[i+1]
            if sys.argv[i] == "-m" or sys.argv[i] == "--methods":
                arg_methods_group = sys.argv[i+1]
            if sys.argv[i] == "-t" or sys.argv[i] == "--time":
                arg_time_status = sys.argv[i+1]
                arg_str = "-t"
                arg_time_status = check_i_value(arg_time_status, arg_str)
    benchstats(arg_benchmark_path, arg_ps_results, arg_ps_output_dir, arg_methods_group, arg_time_status)
