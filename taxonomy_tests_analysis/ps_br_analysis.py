import os
import sys
import copy
import math
import shutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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


def crfiles(ps_dir, ps_results, analysis_results_dict, time_dir):
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
                                kraken_dest_path = "{}/mdb_{}_sample_{}_{}".format(time_dir, method_db, sample_id, kn)
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
                    metabinner_dest_path_4 = "{}/mdb_{}_sample_{}_time_analysis.tsv".format(time_dir, method_db, sample_id)
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
                    comebin_dest_path_4 = "{}/mdb_{}_sample_{}_time_analysis.tsv".format(time_dir, method_db, sample_id)
                    if os.path.exists(comebin_source_path_1):
                        shutil.copyfile(comebin_source_path_1, comebin_dest_path_1)
                    if os.path.exists(comebin_source_path_2):
                        shutil.copyfile(comebin_source_path_2, comebin_dest_path_2)
                    if os.path.exists(comebin_source_path_3):
                        shutil.copyfile(comebin_source_path_3, comebin_dest_path_3)
                    if os.path.exists(comebin_source_path_4):
                        shutil.copyfile(comebin_source_path_4, comebin_dest_path_4)


def ps_kraken_analyze(ps_kraken_dir, filter_name):
    kraken_info_dict = {}
    kraken_files = os.listdir(ps_kraken_dir)
    # Kraken2
    for kf in kraken_files:
        if filter_name in kf:
            kf_nosuffix = kf[:-4]
            kf_splited = kf_nosuffix.split("_")
            sample_id = int(kf_splited[1])
            kraken_info_dict[sample_id] = {}
            relpath_kf = "{}/{}".format(ps_kraken_dir, kf)
            kraken_lines = read_file(relpath_kf)
            for line in kraken_lines:
                line_splited = line.split("\t")
                taxid = line_splited[0]
                if (taxid != "Total") and (taxid != "Toal"):
                    abundance = line_splited[1]
                    kraken_info_dict[sample_id][taxid] = abundance
    return kraken_info_dict


def ps_metabinner_analyze(ps_mtbr, filter_name):
    # MetaBinner
    mtbr_info_dict = {}
    mtbr_files = os.listdir(ps_mtbr)
    for mf in mtbr_files:
        if filter_name in mf:
            mf_nosuffix = mf[:-4]
            mf_splited = mf_nosuffix.split("_")
            sample_id = int(mf_splited[1])
            mtbr_info_dict[sample_id] = {}
            relpath_mf = "{}/{}".format(ps_mtbr, mf)
            mtbr_lines = read_file(relpath_mf)
            for line in mtbr_lines:
                line_splited = line.split("\t")
                taxid = line_splited[2]
                rank = line_splited[3]
                abundance = line_splited[5]
                if rank == "species":
                    mtbr_info_dict[sample_id][taxid] = abundance
    return mtbr_info_dict


def ps_comebin_analyze(ps_cmbn, filter_name):
    # COMEBin
    cmbn_info_dict = {}
    cmbn_files = os.listdir(ps_cmbn)
    for cf in cmbn_files:
        if filter_name in cf:
            cf_nosuffix = cf[:-4]
            cf_splited = cf_nosuffix.split("_")
            sample_id = int(cf_splited[1])
            cmbn_info_dict[sample_id] = {}
            relpath_cf = "{}/{}".format(ps_cmbn, cf)
            cmbn_lines = read_file(relpath_cf)
            for line in cmbn_lines:
                line_splited = line.split("\t")
                taxid = line_splited[2]
                rank = line_splited[3]
                abundance = line_splited[5]
                if rank == "species":
                    cmbn_info_dict[sample_id][taxid] = abundance
    return cmbn_info_dict


def basic_stats(br_info_dict, target_info_dict, target_stats_dict, sample_id):
    tp = 0
    tn = 0
    fn = 0
    fp = 0
    sensitivity = 0
    specificity = 0
    precision = 0
    f1_score = 0
    accuracy = 0
    jac_index = 0
    sorensen_dice = 0
    bray_curtis_dis = 0
    diff_sum_wr = 0
    mean_absolute_error = 0
    squared_diff_sum_wr = 0
    mean_squared_error = 0
    root_mean_squared_error = 0
    gold_si_info = br_info_dict[sample_id]
    target_si_info = target_info_dict[sample_id]
    # Kraken2
    # Create a group of all the items from gold and predicted items.
    unique_items = []
    for item in gold_si_info:
        if item not in unique_items:
            unique_items.append(item)
    for item in target_si_info:
        if item not in unique_items:
            unique_items.append(item)
    # Parse the list of unique items and check the existance of each item in each of the lists.
    # common_items: Common species
    # group_val_items: Gold species
    # group_val_unique_items: Unique gold species
    # group_pred_items: Predicted species
    # group_pred_unique_items: Unique predicted species
    common_items = []
    group_val_items = []
    group_val_unique_items = []
    group_pred_items = []
    group_pred_unique_items = []
    for item in unique_items:
        if (item in gold_si_info) and (item in target_si_info):
            tp += 1
            common_items.append(item)
            group_val_items.append(item)
            group_pred_items.append(item)
        elif (item not in gold_si_info) and (item not in target_si_info):
            tn += 1
        elif (item in gold_si_info) and (item not in target_si_info):
            fn += 1
            group_val_items.append(item)
            group_val_unique_items.append(item)
        else:
            fp += 1
            group_pred_items.append(item)
            group_pred_unique_items.append(item)
    # Sensitivity
    if (tp + fn) != 0:
        sensitivity =  tp / (tp + fn)
        sensitivity = round(sensitivity, 2)
    # Specificity
    if (tn + fp) != 0:
        specificity =  tn / (tn + fp)
        specificity = round(specificity, 2)
    # Precision
    if (tp + fp) != 0:
        precision = tp / (tp + fp)
        precision = round(precision, 2)
    # F1 Score
    if (precision + sensitivity) != 0:
        f1_score = 2 * ((precision * sensitivity) / (precision + sensitivity))
        f1_score = round(f1_score, 2)
    # Accuracy
    if (tp + tn + fp + fn) != 0:
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        accuracy = round(accuracy, 2)
    # Jaccard index:
    # The true positive hits are the number of the items of the intersection of the two gruops.
    # The number of unique items of the two groups is the number of the items of their union.
    unique_spec_num = len(unique_items)
    if unique_spec_num != 0:
        jac_index = tp / unique_spec_num
        jac_index = round(jac_index, 2)
    # Sorensen-Dice Index
    # The true positive hits are the number of the items of the intersection of the two gruops.
    val_species_num = len(gold_si_info)
    pred_species_num = len(target_si_info)
    if (val_species_num + pred_species_num) != 0:
        sorensen_dice = (2 * tp) / (val_species_num + pred_species_num)
        sorensen_dice = round(sorensen_dice, 2)
    # Bray-Curtis Dissimilarity
    # For each species from both groups, take the minimum percentage (0 if the species is not present in one of the groups.)
    # Sum these minimum percentagies, the sum the actual percentagies of the same species for each group seperatly.
    # Then compute the metric.
    min_sum = 0
    val_sum = 0
    pred_sum = 0
    for item in common_items:
        gold_abu_prec = float(br_info_dict[sample_id][item])
        pred_abu_perc = float(target_info_dict[sample_id][item])
        min_abu_perc = min(gold_abu_prec, pred_abu_perc)
        min_sum += min_abu_perc
    for item in group_val_items:
        gold_abu_prec = float(br_info_dict[sample_id][item])
        val_sum += gold_abu_prec
    for item in group_pred_items:
        pred_abu_perc = float(target_info_dict[sample_id][item])
        pred_sum += pred_abu_perc
    if (val_sum + pred_sum) != 0:
        bray_curtis_dis = (2 * min_sum) / (val_sum + pred_sum)
        bray_curtis_dis = round(bray_curtis_dis, 2)
    # Mean Absolute Error / MAE
    # For common species their percentagies is compared directly.
    # For species unique to each of the groups, the percentage of that species for the other group is considered 0.
    unique_spec_num = len(unique_items)
    unique_val_num = len(group_val_unique_items)
    unique_pred_num = len(group_pred_unique_items)
    if unique_spec_num != 0:
        diff_sum = 0
        for item in common_items:
            gold_abu_prec = float(br_info_dict[sample_id][item])
            pred_abu_perc = float(target_info_dict[sample_id][item])
            diff = gold_abu_prec - pred_abu_perc
            abs_diff = abs(diff)
            diff_sum += abs_diff
        for item in group_val_unique_items:
            gold_abu_prec = float(br_info_dict[sample_id][item])
            diff = gold_abu_prec
            abs_diff = abs(diff)
            diff_sum += abs_diff
        for item in group_pred_unique_items:
            pred_abu_perc = float(target_info_dict[sample_id][item])
            diff = pred_abu_perc
            abs_diff = abs(diff)
            diff_sum += abs_diff
        mean_absolute_error = diff_sum / unique_spec_num
        mean_absolute_error = round(mean_absolute_error, 2)
        diff_sum_wr = round(diff_sum, 2)
    # Mean Squared Error (MSE) and Root Mean Squared Error (RMSE)
    # Same methodlogy as for MAE.
    unique_spec_num = len(unique_items)
    unique_val_num = len(group_val_unique_items)
    unique_pred_num = len(group_pred_unique_items)
    if unique_spec_num != 0:
        squared_diff_sum = 0
        for item in common_items:
            gold_abu_prec = float(br_info_dict[sample_id][item])
            pred_abu_perc = float(target_info_dict[sample_id][item])
            diff = gold_abu_prec - pred_abu_perc
            abs_diff = abs(diff)
            abs_diff = abs_diff**2
            squared_diff_sum += abs_diff
        for item in group_val_unique_items:
            gold_abu_prec = float(br_info_dict[sample_id][item])
            diff = gold_abu_prec
            abs_diff = abs(diff)
            abs_diff = abs_diff**2
            squared_diff_sum += abs_diff
        for item in group_pred_unique_items:
            pred_abu_perc = float(target_info_dict[sample_id][item])
            diff = pred_abu_perc
            abs_diff = abs(diff)
            abs_diff = abs_diff**2
            squared_diff_sum += abs_diff
        mean_squared_error = squared_diff_sum / unique_spec_num
        mean_squared_error = round(mean_squared_error, 2)
        squared_diff_sum_wr = round(squared_diff_sum, 2)
        root_mean_squared_error = math.sqrt(mean_squared_error)
        root_mean_squared_error = round(root_mean_squared_error, 2)    
    # Store information in the dictionary
    target_stats_dict[sample_id] = {
        "Gold species number": val_species_num,
        "Predicted species number": pred_species_num,
        "Common species (intersection): ": tp,
        "Unique species from both groups (union)": unique_spec_num,
        "Unique to gold group": unique_val_num,
        "Unique to predicted group": unique_pred_num,
        "True Positive (TP)": tp,
        "True Negative (TN)": tn,
        "False Positive (FP)": fp,
        "False Negative (FN)": fn,
        "Sensitivity": sensitivity,
        "Specificty": specificity,
        "Precision": precision,
        "Accuracy": accuracy,
        "F1 score": f1_score,
        "Jaccard Index": jac_index,
        "Sorensen-Dice Index": sorensen_dice,
        "Bray-Curtis Dissimilarity": bray_curtis_dis,
        "Sum of absolute differences": diff_sum_wr,
        "MAE": mean_absolute_error,
        "Sum of squared absolute differencies": squared_diff_sum_wr,
        "MSE": mean_squared_error,
        "RMSE": root_mean_squared_error
    }
    # Write the information in a csv file.
    return target_stats_dict, common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items


def write_list(path_tow, list_tow):
    file_tow = open(path_tow, "w")
    for item in list_tow:
        file_tow.write("{}\n".format(item))
    file_tow.close()


def write_species(common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items, stats_dir_path, label, si):
    print("Writting information about the common and unique species for sample {} and method {}...".format(si, label))
    # File paths
    common_path = "{}/common_{}_sample_{}.txt".format(stats_dir_path, label, si)
    group_val_path = "{}/val_{}_sample_{}.txt".format(stats_dir_path, label, si)
    group_val_unique_path = "{}/val_unique_{}_sample_{}.txt".format(stats_dir_path, label, si)
    group_pred_path = "{}/pred_{}_sample_{}.txt".format(stats_dir_path, label, si)
    group_pred_unique_path = "{}/pred_unique_{}_sample_{}.txt".format(stats_dir_path, label, si)
    sample_synopsis_path = "{}/synopsis_{}_sample_{}.csv".format(stats_dir_path, label, si)
    write_list(common_path, common_items)
    write_list(group_val_path, group_val_items)
    write_list(group_val_unique_path, group_val_unique_items)
    write_list(group_pred_path, group_pred_items)
    write_list(group_pred_unique_path, group_pred_unique_items)
    # Dictionary with the species. It is needed to create the dataframe.
    species_dict = {
        "common_species": common_items,
        "validation_species": group_val_items,
        "validation_species_unique": group_val_unique_items,
        "predicted_species": group_pred_items,
        "predicted_species_unique": group_pred_unique_items
    }
    # Create the pandas dataframe for the synopsis data. Each empty cell is filled with the NaN value.
    col_labels = ["common_species", "validation_species", "validation_species_unique", "predicted_species", "predicted_species_unique"]
    syn_df = pd.DataFrame(columns=col_labels)
    syn_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in species_dict.items()]))
    # CSV file for one sample and one analysis with all the stats for the species
    syn_df.to_csv(sample_synopsis_path, index=False)


def comp_stats(br_info_dict, comp_dict, stats_dir_path, label):
    print("Comparing results with the validation data...\n")
    if os.path.exists(stats_dir_path):
        shutil.rmtree(stats_dir_path)
    os.mkdir(stats_dir_path)
    # Benchmark | Sample | Result
    #    Yes    |   Yes  |   TP
    #    No     |   No   |   TN
    #    Yes    |   No   |   FN
    #    No     |   Yes  |   FP
    # Sensitivity / Recall / True positive rate = TP / (TP + FN)
    # Specifcity = TN / (TN + FP) = 0
    # Precision = TP / (TP + FP)
    # F1 Score = 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
    # Accuracy = (TP + TN) / (TP + TN + FP + FN) = TP / (TP + FN)
    # Jaccard index = common_species / unique_species
    # Sorensen-Dice Index = 2 * common_species / (species_val + species_pred)
    #
    # Note: The results of MetaBinner and COMEBin may have more than one species for one bin.
    # All species for each bin are considered in the results. Also all species of the same bin
    # habe the same abundance.
    # 1. Parsing the species of the benchmark dataset and checking for their existance in the results of ProteoSeeker. Outputs TP and FN.
    # 2. Parsing the species of ProteoSeeker and checking for their existance in the results of the benchmark dataset. Outputs TP and FP.
    # No TN can come out of our analysis. No species are predicted by ProteoSeeker or the benchmark dataset to be absent from the sample
    # because none of the tools make negative predictions.
    sample_ids = []
    sample_ids_dupl = list(br_info_dict.keys())
    for item in sample_ids_dupl:
        if item in sample_ids:
            print("Error. Duplicate sample ID. Exiting.")
            exit()
    sample_ids = copy.deepcopy(sample_ids_dupl)
    sample_ids = sorted(sample_ids, key=int)
    # Compute sensitivity
    comb_stats_dict = {}
    for si in sample_ids:
        if si in comp_dict.keys():
            comb_stats_dict, common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items = basic_stats(br_info_dict, comp_dict, comb_stats_dict, si)
            write_species(common_items, group_val_items, group_val_unique_items, group_pred_items, group_pred_unique_items, stats_dir_path, label, si)
    print()
    return comb_stats_dict


def design_stacked_plots(comb_info_dict, metric_label_dict, legend_label_dict, plot_dir_parh):
    print("Plotting...\n")
    if os.path.exists(plot_dir_parh):
        shutil.rmtree(plot_dir_parh)
    os.mkdir(plot_dir_parh)
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
    metric_dict = {}
    sns.set_theme(style="whitegrid")
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
    pandas_dict = {}
    for key_metric in metric_dict.keys():
        # Initialize an empty Pandas dataframe.
        col_labels = ["sample", "method", "value"]
        metric_df = pd.DataFrame(columns=col_labels)
        row_index = 0
        for item in metric_dict[key_metric]:
            metric_df.loc[row_index] = item
            row_index += 1
        pandas_dict[key_metric] = metric_df
    for key_metric in pandas_dict.keys():
        df_metric = pandas_dict[key_metric]
        # File path
        csv_path = "{}/{}.csv".format(plot_dir_parh, key_metric)
        # Storing the pandas dataframe as a csv.
        df_metric.to_csv(csv_path, index=False)
        # Figure
        figure, cur_axis = plt.subplots(figsize=(20, 12))
        # Plotting the pandas dataframe: Grouped barplot
        sns.barplot(data=df_metric, x="sample", y="value", hue="method", palette="dark", alpha=.6, ax=cur_axis)
        # Place vertical lines to divie the sample barplots.
        for i in range(1, 19):
            cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
        # Legend title
        legend_title = "Taxonomy method\nand database"
        # Legend
        cur_axis.legend(title=legend_title, loc='center left', bbox_to_anchor=(1, 0.5))
        # Replace the legend's labels.
        for legend_text in cur_axis.get_legend().get_texts():
            cur_legend_label = legend_text.get_text()
            if cur_legend_label in legend_label_dict.keys():
                new_legnd_label = legend_label_dict[cur_legend_label]
                legend_text.set_text(new_legnd_label)
        # Layout
        plt.tight_layout(rect=[0, 0, 0.99, 0.99])
        # Axis labels
        metric_label = metric_label_dict[key_metric]
        x_axis_label = "Sample ID"
        y_axis_label = metric_label
        cur_axis.set_xlabel(x_axis_label)
        cur_axis.set_ylabel(y_axis_label)
        axis_title_label = "{} vs Sample ID".format(metric_label)
        cur_axis.set_title(axis_title_label)
        # File paths
        plot_png_path = "{}/{}.png".format(plot_dir_parh, key_metric)
        plot_jpg_path = "{}/{}.jpg".format(plot_dir_parh, key_metric)
        # Store the plots
        figure.savefig(plot_png_path)
        figure.savefig(plot_jpg_path)
        plt.close()
    return pandas_dict


def design_metric_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict, legend_label_dict, plot_dir_parh):
    # Loop:
    # key_metric = 1
    # cur_metric_group = ["True Positive (TP)", "False Positive (FP)", "False Negative (FN)"]
    # row_num = 3
    # Create figure
    #
    # item = "True Positive (TP)"
    # df_metric = Dataframe for "True Positive (TP)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # item = "False Positive (FP)"
    # df_metric = Dataframe for "False Positive (FP)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # item = "False Negative (FN)"
    # df_metric = Dataframe for "False Negative (FN)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # Save figure
    #
    # ...
    for key_group_metric in metric_group_dict.keys():
        cur_metric_group = metric_group_dict[key_group_metric]
        row_num = len(cur_metric_group)
        # Create the figure to store the subplots.
        figure, axis = plt.subplots(row_num, 1, figsize=(20, 5 * row_num)) 
        row_fig_index = 0
        for item in cur_metric_group:
            cur_axis = axis[row_fig_index]
            df_metric = pandas_dict[item]
            sns.barplot(ax=cur_axis, data=df_metric, x="sample", y="value", hue="method", palette="dark", alpha=.6)
            # Removes the right, top, left and bottom axis.
            sns.despine(ax=cur_axis, left=True)
            # Labels for x axis, y axis and title.
            metric_label = metric_label_dict[item]
            x_label_label = "Sample ID"
            y_axis_label = metric_label
            axis_title_label = "{} vs Sample ID".format(metric_label)
            cur_axis.set(xlabel=x_label_label, ylabel=y_axis_label)
            cur_axis.set_title(axis_title_label, pad=20, loc='center', fontsize=12)
            # Remove the legend, if for the last plot.
            if row_fig_index < row_num:
                cur_axis.legend().remove()
            # Replace the y axis ticks.
            cur_axis.set_yticks(list(cur_axis.get_yticks()))
            # Place vertical lines to divie the sample barplots.
            for i in range(1, 19):
                cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
            # Increasing the row index that shows which subplot will be added in the figure.
            row_fig_index += 1
        # Add one legend at the bottom. Make the legend as much horizontal as possible.
        # Legend
        legend_obj = cur_axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=11, title='Taxonomy Method and Database')
        legend_title = "Taxonomy method and database"
        legend_obj.set_title(legend_title, prop={'size': 10})
        # Replace the legend's labels.
        for legend_text in legend_obj.get_texts():
            cur_legend_label = legend_text.get_text()
            if cur_legend_label in legend_label_dict.keys():
                new_legnd_label = legend_label_dict[cur_legend_label]
                legend_text.set_text(new_legnd_label)
        # Padding from left, bottom, right, top
        plt.tight_layout()
        # Space between the plots of the figure.
        plt.subplots_adjust(hspace=0.4)
        # File paths
        group_metrics_path_png = "{}/group_metrics_{}.png".format(plot_dir_parh, key_group_metric)
        group_metrics_path_jpg = "{}/group_metrics_{}.jpg".format(plot_dir_parh, key_group_metric)
        # Saving the figure.
        figure.savefig(group_metrics_path_png)
        figure.savefig(group_metrics_path_jpg)
        # Closing the figure.
        plt.close()


def design_metric_sample_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict, legend_label_dict, sample_group_dict, sample_group_label_dict, plot_dir_parh):
    # Loop:
    # key_group = 1
    # sample_id_group = [2, 8, 9, 16, 17, 18, 19]
    #
    # key_metric = 1
    # cur_metric_group = ["True Positive (TP)", "False Positive (FP)", "False Negative (FN)"]
    # row_num = 3
    # Create figure
    #
    # item = "True Positive (TP)"
    # df_metric = Dataframe for "True Positive (TP)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # item = "False Positive (FP)"
    # df_metric = Dataframe for "False Positive (FP)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # item = "False Negative (FN)"
    # df_metric = Dataframe for "False Negative (FN)"
    # Convert dataframe to plot
    # Add plot to specific positio in the figure.
    #
    # Save figure
    # 
    # key_group = 2
    # sample_id_group = [5, 6, 11]
    # ...
    # Filter specific samples.
    for key_group in sample_group_dict.keys():
        sample_id_group = sample_group_dict[key_group]
        # Group label
        group_label = sample_group_label_dict[key_group]
        for key_group_metric in metric_group_dict.keys():
            cur_metric_group = metric_group_dict[key_group_metric]
            row_num = len(cur_metric_group)
            # Create the figure to store the subplots.
            figure, axis = plt.subplots(row_num, 1, figsize=(20, 5 * row_num)) 
            row_fig_index = 0
            # Checks whether the dataframes for the metrics will be empty for the current sample group.
            # If for one metric the filtered dataframe is empty then it will also be for the rest of the metrics.
            df_empty_status = False
            for item in cur_metric_group:
                df_metric = pandas_dict[item]
                df_metric_filtered = df_metric[df_metric['sample'].isin(sample_id_group)]
                if df_metric_filtered.empty:
                    df_empty_status = True
                    break
                # Placing the plot in the figure.
                cur_axis = axis[row_fig_index]
                sns.barplot(ax=cur_axis, data=df_metric_filtered, x="sample", y="value", hue="method", palette="dark", alpha=.6)
                # Removes the right, top, left and bottom axis.
                sns.despine(ax=cur_axis, left=True)
                # Labels for x axis, y axis and title.
                metric_label = metric_label_dict[item]
                x_label_label = "Sample ID"
                y_axis_label = metric_label
                axis_title_label = "{} vs Sample ID for {}".format(metric_label, group_label)
                cur_axis.set(xlabel=x_label_label, ylabel=y_axis_label)
                cur_axis.set_title(axis_title_label, pad=20, loc='center', fontsize=12)
                # Remove the legend, if for the last plot.
                if row_fig_index < row_num:
                    cur_axis.legend().remove()
                # Replace the y axis ticks.
                cur_axis.set_yticks(list(cur_axis.get_yticks()))
                # Place vertical lines to divie the sample barplots.
                for i in range(1, len(sample_id_group)):
                    cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
                # Increasing the row index that shows which subplot will be added in the figure.
                row_fig_index += 1
            if df_empty_status:
                break
            # Add one legend at the bottom. Make the legend as much horizontal as possible.
            # Legend
            legend_obj = cur_axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=11, title='Taxonomy Method and Database')
            legend_title = "Taxonomy method and database"
            legend_obj.set_title(legend_title, prop={'size': 10})
            # Replace the legend's labels.
            for legend_text in legend_obj.get_texts():
                cur_legend_label = legend_text.get_text()
                if cur_legend_label in legend_label_dict.keys():
                    new_legnd_label = legend_label_dict[cur_legend_label]
                    legend_text.set_text(new_legnd_label)
            # Padding from left, bottom, right, top
            plt.tight_layout()
            # Space between the plots of the figure.
            plt.subplots_adjust(hspace=0.4)
            # File paths
            group_metrics_path_png = "{}/group_{}_metrics_{}.png".format(plot_dir_parh, group_label, key_group_metric)
            group_metrics_path_jpg = "{}/group_{}_metrics_{}.jpg".format(plot_dir_parh, group_label, key_group_metric)
            # Saving the figure.
            figure.savefig(group_metrics_path_png)
            figure.savefig(group_metrics_path_jpg)
            # Closing the figure.
            plt.close()


def check_add(add_item, target_num):
    if add_item != "None":
        add_item = float(add_item)
        target_num += add_item
    return target_num


def collect_times(time_dir):
    pre_time_dict = {}
    time_file_names = os.listdir(time_dir)
    for tfn in time_file_names:
        tfn_spliited = tfn.split("_")
        method_dp = tfn_spliited[1]
        sample_id = tfn_spliited[3]
        # Dictionary information
        if sample_id not in pre_time_dict.keys():
            pre_time_dict[sample_id] = {}
        pre_time_dict[sample_id][method_dp] = {}
        # File path
        tfn_path = "{}/{}".format(time_dir, tfn)
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
    orderded_methods = ["k8", "k16", "k72", "mnr", "cnr"]
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
            k8_cd_hit_time = time_dict[key_sample]["k8"]["cd_hit_time"]
            if (k8_fastqc_init_time == "None") or (k8_preproc_time == "None") or (k8_fastqc_final_time == "None") or (k8_assembly_time == "None") or (k8_gene_pred_time == "None") or (k8_cd_hit_time == "None"):
                sum_common_times = None
            else:
                k8_fastqc_init_time = float(k8_fastqc_init_time)
                k8_preproc_time = float(k8_preproc_time)
                k8_fastqc_final_time = float(k8_fastqc_final_time)
                k8_assembly_time = float(k8_assembly_time)
                k8_gene_pred_time = float(k8_gene_pred_time)
                k8_cd_hit_time = float(k8_cd_hit_time)
                sum_common_times = k8_fastqc_init_time + k8_preproc_time + k8_fastqc_final_time + k8_assembly_time + k8_gene_pred_time
        # Perform the additions
        if sum_common_times is None:
            print("Error in computing the additional times. Continuing to the next sample.\n")
            for key_method in time_dict[key_sample].keys():
                time_tool_dict[key_sample][key_method] = None
            continue
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
    for key_sample in time_dict.keys():
        if "k8" in time_dict[key_sample].keys():
            k8_fastqc_init_time = time_dict[key_sample]["k8"]["fastqc_initial_time"]
            k8_preproc_time = time_dict[key_sample]["k8"]["preprocessing_time"]
            k8_fastqc_final_time = time_dict[key_sample]["k8"]["fastqc_final_time"]
            k8_assembly_time = time_dict[key_sample]["k8"]["assembly_time"]
            k8_gene_pred_time = time_dict[key_sample]["k8"]["gene_prediction_time"]
            #k8_cd_hit_time = time_dict[key_sample]["k8"]["cd_hit_time"]
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


def create_df_stacked(time_dict):
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
    all_col_labels = ["sample", "method", "fastqc", "preprocessing", "assembly", "gene_prediction", "cd_hit", "binning", "taxonomy", "bowtie", "results"]
    df_full_time_all = pd.DataFrame(columns=all_col_labels)
    df_time_dict = {}
    row_index_all = 0
    for key_sample in time_dict.keys():
        # Initialize an empty Pandas dataframe.
        col_labels = ["method", "fastqc", "preprocessing", "assembly", "gene_prediction", "cd_hit", "binning", "taxonomy", "bowtie", "results"]
        time_metric_df = pd.DataFrame(columns=col_labels)
        row_index = 0
        for key_method in time_dict[key_sample].keys():
            df_total_time_str = time_dict[key_sample][key_method]["tool_time"]
            #df_sra_time_str = time_dict[key_sample][key_method]["sra_time"]
            #df_dbs_time_str = time_dict[key_sample][key_method]["sra_time"]
            df_fastqc_init_time_str = time_dict[key_sample][key_method]["fastqc_initial_time"]
            df_preproc_time_str = time_dict[key_sample][key_method]["preprocessing_time"]
            df_fastqc_final_time_str = time_dict[key_sample][key_method]["fastqc_final_time"]
            df_assembly_time_str = time_dict[key_sample][key_method]["assembly_time"]
            df_gene_pred_time_str = time_dict[key_sample][key_method]["gene_prediction_time"]
            #df_gene_anno_time_str = time_dict[key_sample][key_method]["gene_annotation_time"]
            df_cd_hit_time_str = time_dict[key_sample][key_method]["cd_hit_time"]
            df_kraken_tax_time_str = time_dict[key_sample][key_method]["kraken_time"]
            # df_kraken_tax_spec_time_str = time_dict[key_sample][key_method]["kraken_specific_time"]
            df_cm_binning_time_str = time_dict[key_sample][key_method]["binning_time"]
            df_bowtie_time_str = time_dict[key_sample][key_method]["bowtie_time"]
            # df_hmmer_spec_time_str = time_dict[key_sample][key_method]["hmmer_spec_time"]
            df_hmmer_tax_spec_time_str = time_dict[key_sample][key_method]["hmmer_spec_taxonomy_time"]
            # df_blastp_1_time_str = time_dict[key_sample][key_method]["blastp_fpd_no_doms_time"]
            # df_blastp_2_time_str = time_dict[key_sample][key_method]["blastp_fpd_swiss_doms_time"]
            df_blastp_3_time_str = time_dict[key_sample][key_method]["blastp_fpd_swiss_taxonomy_time"]
            df_cm_analysis_time_str = time_dict[key_sample][key_method]["bin_analysis_cm_time"]
            df_cm_tax_time_str = time_dict[key_sample][key_method]["bin_taxonomy_cm_time"]
            df_kraken_binning_time_str = time_dict[key_sample][key_method]["kraken_binning_time"]
            # df_hmmer_broad_time_str = time_dict[key_sample][key_method]["hmmer_broad"]
            # df_topology_time_str = time_dict[key_sample][key_method]["topology_time"]
            # df_motifs_time_str = time_dict[key_sample][key_method]["motifs_time"]
            # df_fam_pred_time_str = time_dict[key_sample][key_method]["family_prediction_time"]
            df_info_time_str = time_dict[key_sample][key_method]["info_collection_time"]
            df_results_time_str = time_dict[key_sample][key_method]["results_time"]
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
            df_total_time_m, df_total_time_h = convert_to_minutes(df_total_time)
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
        df_time_dict[key_sample] = time_metric_df
    return df_time_dict, df_full_time_all


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
            total_time = total_time / 60
            total_time_m = round(total_time, 2)
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
            total_time = total_time / 60
            total_time_m = round(total_time, 2)
            key_sample_int = int(key_sample)
            temp_list = [key_sample_int, key_method, total_time]
            total_time_all_df.loc[row_index] = temp_list
            row_index += 1
    return df_total_time_dict, total_time_all_df


def count_megahit_time(df_full_time_dict):
    # Counter to sum the time up untill the stage of megahit.
    megahit_time_dict = {}
    for key_sample in df_full_time_dict.keys():
        full_time_metric_df = df_full_time_dict[key_sample]
        sums_up_to_assembly_list = full_time_metric_df.apply(lambda row: row[:"assembly"].sum(), axis=1)
        sums_up_to_assembly_list = sums_up_to_assembly_list.tolist()
        # Check if all sums are the same. If not exit with error.
        item_0 = sums_up_to_assembly_list[0]
        for item in sums_up_to_assembly_list:
            if item_0 != item:
                print("Error. Different sum of time up untill the megahit stage. Exiting.\n")
                exit()
        item_0 = float(item_0)
        sum_up_to_assembly = round(item_0, 2)
        megahit_time_dict[key_sample] = sum_up_to_assembly
    return megahit_time_dict


def plot_full_sample_time(df_full_time_dict, megahit_time_dict, time_dir, stats_dir_path):
    print("Full time plots for each sample...\n")
    # Grouped barplots
    for key_sample in df_full_time_dict.keys():
        # Get the dataframe from the dictionary.
        full_time_metric_df = df_full_time_dict[key_sample]
        # Get the time for megahit.
        cur_megahit_time = megahit_time_dict[key_sample]
        # Create the filename.
        full_time_df_png_path = "{}/sample_{}_full_time.png".format(time_dir, key_sample)
        full_time_df_jpg_path = "{}/sample_{}_full_time.jpg".format(time_dir, key_sample)
        stacked_csv_path = "{}/sample_{}_full_time.csv".format(stats_dir_path, key_sample)
        # Save to CSV
        full_time_metric_df.to_csv(stacked_csv_path, index=True)
        # Create the plot
        ax = full_time_metric_df.plot(kind='bar', stacked=True)
        # Place a horizontal line at the sum of time up to megahit.
        ax.axhline(y=cur_megahit_time, color='r', linestyle='--', linewidth=1)
        # Labels
        plt.xlabel('Taxonomy mode route and database', fontsize=10)
        plt.ylabel('Stage times (m)', fontsize=10)
        plt.title('Execution time vs Taxonomy mode route and database', fontsize=10)
        # Place the legend outside the plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Layout
        plt.tight_layout(rect=[0, 0, 0.8, 1])
        # Save the plot
        plt.savefig(full_time_df_png_path)
        plt.savefig(full_time_df_jpg_path)
        # Close the plot
        plt.close()


def plot_total_sample_time(df_total_time_dict, megahit_time_dict, time_dir, stats_dir_path):
    print("Total time plot for each sample...\n")
    # Create a barplot for each sample.
    for key_sample in df_total_time_dict.keys():
        # Get the dataframe from the dictionary.
        total_time_sample_df = df_total_time_dict[key_sample]
        # Get the time for megahit.
        cur_megahit_time = megahit_time_dict[key_sample]
        # Create the filename.
        total_time_sample_df_png_path = "{}/sample_{}_total_time.png".format(time_dir, key_sample)
        total_time_sample_df_jpg_path = "{}/sample_{}_total_time.jpg".format(time_dir, key_sample)
        stacked_csv_path = "{}/sample_{}_total_time.csv".format(stats_dir_path, key_sample)
        # Save to CSV
        total_time_sample_df.to_csv(stacked_csv_path, index=True)
        # Create a new figure
        plt.figure()
        # Create the plot
        sns.barplot(x='method', y='total_time', data=total_time_sample_df)
        # Place a horizontal line at the sum of time up to megahit.
        plt.axhline(y=cur_megahit_time, color='r', linestyle='--', linewidth=1)
        # Labels
        plt.xlabel('Analysis type', fontsize=10)
        plt.ylabel('Total time (m)', fontsize=10)
        plt.title('Total execution time vs Analysis type', fontsize=10)
        # Save the plot
        plt.savefig(total_time_sample_df_png_path)
        plt.savefig(total_time_sample_df_jpg_path)
        # Close the plot
        plt.close()


def plot_total_all_time(total_time_all_df, time_dir, stats_dir_path):
    print("Total time plot for all samples...\n")
    # Create a grouped barplot for the total times
    total_time_all_csv_path = "{}/total_time_all.csv".format(stats_dir_path)
    plot_png_path = "{}/total_time_all.png".format(time_dir)
    plot_jpg_path = "{}/total_time_all.jpg".format(time_dir)
    # Storing the pandas dataframe as a csv.
    total_time_all_df.to_csv(total_time_all_csv_path, index=False)
    # Plotting the pandas dataframe: Barplot
    total_time_all_barplot = sns.catplot(data=total_time_all_df, kind="bar", x="sample", y="total_time", hue="method", palette="dark", alpha=.6, height=6)
    total_time_all_barplot.despine(left=True)
    y_axis_label = ""
    total_time_all_barplot.set_axis_labels("", y_axis_label)
    total_time_all_barplot.legend.set_title("")
    total_time_all_barplot.savefig(plot_png_path)
    total_time_all_barplot.savefig(plot_jpg_path)
    plt.close()


def plot_full_group_sample_time(df_full_time_all, time_dir, stats_dir_path, time_sample_group_dict, time_sample_group_label_dict, methods_group):
    print("Full time plot for all samples grouped...\n")
    # The number of the methods to plot.
    methods_num = len(methods_group)
    # Save the dataframe to a CSV file.
    full_time_sample_grouped_time_csv_path = "{}/full_time_sample_grouped_time.csv".format(stats_dir_path)
    df_full_time_all.to_csv(full_time_sample_grouped_time_csv_path, index=True)
    # Create the figure to store the subplots.
    row_num = 5
    figure, axis = plt.subplots(row_num, 1, figsize=(20, 5 * row_num)) 
    row_fig_index = 0
    # Filter specific samples.
    for key_group in time_sample_group_dict.keys():
        sample_id_group = time_sample_group_dict[key_group]
        # Group label
        group_label = time_sample_group_label_dict[key_group]
        # Filtering the dataframe.
        df_time_filtered = df_full_time_all[df_full_time_all['sample'].isin(sample_id_group)]
        if df_time_filtered.empty:
            break
        # Create the plot
        cur_axis = axis[row_fig_index]
        df_time_filtered.set_index(['sample', 'method'], inplace=True)
        df_time_filtered.plot(kind='bar', stacked=True, ax=cur_axis, alpha=.6)
        cur_axis.set_xticklabels([f"{sample}-{method}" for sample, method in df_time_filtered.index], rotation=45, ha="right")
        # Labels for x axis, y axis and title.
        x_label_label = "Sample ID"
        y_axis_label = "Time (m)"
        axis_title_label = "{} vs Sample ID for {}".format(y_axis_label, group_label)
        cur_axis.set(xlabel=x_label_label, ylabel=y_axis_label)
        cur_axis.set_title(axis_title_label, pad=20, loc='center', fontsize=12)
        # Remove the legend, if for the last plot.
        if row_fig_index < row_num:
            cur_axis.legend().remove()
        # Replace the y axis ticks.
        cur_axis.set_yticks(list(cur_axis.get_yticks()))
        # Place vertical lines to divie the sample barplots.
        max_ver_line_col = len(sample_id_group) * methods_num
        for i in range(methods_num, max_ver_line_col, methods_num):
            cur_axis.axvline(x=i - 0.5, color='gray', linestyle='--', linewidth=0.8)
        # Increasing the row index that shows which subplot will be added in the figure.
        row_fig_index += 1
    # Legend
    legend_obj = cur_axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=11, title='Taxonomy Method and Database')
    legend_title = "Taxonomy method and database"
    legend_obj.set_title(legend_title, prop={'size': 10})
    # Padding from left, bottom, right, top.
    plt.tight_layout()
    # Space between the plots of the figure.
    plt.subplots_adjust(hspace=0.4)
    # File paths.
    full_time_grouped_path_png = "{}/full_time_sample_grouped_time.png".format(time_dir)
    full_time_grouped_path_jpg = "{}/full_time_sample_grouped_time.jpg".format(time_dir)
    # Saving the figure.
    figure.savefig(full_time_grouped_path_png)
    figure.savefig(full_time_grouped_path_jpg)
    # Closing the figure.
    plt.close()


def plot_full_group_method_time(total_time_all_df, time_dir, time_sample_group_dict, time_sample_group_label_dict, megahit_time_dict):
    print("Full time plot for all methods grouped...\n")
    # Create the figure to store the subplots.
    col_num = 5
    figure, axis = plt.subplots(1, col_num, figsize=(40, 20))
    row_fig_index = 0
    # Filter specific samples.
    for key_group in time_sample_group_dict.keys():
        sample_id_group = time_sample_group_dict[key_group]
        # Group label
        group_label = time_sample_group_label_dict[key_group]
        # Filtering the dataframe.
        total_time_all_filtered_df = total_time_all_df[total_time_all_df['sample'].isin(sample_id_group)]
        cur_axis = axis[row_fig_index]
        if total_time_all_filtered_df.empty:
            break
        sns.barplot(ax=cur_axis, data=total_time_all_filtered_df, x="sample", y="total_time", hue="method", palette="dark", alpha=.6)
        # Removes the right and top axis.
        sns.despine(ax=cur_axis)
        # Replace the y axis ticks.
        cur_axis.set_yticks(list(cur_axis.get_yticks()))
        # Remove the legend, if for the last plot.
        if row_fig_index < col_num:
            cur_axis.legend().remove()
        # Labels
        cur_axis.set_xlabel(group_label)
        cur_axis.set_ylabel("Total Time")
        # Place vertical lines to divie the sample barplots.
        row_fig_index += 1
    # Legend
    legend_obj = cur_axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=11, title='Taxonomy Method and Database')
    legend_title = "Taxonomy method and database"
    legend_obj.set_title(legend_title, prop={'size': 10})
    # Padding from left, bottom, right, top
    plt.tight_layout()
    full_time_grouped_path_png = "{}/full_time_method_grouped_time.png".format(time_dir)
    full_time_grouped_path_jpg = "{}/full_time_method_grouped_time.jpg".format(time_dir)
    # Saving the figure.
    figure.savefig(full_time_grouped_path_png)
    figure.savefig(full_time_grouped_path_jpg)
    # Closing the figure.
    plt.close()


def plot_total_size(df_total_time_dict, sample_size_dict, methods_group, sample_speciesnum_dict, time_dir):
    # Create a dictionary where the first key is the method and the second key the sample ID and its value the total time.
    method_sample_time_dict = {}
    for key_sample in df_total_time_dict.keys():
        temp_df = df_total_time_dict[key_sample]
        for mg in methods_group:
            if mg not in method_sample_time_dict.keys():
                method_sample_time_dict[mg] = {}
            mg_total_time = temp_df.loc[temp_df['method'] == mg, 'total_time'].iloc[0]
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
    # Create the scatter plots for size vs time.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        # Create figure
        plt.figure(figsize=(10, 10))
        cur_axis = plt.gca()
        # Scatter plot
        sns.scatterplot(data=plot_df, x="size", y="total_time", ax=cur_axis)
        # Regression plot
        # A line with an area of 95% confidence level can be added by removing the ci=None argument.
        sns.regplot(data=plot_df, x='size', y='total_time', scatter=False, color='blue', line_kws={"linewidth":1}, ci=None, ax=cur_axis)
        # Set the minimum value for the x and y axis.
        cur_axis.set_xlim(left=0)
        cur_axis.set_ylim(bottom=0)
        # Labels
        plt.xlabel('Sample size (GBs)')
        plt.ylabel('Total time (m)')
        plt.title(("For {}: Sample size vs Total time".format(key_method)))
        plt.tight_layout()
        size_time_png_path = "{}/{}_size_time.png".format(time_dir, key_method)
        size_time_jpg_path = "{}/{}_size_time.png".format(time_dir, key_method)
        plt.savefig(size_time_png_path)
        plt.savefig(size_time_jpg_path)
        plt.close()
    # Create the scatter plots for species number vs time.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        # Create figure
        plt.figure(figsize=(10, 10))
        cur_axis = plt.gca()
        # Returns a list of colors or continuous colormap defining a palette.
        uni_colors = sns.color_palette('hsv', len(plot_df['species_number'].unique()))
        sns.scatterplot(data=plot_df, x="species_number", y="total_time", hue="species_number", palette=uni_colors, ax=cur_axis)
        # Regression plot
        # A line with an area of 95% confidence level can be added by removing the ci=None argument.
        sns.regplot(data=plot_df, x='species_number', y='total_time', scatter=False, color='blue', line_kws={"linewidth":1}, ci=None, ax=cur_axis)
        # Set the minimum value for the x and y axis.
        cur_axis.set_xlim(left=0)
        cur_axis.set_ylim(bottom=0)
        # Legend
        plt.legend(title='Species Number', bbox_to_anchor=(1, 0.5), loc='center left')
        # Labels
        plt.xlabel('Species number')
        plt.ylabel('Total time (m)')
        plt.title(("For {}: Species number vs Total time".format(key_method)))
        plt.tight_layout()
        size_time_png_path = "{}/{}_species_time.png".format(time_dir, key_method)
        size_time_jpg_path = "{}/{}_species_time.png".format(time_dir, key_method)
        plt.savefig(size_time_png_path)
        plt.savefig(size_time_jpg_path)
        plt.close()
    # Create the scatter plots for species number vs mean time for each group of species number.
    for key_method in method_df_dict.keys():
        plot_df = method_df_dict[key_method]
        mean_plot_df = plot_df.groupby('species_number')['total_time'].mean().reset_index()
        # Create figure
        plt.figure(figsize=(10, 10))
        cur_axis = plt.gca()
        # Scatter plot
        sns.scatterplot(data=mean_plot_df, x="species_number", y="total_time", hue="species_number", ax=cur_axis)
        # Regression plot
        # A line with an area of 95% confidence level can be added by removing the ci=None argument.
        sns.regplot(data=mean_plot_df, x='species_number', y='total_time', scatter=False, color='blue', line_kws={"linewidth":1}, ci=None, ax=cur_axis)
        # Set the minimum value for the x and y axis.
        cur_axis.set_xlim(left=0)
        cur_axis.set_ylim(bottom=0)
        # Legend
        plt.legend(title='Species Number', bbox_to_anchor=(1, 0.5), loc='center left')
        # Labels
        plt.xlabel('Species number')
        plt.ylabel('Mean total time (m)')
        plt.title(("For {}: Species number vs Mean total time".format(key_method)))
        plt.tight_layout()
        size_time_png_path = "{}/{}_species_mean_time.png".format(time_dir, key_method)
        size_time_jpg_path = "{}/{}_species_mean_time.png".format(time_dir, key_method)
        plt.savefig(size_time_png_path)
        plt.savefig(size_time_jpg_path)
        plt.close()


def time_analysis(time_dir, stats_dir_path, time_sample_group_dict, time_sample_group_label_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict):
    print("Perforing analysis of the execution time...\n")
    # Collect the time periods
    time_dict = collect_times(time_dir)

    # Seperate the time stages
    time_full_dict, time_tool_dict = time_stages(time_dict)

    # Fill the full dict with the missing values.
    time_full_dict = fill_times(time_dict, time_full_dict)

    # Create the dataframe for the complete time slots.
    df_full_time_dict, df_full_time_all = create_df_stacked(time_full_dict)

    # Create the dataframe for the relative time slots.
    plain_time_metric_df = create_df_stacked(time_dict)
    
    # Create the dataframe for the total time.
    df_total_time_dict, total_time_all_df = create_df_total_time(time_tool_dict)

    # Count the time up to the stage of megahit for each sample.
    megahit_time_dict = count_megahit_time(df_full_time_dict)

    # Plot the dataframe
    plot_full_sample_time(df_full_time_dict, megahit_time_dict, time_dir, stats_dir_path)
    plot_total_sample_time(df_total_time_dict, megahit_time_dict, time_dir, stats_dir_path)
    plot_total_all_time(total_time_all_df, time_dir, stats_dir_path)
    plot_full_group_sample_time(df_full_time_all, time_dir, stats_dir_path, time_sample_group_dict, time_sample_group_label_dict, methods_time_group)
    plot_full_group_method_time(total_time_all_df, time_dir, time_sample_group_dict, time_sample_group_label_dict, megahit_time_dict)

    # Plot total time vs size
    plot_total_size(df_total_time_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict, time_dir)


def benchstats(benchmark_path="12864_2022_8803_MOESM1_ESM.txt", ps_results="", ps_output_dir="analysis_plots_ps", methods_group="k8,k16,k72,k8_ng,k16_ng,k72_ng,cnr,mnr", time_status=True):
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
    stats_dir_path = "{}/stats".format(ps_output_dir)
    # Directories for the results
    ps_kraken_8 = "{}/kraken_8".format(ps_dir)
    ps_kraken_16 = "{}/kraken_16".format(ps_dir)
    ps_kraken_72 = "{}/kraken_72".format(ps_dir)
    # Add the filtered results for Kraken.
    ps_mtbr_nr = "{}/metabinner_nr".format(ps_dir)
    ps_cmbn_nr = "{}/comebin_nr".format(ps_dir)

    analysis_results_dict = {
        "k8": "{}/kraken_8".format(ps_dir),
        "k16": "{}/kraken_16".format(ps_dir),
        "k72": "{}/kraken_72".format(ps_dir),
        "cnr": "{}/comebin_nr".format(ps_dir),
        "mnr": "{}/metabinner_nr".format(ps_dir),
    }
    
    # Analyze the bencharking file.
    if benchmark_path == "":
        print("Error. No input file was found. Exiting.\n")
        exit()

    # Determine the methods for which to get statistics and plots
    # methods_group = "k8,k16,k72,k8_ng,k16_ng,k72_ng,cnr,mnr"
    # methods_group = "k8,k16,k72,k8_0c1,k16_0c1,k72_0c1,k8_ng,k16_ng,k72_ng,cnr,mnr"
    # methods_group = "k8,k16,k72,k8_1c0,k16_1c0,k72_1c0,k8_ng,k16_ng,k72_ng,cnr,mnr"
    # methods_group = "k8,k16,k72,k8_10,k16_10,k72_10,k8_ng,k16_ng,k72_ng,cnr,mnr"
    # methods_group = "k8,k16,k72,k8_100,k16_100,k72_100,k8_ng,k16_ng,k72_ng,cnr,mnr"
    # methods_group = "k8,k16,k72,k8_0c1,k16_0c1,k72_0c1,k8_1c0,k16_1c0,k72_1c0,k8_10,k16_10,k72_10,k8_100,k16_100,k72_100,k8_ng,k16_ng,k72_ng"
    methods_group = methods_group.split(",")

    # Analyzing the benchmark dataset.
    br_info_dict = br_analysis(benchmark_path, benchmark_mod_path, benchmark_mod_spabu_path)
    
    if not br_info_dict:
        print("Error. No information from file 1. Exiting.")
        exit()

    # Copy the files with the results in the proper directories and rename them.
    crfiles(ps_dir, ps_results, analysis_results_dict, time_dir)

    # Collect the information from the results of ProteoSeeker
    # Kraken2
    filter_name_k_base = "_kraken_species.tsv"
    if "k8" in methods_group:
        kraken_8_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_base)
    if "k16" in methods_group:
        kraken_16_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_base)
    if "k72" in methods_group:
        kraken_72_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_base)
    # Kraken2: non-gut
    filter_name_k_f_base = "thr_-1_"
    if "k8_ng" in methods_group:
        kraken_8_ng_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base)
    if "k16_ng" in methods_group:
        kraken_16_ng_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base)
    if "k72_ng" in methods_group:
        kraken_72_ng_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base)
    # Kraken2: 0.1%
    filter_name_k_f_base = "thr_0.1_"
    if "k8_0c1" in methods_group:
        kraken_8_0c1_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base)
    if "k16_0c1" in methods_group:
        kraken_16_0c1_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base)
    if "k72_0c1" in methods_group:
        kraken_72_0c1_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base)
    # Kraken2: 1%
    filter_name_k_f_base = "thr_1.0_"
    if "k8_1c0" in methods_group:
        kraken_8_1c0_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base)
    if "k16_1c0" in methods_group:
        kraken_16_1c0_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base)
    if "k72_1c0" in methods_group:
        kraken_72_1c0_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base)
    # Kraken2: 10
    filter_name_k_f_base = "thr_10_"
    if "k8_10" in methods_group:
        kraken_8_10_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base)
    if "k16_10" in methods_group:
        kraken_16_10_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base)
    if "k72_10" in methods_group:
        kraken_72_10_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base)
    # Kraken2: 100
    filter_name_k_f_base = "thr_100_"
    if "k8_100" in methods_group:
        kraken_8_100_info_dict = ps_kraken_analyze(ps_kraken_8, filter_name_k_f_base)
    if "k16_100" in methods_group:
        kraken_16_100_info_dict = ps_kraken_analyze(ps_kraken_16, filter_name_k_f_base)
    if "k72_100" in methods_group:
        kraken_72_100_info_dict = ps_kraken_analyze(ps_kraken_72, filter_name_k_f_base)
    # COMEBin nr
    filter_name_c = "_b_summary_info_comebin.tsv"
    if "cnr" in methods_group:
        cmbn_nr_info_dict = ps_comebin_analyze(ps_cmbn_nr, filter_name_c)
    # MetaBinner nr
    filter_name_m =  "_b_summary_info_metabinner.tsv"
    if "mnr" in methods_group:
        mtbr_nr_info_dict = ps_metabinner_analyze(ps_mtbr_nr, filter_name_m)

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
    # Compare the results of ProteoSeeker with the benchmarking dataset.
    # Kraken2
    if "k8" in methods_group:
        label = "kraken_8"
        kraken_8_stats_dict = comp_stats(br_info_dict, kraken_8_info_dict, stats_dir_path, label)
    if "k16" in methods_group:
        label = "kraken_16"
        kraken_16_stats_dict = comp_stats(br_info_dict, kraken_16_info_dict, stats_dir_path, label)
    if "k72" in methods_group:
        label = "kraken_72"
        kraken_72_stats_dict = comp_stats(br_info_dict, kraken_72_info_dict, stats_dir_path, label)
    # Kraken2 0.1%
    if "k8_0c1" in methods_group:
        label = "kraken_8_0.1%"
        kraken_8_0c1_stats_dict = comp_stats(br_info_dict, kraken_8_0c1_info_dict, stats_dir_path, label)
    if "k16_0c1" in methods_group:
        label = "kraken_16_0.1%"
        kraken_16_0c1_stats_dict = comp_stats(br_info_dict, kraken_16_0c1_info_dict, stats_dir_path, label)
    if "k72_0c1" in methods_group:
        label = "kraken_72_0.1%"
        kraken_72_0c1_stats_dict = comp_stats(br_info_dict, kraken_72_0c1_info_dict, stats_dir_path, label)
    # Kraken2 1.0%
    if "k8_1c0" in methods_group:
        label = "kraken_8_1.0%"
        kraken_8_1c0_stats_dict = comp_stats(br_info_dict, kraken_8_1c0_info_dict, stats_dir_path, label)
    if "k16_1c0" in methods_group:
        label = "kraken_16_1.0%"
        kraken_16_1c0_stats_dict = comp_stats(br_info_dict, kraken_16_1c0_info_dict, stats_dir_path, label)
    if "k72_1c0" in methods_group:
        label = "kraken_72_1.0%"
        kraken_72_1c0_stats_dict = comp_stats(br_info_dict, kraken_72_1c0_info_dict, stats_dir_path, label)
    # Kraken2 10
    if "k8_10" in methods_group:
        label = "kraken_8_10"
        kraken_8_10_stats_dict = comp_stats(br_info_dict, kraken_8_10_info_dict, stats_dir_path, label)
    if "k16_10" in methods_group:
        label = "kraken_16_10"
        kraken_16_10_stats_dict = comp_stats(br_info_dict, kraken_16_10_info_dict, stats_dir_path, label)
    if "k72_10" in methods_group:
        label = "kraken_72_10"
        kraken_72_10_stats_dict = comp_stats(br_info_dict, kraken_72_10_info_dict, stats_dir_path, label)
    # Kraken2 100
    if "k8_100" in methods_group:
        label = "kraken_8_100"
        kraken_8_100_stats_dict = comp_stats(br_info_dict, kraken_8_100_info_dict, stats_dir_path, label)
    if "k16_100" in methods_group:
        label = "kraken_16_100"
        kraken_16_100_stats_dict = comp_stats(br_info_dict, kraken_16_100_info_dict, stats_dir_path, label)
    if "k72_100" in methods_group:
        label = "kraken_72_100"
        kraken_72_100_stats_dict = comp_stats(br_info_dict, kraken_72_100_info_dict, stats_dir_path, label)
    # Kraken2 non-gut
    if "k8_ng" in methods_group:
        label = "kraken_8_ng"
        kraken_8_ng_stats_dict = comp_stats(br_info_dict, kraken_8_ng_info_dict, stats_dir_path, label)
    if "k16_ng" in methods_group:
        label = "kraken_16_ng"
        kraken_16_ng_stats_dict = comp_stats(br_info_dict, kraken_16_ng_info_dict, stats_dir_path, label)
    if "k72_ng" in methods_group:
        label = "kraken_72_ng"
        kraken_72_ng_stats_dict = comp_stats(br_info_dict, kraken_72_ng_info_dict, stats_dir_path, label)
    # Kraken2 1%
    # COMEBin
    if "cnr" in methods_group:
        label = "comebine_nr"
        cmbn_nr_stats_dict = comp_stats(br_info_dict, cmbn_nr_info_dict, stats_dir_path, label)
    # MetaBinner
    if "mnr" in methods_group:
        label = "metabinner_nr"
        mtbr_nr_stats_dict = comp_stats(br_info_dict, mtbr_nr_info_dict, stats_dir_path, label)

    # Combine the information from all the dictionaries.
    comb_info_dict = {}
    # Kraken unfiltered
    if kraken_8_stats_dict:
        comb_info_dict["kraken_8"] = kraken_8_stats_dict
    if kraken_16_stats_dict:
        comb_info_dict["kraken_16"] = kraken_16_stats_dict
    if kraken_72_stats_dict:
        comb_info_dict["kraken_72"] = kraken_72_stats_dict
    # Kraken 0.1%
    if kraken_8_0c1_stats_dict:
        comb_info_dict["kraken_8_0c1"] = kraken_8_0c1_stats_dict
    if kraken_16_0c1_stats_dict:
        comb_info_dict["kraken_16_0c1"] = kraken_16_0c1_stats_dict
    if kraken_72_0c1_stats_dict:
        comb_info_dict["kraken_72_0c1"] = kraken_72_0c1_stats_dict
    # Kraken 1.0%
    if kraken_8_1c0_stats_dict:
        comb_info_dict["kraken_8_1c0"] = kraken_8_1c0_stats_dict
    if kraken_16_1c0_stats_dict:
        comb_info_dict["kraken_16_1c0"] = kraken_16_1c0_stats_dict
    if kraken_72_1c0_stats_dict:
        comb_info_dict["kraken_72_1c0"] = kraken_72_1c0_stats_dict
    # Kraken 10
    if kraken_8_10_stats_dict:
        comb_info_dict["kraken_8_10"] = kraken_8_10_stats_dict
    if kraken_16_10_stats_dict:
        comb_info_dict["kraken_16_10"] = kraken_16_10_stats_dict
    if kraken_72_10_stats_dict:
        comb_info_dict["kraken_72_10"] = kraken_72_10_stats_dict
    # Kraken 100
    if kraken_8_100_stats_dict:
        comb_info_dict["kraken_8_100"] = kraken_8_100_stats_dict
    if kraken_16_100_stats_dict:
        comb_info_dict["kraken_16_100"] = kraken_16_100_stats_dict
    if kraken_72_100_stats_dict:
        comb_info_dict["kraken_72_100"] = kraken_72_100_stats_dict
    # Kraken non-gut
    if kraken_8_ng_stats_dict:
        comb_info_dict["kraken_8_ng"] = kraken_8_ng_stats_dict
    if kraken_16_ng_stats_dict:
        comb_info_dict["kraken_16_ng"] = kraken_16_ng_stats_dict
    if kraken_72_ng_stats_dict:
        comb_info_dict["kraken_72_ng"] = kraken_72_ng_stats_dict
    # COMEBin nr
    if cmbn_nr_stats_dict:
        comb_info_dict["comebin_nr"] = cmbn_nr_stats_dict
    # MetaBinner nr
    if mtbr_nr_stats_dict:
        comb_info_dict["metabinner_nr"] = mtbr_nr_stats_dict

    # Create plots
    # Metric groups
    metric_group_1 = ["True Positive (TP)", "False Positive (FP)", "False Negative (FN)"]
    metric_group_2 = ["Sensitivity", "Precision", "Accuracy"]
    metric_group_3 = ["F1 score", "Jaccard Index", "Sorensen-Dice Index", "Bray-Curtis Dissimilarity"]
    metric_group_4 = ["MAE", "MSE", "RMSE"]
    metric_group_dict = {
        1: metric_group_1,
        2: metric_group_2,
        3: metric_group_3,
        4: metric_group_4
    }
    metric_label_dict = {
        "Gold species number": "Standard Species Abundance",
        "Predicted species number": "Predicted Species Abundance",
        "Common species (intersection): ": "Common Species Abundance",
        "Unique species from both groups (union)": "Unique Stdandard and Precietd Species Abundance",
        "Unique to gold group": "Unique Standard Species Abundance",
        "Unique to predicted group": "Unique Predicted Species Abundance",
        "True Positive (TP)": "True Positive",
        "True Negative (TN)": "True Negative",
        "False Positive (FP)": "False Positive",
        "False Negative (FN)": "False Negative",
        "Sensitivity": "Sensitivity",
        "Specificty": "Specificity",
        "Precision": "Precision",
        "Accuracy": "Accuracy",
        "F1 score": "F1 Score",
        "Jaccard Index": "Jaccard Index",
        "Sorensen-Dice Index": "Sorensen-Dice Index",
        "Bray-Curtis Dissimilarity": "Bray-Curtis Dissimilarity",
        "Sum of absolute differences": "Sum of Absolute Differences",
        "MAE": "Mean Absolute Error",
        "Sum of squared absolute differencies": "Sum of Squared Absolute ifferencies",
        "MSE": "Mean Squared Error",
        "RMSE": "Root Mean Squared Error"
    }
    legend_label_dict = {
        "kraken_8": "Kraken2 db:8",
        "kraken 16": "Kraken2 db:16",
        "kraken 72": "Kraken2 db:72",
        "kraken_8_0c1": "Kraken2 db:8 0.1%",
        "kraken_16_0c1": "Kraken2 db:16",
        "kraken_72_0c1": "Kraken2 db:72",
        "kraken_8_1c0": "Kraken2 db:8 1.0%",
        "kraken_16_1c0": "Kraken2 db:16 1.0%",
        "kraken_72_1c0": "Kraken2 db:72 1.0%",
        "kraken_8_10": "Kraken2 db:8 10",
        "kraken_16_10": "Kraken2 db:16 10",
        "kraken_72_10": "Kraken2 db:72 10",
        "kraken_8_100": "Kraken2 db:8 100",
        "kraken_16_100": "Kraken2 db:16 100",
        "kraken_72_100": "Kraken2 db:72 100",
        "kraken_8_ng": "Kraken2 db:8 non-gut",
        "kraken_16_ng": "Kraken2 db:8 non-gut",
        "kraken_72_ng": "Kraken2 db:8 non-gut",
        "comebin_nr": "COMEBin db:nr",
        "metabinner_nr": "MetaBinner db:nr"
    }
    # Sample groups.
    sample_groups_1 = [2, 8, 9, 16, 17, 18, 19]
    sample_groups_2 = [8, 18, 19]
    sample_groups_3 = [2, 9, 16, 17]
    sample_groups_4 = [5, 6, 11]
    sample_groups_5 = [7, 13, 15]
    sample_groups_6 = [1, 4, 10]
    sample_groups_7 = [3, 12, 14]
    sample_groups_8 = [2, 7, 9, 10, 11, 14, 16, 17, 19]
    sample_groups_9 = [7, 10, 11, 14, 19]
    sample_groups_10 = [2, 16, 17]
    sample_groups_11 = [1, 5, 12, 13, 18]
    sample_groups_12 = [3, 4, 6, 8, 15]
    sample_group_dict = {
        1: sample_groups_1,
        2: sample_groups_2,
        3: sample_groups_3,
        4: sample_groups_4,
        5: sample_groups_5,
        6: sample_groups_6,
        7: sample_groups_7,
        8: sample_groups_8,
        9: sample_groups_9,
        10: sample_groups_10,
        11: sample_groups_11,
        12: sample_groups_12
    }
    # Group labels
    sample_group_label_dict = {
        1: "10 Species",
        2: "10 Species - Simulated",
        3: "10 Species - Cultured",
        4: "40 Species",
        5: "120 Species",
        6: "500 Species",
        7: "1000 Species",
        8: "No bias",
        9: "No bias - Simulated",
        10: "No bias - Cultured",
        11: "AT-rich bias",
        12: "GC-rich bias"
    }
    pandas_dict = design_stacked_plots(comb_info_dict, metric_label_dict, legend_label_dict, plot_dir_parh)
    design_metric_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict, legend_label_dict, plot_dir_parh)
    design_metric_sample_grouped_plots(pandas_dict, metric_group_dict, metric_label_dict, legend_label_dict, sample_group_dict, sample_group_label_dict, plot_dir_parh)

    if time_status:
        # Time groups
        time_sample_groups_1 = [8, 18, 19]
        time_sample_groups_2 = [5, 6, 11]
        time_sample_groups_3 = [7, 13, 15]
        time_sample_groups_4 = [1, 4, 10]
        time_sample_groups_5 = [3, 12, 14]
        time_sample_group_dict = {
            1: time_sample_groups_1,
            2: time_sample_groups_2,
            3: time_sample_groups_3,
            4: time_sample_groups_4,
            5: time_sample_groups_5
        }
        time_sample_group_label_dict = {
            1: "10 Species",
            2: "40 Species",
            3: "120 Species",
            4: "500 Species",
            5: "1000 Species",
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
        time_analysis(time_dir, stats_dir_path, time_sample_group_dict, time_sample_group_label_dict, sample_size_dict, methods_time_group, sample_speciesnum_dict)


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
