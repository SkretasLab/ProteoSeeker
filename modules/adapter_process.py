import os
import shutil
# ProteoSeeker modules
import supportive_functions


def adapter_file(output_path_fastqc, adapters_status, adapters_path, adapters_ap_path):
    print("\nAnalyzing the file with the adapters...")
    # adapters_status
    # identify:The overrepresented sequences are added to the already existing adatpers to be searched.
    # fastqc: Only the overrepresented sequences are used as adatpers.
    # Multiple files can be given as inputs, so ultiple FastQC HMTL files may exist, thus multiple adapters can be found by each of these files. That's why the counter is initiazlied here.
    # When overrepresented sequences are not found then the sentence in the HTML file is a header.
    # When overrepresented sequences are found then these sequences are shown in a table.
    over_sequences = []
    if (adapters_status == "ide") or (adapters_status == "fas"):
        fastqc_files = os.listdir(output_path_fastqc)
        for fqc_file in fastqc_files:
            if fqc_file[-5:] == ".html":
                fqc_lines = supportive_functions.read_file(output_path_fastqc + fqc_file)
                for fqc_line in fqc_lines:
                    if "Overrepresented sequences</h2><p>" in fqc_line:
                        splited_fq_1 = fqc_line.split("Overrepresented sequences</h2><p>")
                        splited_fq_2 = splited_fq_1[1].split("</p>")
                        over_info = splited_fq_2[0]
                        if over_info == "No overrepresented sequences":
                            print("No overrepresented sequences for the reads of the file: {}{}.\n".format(output_path_fastqc, fqc_file))
                    if "Overrepresented sequences</h2><table>" in fqc_line:
                        splited_fq_1 = fqc_line.split("Overrepresented sequences</h2><table>")
                        splited_fq_2 = splited_fq_1[1].split("<tr><td>")
                        for ov_seq_info in splited_fq_2[1:]:
                            over_seq_pre = ov_seq_info.split("</td><td>")[0]
                            over_sequences.append(over_seq_pre)
                            if "</table>" in ov_seq_info:
                                break
        # If no overrepresented sequences were found and "ide" was selected, then the adapters file will contain the already known adapters.
        # If no overrepresented sequences were found and "fas" was selected, then the adapters file will contain no sequences (empty file).
        if over_sequences:
            over_seq_counter = 1
            if adapters_status == "ide":
                shutil.copy(adapters_path, adapters_ap_path)
                adapters_ap_file = open(adapters_ap_path, "a")
                for over_seq in over_sequences:
                    adapters_ap_file.write(">Overrepresented_Sequence_{}\n".format(over_seq_counter))
                    adapters_ap_file.write("{}\n".format(over_seq))
                    over_seq_counter += 1
                adapters_ap_file.close()
            elif adapters_status == "fas":
                adapters_ap_file = open(adapters_ap_path, "w")
                for over_seq in over_sequences:
                    adapters_ap_file.write(">Overrepresented_Sequence_{}\n".format(over_seq_counter))
                    adapters_ap_file.write("{}\n".format(over_seq))
                    over_seq_counter += 1
                adapters_ap_file.close()
    elif adapters_status == "pre":
        shutil.copy(adapters_path, adapters_ap_path)
    else:
        shutil.copy(adapters_path, adapters_ap_path)