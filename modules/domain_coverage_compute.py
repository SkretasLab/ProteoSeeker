import os
# ProteoSeeker modules
import supportive_functions


def dom_percs(comb_all_domains_proteins, hmmer_profile_lengths_file, hmmer_common_lengths_file):
    print("\nComputing profile coverages...")
    if os.path.exists(comb_all_domains_proteins):
        profiles_lengths_lines = supportive_functions.read_file(hmmer_profile_lengths_file)
        profiles_codes_lengths_dict = {}
        # Collect the length of each profile.
        for line in profiles_lengths_lines:
            line_splited = line.split("\t")
            profile_code = line_splited[1]
            profile_length = int(line_splited[2])
            if profile_code in profiles_codes_lengths_dict.keys():
                print("Profile code duplicate from lengths info found!")
                exit()
            else:
                profiles_codes_lengths_dict[profile_code] = profile_length
        new_file = open(hmmer_common_lengths_file, "w")
        # Parse each protein from the file of combined proteins.
        domains_info_lines = supportive_functions.read_file(comb_all_domains_proteins)
        for i in range(1, len(domains_info_lines)):
            line = domains_info_lines[i]
            splited_line = line.split("\t")
            contig_name = splited_line[0]
            domain_code = splited_line[4]
            profile_common = int(splited_line[9]) - int(splited_line[8])
            # If for some reason the domain code is not found in the dictionary of profile lengths,
            # the profile coverage is undefined and it is written as -1.
            if domain_code in profiles_codes_lengths_dict.keys():
                profile_length_cur = profiles_codes_lengths_dict[domain_code]
                profile_coverage = profile_common / profile_length_cur
                profile_coverage = round(profile_coverage, 2)
                profile_coverage_perc = profile_coverage * 100
                profile_coverage_perc = int(profile_coverage_perc)
            else:
                profile_coverage_perc = -1
            new_file.write("{}\t{}%\n".format(contig_name, profile_coverage_perc))
        new_file.close()