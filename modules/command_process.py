import subprocess


def command_run(phrase_1, phrase_2, title_1, title_2, capture_status, shell_status, pr_status, input_log_file, output_log_file):
    # Process
    input_log_file.write("{}\n".format(title_1))
    if pr_status:
        print("\n{}".format(phrase_1))
    input_log_file.write("{}\n\n".format(phrase_1))
    if shell_status:
        command = phrase_1
    else:
        command = phrase_1.split(" ")
    proc = subprocess.run(command, capture_output=capture_status, shell=shell_status)
    stdout_raw = proc.stdout
    stdout_str = stdout_raw.decode("utf-8")
    stdout_lines = stdout_str.split("\n")
    for line in stdout_lines:
        output_log_file.write("{}\n".format(line))
    # Version, if available through the command-line tool.
    input_log_file.write("{}\n".format(title_2))
    if phrase_2 and phrase_2 != "":
        input_log_file.write("{}\n\n".format(phrase_2))
        if shell_status:
            command = phrase_2
        else:
            command = phrase_2.split(" ")
        proc = subprocess.run(command, capture_output=capture_status, shell=shell_status)
        stdout_raw = proc.stdout
        stdout_str = stdout_raw.decode("utf-8")
        stdout_lines = stdout_str.split("\n")
        output_log_file.write("Version information:\n")
        for line in stdout_lines:
            output_log_file.write("{}\n".format(line))
    output_log_file.write("\n{}\n\n\n".format(100*"-"))