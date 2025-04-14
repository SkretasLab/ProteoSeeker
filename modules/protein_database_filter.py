import os
import sys
import copy
import math
import time
import multiprocessing


def utf8len(str_count):
    return len(str_count.encode('utf-8'))


def read_chunk(pd_path, proc_id, byte_start, byte_end, file_size, db_info_dict, return_dict):
    cur_line = "Empty"
    first_pass = True
    byte_count = byte_start
    last_line_size = 0
    allow_end = True
    first_overhead = True
    for key in db_info_dict.keys():
        pnr_db_file_name = db_info_dict[key][0]
        pnr_db_file_name = "{}_part_{}.fasta".format(pnr_db_file_name, proc_id)
        new_file = open(pnr_db_file_name, "w", encoding="utf-8")
        db_info_dict[key][3] = new_file
    with open(pd_path, encoding='utf-8') as fh:
        while (byte_count < byte_end) or (cur_line[0] != ">") or (not allow_end):
            if first_pass:
                fh.seek(byte_start)
                if byte_start != 0:
                    pr_char = fh.read(1)
                    byte_count += 1
                    if byte_count == file_size:
                        break
                    cur_char = fh.read(1)
                    byte_count += 1
                    if byte_count == file_size:
                        break
                    while cur_char != ">" or pr_char != "\n":
                        pr_char = copy.deepcopy(cur_char)
                        cur_char = fh.read(1)
                        byte_count += 1
                        if byte_count == file_size:
                            break
                    cur_line_pre = fh.readline()
                    cur_line = ">{}".format(cur_line_pre)
                    mod_byte_start = byte_count
                    byte_count += utf8len(cur_line_pre)
                    print("Modified starting byte: {}".format(mod_byte_start))
                    return_dict[proc_id][0] = mod_byte_start
                else:
                    mod_byte_start = byte_start
                    cur_line = fh.readline()
                    byte_count += utf8len(cur_line)
                    print("Modified starting byte: {}".format(mod_byte_start))
                    return_dict[proc_id][0] = mod_byte_start
                    if byte_count == file_size:
                        break
                first_pass = False
            else:
                cur_line = fh.readline()
                byte_count += utf8len(cur_line)
                if byte_count == file_size:
                    break
            if (byte_count >= byte_end) and (cur_line[0] == ">") and (not allow_end):
                allow_end = True
            if (byte_count >= byte_end) and first_overhead:
                allow_end = False
                first_overhead = False
            if cur_line[0] == ">":
                for key in db_info_dict.keys():
                    db_info_dict[key][2] = False
                    target_names = db_info_dict[key][1]
                    for item in target_names:
                        if item in cur_line:
                            db_info_dict[key][2] = True
                            break
            for key in db_info_dict.keys():
                if db_info_dict[key][2]:
                    new_file = db_info_dict[key][3]
                    new_file.write("{}".format(cur_line))
        last_line_size = utf8len(cur_line)
        # The last line where the process ends must be removed from the byte count, EXCEPT in the case of EOF.
        if byte_count == file_size:
            mod_byte_end = byte_count
        else:
            mod_byte_end = byte_count - last_line_size
        print("Modified ending byte: {}".format(mod_byte_end))
        return_dict[proc_id][1] = mod_byte_end
    for key in db_info_dict.keys():
        new_file = db_info_dict[key][3]
        new_file.close()


def prfilter(pd_path="", info_file_path="", db_info_dict={}, log_file_path="", thread_num=2):
    start_time = time.time()
    
    if not log_file_path:
        log_file_path = "pd_filter_log.txt"
    log_file = open(log_file_path, "w")

    if info_file_path:
        if os.path.exists(info_file_path):
            r = 0
            info_file = open(info_file_path, "r")
            info_lines = info_file.readlines()
            info_file.close()
            db_path_line = True
            db_info_dict[r] = [None, [], False, None]
            for line in info_lines:
                line = line.rstrip()
                if db_path_line:
                    db_info_dict[r][0] = line
                    db_path_line = False
                elif line:
                    db_info_dict[r][1].append(line)
                else:
                    r += 1
                    db_path_line = True
        else:
            print("\nThe input info path was not found. Exiting.")
            log_file.write("The input info path was not found. Exiting.\n")
            exit()
    else:
        if not db_info_dict:
            print("\nNo input info dictionary was provided. Exiting.")
            log_file.write("No input info dictionary was provided. Exiting.\n")
            exit()

    # nr = 363213505071 (06/03/2024)
    file_size =  os.path.getsize(pd_path)
    chunk = file_size / thread_num
    chunk_size = math.ceil(chunk)
    first_iteration = True
    byte_end = -1
    starting_points = []
    ending_points = []
    while byte_end != file_size:
        # Make sure that the ending line is the last line a of protein sequence.
        if first_iteration:
            byte_start = 0
            byte_end = chunk_size
            first_iteration = False
        else:
            byte_start += chunk_size
            byte_end += chunk_size
        if byte_end > file_size:
            byte_end = file_size
        starting_points.append(byte_start)
        ending_points.append(byte_end)

    chunk_num = len(starting_points)
    file_size_gb = file_size / 1073741824
    file_size_gb = round(file_size_gb , 2)
    chunk_size_gb = chunk_size / 1073741824
    chunk_size_gb = round(chunk_size_gb , 2)
    print("File_size: {}".format(file_size))
    print("File_size: {} GB".format(file_size_gb))
    print("Chunk size: {} B".format(chunk_size))
    print("Chunk size: {} GB".format(chunk_size_gb))
    print("Number of chunks: {}".format(chunk_num))
    print("Starting bytes:")
    print(starting_points)
    print("Ending bytes:")
    print(ending_points)
    log_file.write("File_size: {}\n".format(file_size))
    log_file.write("File_size: {} GB\n".format(file_size_gb))
    log_file.write("Chunk size: {} B\n".format(chunk_size))
    log_file.write("Chunk size: {} GB\n".format(chunk_size_gb))
    log_file.write("Number of chunks: {}\n".format(chunk_num))
    log_file.write("Starting bytes: {}\n".format(starting_points))
    log_file.write("Ending bytes: {}\n".format(ending_points))


    # A dictinoary to return values to write in the log file.
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    # Creating the processes.
    db_limit = len(starting_points)
    proc_group  = []
    for point_index in range(0, db_limit):
        byte_start = starting_points[point_index]
        byte_end = ending_points[point_index]
        return_dict[point_index] = manager.list([None, None])
        proc_args = (pd_path, point_index, byte_start, byte_end, file_size, db_info_dict, return_dict)
        proc_id = multiprocessing.Process(target=read_chunk, args=proc_args)
        proc_group.append(proc_id)

    # Starting the proceses and waiting for them to finish.
    # 3600 seconds = 1 hour
    proc_wait_time = 3600
    for proc_index in range(0, db_limit):
        proc_id = proc_group[proc_index]
        # Start process
        print("Starting process: {}".format(proc_index))
        log_file.write("Starting process: {}\n".format(proc_index))
        proc_id.start()
    for proc_index in range(0, db_limit):
        proc_id = proc_group[proc_index]
        # There is no problem if the process to be joined has already finished.
        proc_id.join(proc_wait_time)
        print("Ended process: {}".format(proc_index))
        print()
        log_file.write("Ended process: {}\n".format(proc_index))
        # If the process is still alive, it means that it did not finish in the specified period of time and thus it has to be terminated.
        proc_status = proc_id.is_alive()
        if proc_status:
            print("Process is still running. Terminating process with ID {}".format(proc_id))
            proc_id.terminate()

    # Write information from the processes in the log file.
    for key_proc in return_dict.keys():
        log_file.write("Modified starting byte: {}\n".format(return_dict[key_proc][0]))
        log_file.write("Modified ending byte: {}\n".format(return_dict[key_proc][1]))

    end_time = time.time()
    elpased_time = end_time - start_time
    elpased_time_min = elpased_time / 60
    elpased_time_hour = elpased_time_min / 60
    elpased_time = round(elpased_time, 3)
    elpased_time_min = round(elpased_time_min, 3)
    elpased_time_hour = round(elpased_time_hour, 3)
    print("Time elapsed: {} seconds / {} minutes / {} hours".format(elpased_time, elpased_time_min, elpased_time_hour))
    log_file.write("Time elapsed: {} seconds / {} minutes / {} hours\n".format(elpased_time, elpased_time_min, elpased_time_hour))

    log_file.close()


if __name__ == "__main__":
    arg_pd_path = ""
    arg_info_file_path = ""
    arg_db_info_dict = {}
    arg_log_file_path = ""
    arg_thread_num = 2
    if len(sys.argv) > 1:
        arg_input_command = "{}".format(sys.argv[0])
        for i in range(1, len(sys.argv), 2):
            if sys.argv[i] == "-i" or sys.argv[i] == "--input-db":
                arg_pd_path = sys.argv[i+1]
            elif sys.argv[i] == "-b" or sys.argv[i] == "--input-info":
                arg_info_file_path = sys.argv[i+1]
            elif sys.argv[i] == "-l" or sys.argv[i] == "--log-file":
                arg_log_file_path = sys.argv[i+1]
            elif sys.argv[i] == "-t" or sys.argv[i] == "--threads":
                arg_thread_num = int(sys.argv[i+1])
    prfilter(arg_pd_path, arg_info_file_path, arg_db_info_dict, arg_log_file_path, arg_thread_num)
