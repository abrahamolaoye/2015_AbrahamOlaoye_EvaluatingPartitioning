__author__ = 'Abraham'

from glob import *
from random import choice
from constants import *


def collect_list():
    pathname = script_dir + 'OGs/*.fa.out'
    file_addresses = [f for f in glob(pathname)]  # file addresses
    file_list = []
    if len(file_addresses) > 0:
        for k in range(1, (no_of_samples * sample_size) + 1):
            list_index = choice(range(0, len(file_addresses)))
            file_chosen = file_addresses[list_index]
            file_list.append(file_chosen)

    return file_list


def find_correct_sample(filename):
    index = list_of_files.index(filename) + 1
    sample_no = 0
    if index % sample_size == 0:
        sample_no = index // sample_size
    else:
        sample_no = (index // sample_size) + 1
    return sample_no


def muscle_job(header_txt):
    filename = "alignment_to_partitions.sh"
    jobs = ""
    for input_seq in list_of_files:
        sample_no = find_correct_sample(input_seq)  # need to identify correct folder to put output into
        output_seqs = input_seq.replace("OGs", "OGs/Sample_%i" % sample_no)
        command = "muscle -in %s -out %s" % (input_seq, output_seqs)  # run muscle on file, store in folder
        jobs += command + "\n"

    with open(r'%s' % filename, 'w') as newJob:
        newJob.write(header_txt + jobs)


muscle_job()
