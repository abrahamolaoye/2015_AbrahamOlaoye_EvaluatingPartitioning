__author__ = 'Abraham'
from random import *
from glob import *
# specifies the essential variables for use of the pipeline.
sample_size = 40
no_of_samples = 20
script_dir = "/SAN/biosciences/oma/research/Abraham_Olaoye/Pipeline/"
alignments_dir = script_dir + "Data/"


def collect_list():
    addresses = alignments_dir + '*.fas'
    file_addresses = [f for f in glob(addresses)]  # file addresses
    file_list = []
    for m in range(1, no_of_samples + 1):
        file_list.append([])
        for n in range(1, sample_size + 1):
            list_index = choice(range(0, len(file_addresses)))
            file_chosen = file_addresses[list_index]
            while file_chosen in file_list[m - 1]:  # make sure no OG is repeated in a sample
                list_index = choice(range(0, len(file_addresses)))
                file_chosen = file_addresses[list_index]
            file_list[m - 1].append(file_chosen)

    return file_list

list_of_files = collect_list()
