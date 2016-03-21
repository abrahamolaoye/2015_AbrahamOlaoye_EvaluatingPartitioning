__author__ = 'Abraham'
from constants import *


def define_data_blocks(alignment_path):
    # defining blocks of sites for PF
    # this method works for alignments in the Nexus format.
    with open(r'%s' % alignment_path, 'r') as alignment_file:
        line = ""
        while line != "begin sets;\n":
            line = alignment_file.readline()
        blocks = ""
        for j in range(1, sample_size + 1):
            line = alignment_file.readline()
            blocks += line.replace("charset ", "").replace(".fa.out.nex", "").replace(".", "")
    return blocks


def create_configuration_files(alignment_set):
    print("Creating configuration files for each sample...")
    for sample_no in alignment_set.keys():
        with open(r'./OGs/Sample_%i/partition_finder.cfg' % sample_no, 'w') as new_config:
            filename = alignment_set[sample_no].replace("./OGs/Sample_%i/" % sample_no, "")
            model = "WAG+I+G"
            new_config.write("# ALIGNMENT FILE #\n" +
                             "alignment = %s;\n" % filename +
                             "\n"
                             "#BRANCHLENGTHS: linked | unlinked #\n"
                             "branchlengths = linked;\n"
                             "\n"
                             "# MODELS OF EVOLUTION #\n"
                             "models = %s;\n" % model +
                             "\n"
                             "# MODEL SELECTION #\n"
                             "model_selection = BIC;\n"
                             "\n"
                             "# DATA BLOCKS #\n"
                             "[data_blocks]\n"
                             "\n"
                             "%s\n" % define_data_blocks(alignment_set[sample_no].replace(".phy", ".nex")) +
                             "# SCHEMES #\n"
                             "[schemes]\n"
                             "search = greedy;\n")
