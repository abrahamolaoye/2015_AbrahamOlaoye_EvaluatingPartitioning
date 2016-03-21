__author__ = 'Abraham'
# This script requires as input a folder of aligned OGs, and a number of folders in which the results
# for each sample will be stored (the number is specified in 'constants.py'
# The main outputs are the trees produced by the programs you wish to run, along with the monophyly results for each


from job_file_methods import *


def choose_which_to_run(progs):
    for name in progs.keys():
        if progs[name] == "":
            run_prog = ""
            run_prog = raw_input('Run %s? (y/n): ' % name).lower()

            while run_prog != 'y' and run_prog != 'n':
                run_prog = raw_input('Run %s? (y/n): ' % name).lower()

            progs[name] = run_prog

            if name == "PartitionFinder" and run_prog == 'n':
                print "You won't be able to run RaxML on PF partitions since PartitionFinder won't run."
                progs["RaxML on PF partitions"] = 'n'
            elif name == "MtPAN" and run_prog == 'n':
                print "You won't be able to run RaxML on MtPAN partitions since PartitionFinder won't run."
                progs["RaxML on MtPAN partitions"] = 'n'

    return progs


def move_files_into_samples():
    tasks = ""
    for k in range(1, no_of_samples + 1):
        for og in list_of_files[k-1]:
            tasks += "mv %s %s\n" % (og, script_dir + "OGs/Sample_%i/" % k)

    return tasks


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
            filename = alignment_set[sample_no].replace(script_dir + "OGs/Sample_%i/" % sample_no, "")
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

programs = {"PartitionFinder": "", "MtPAN": "", "RaxML vanilla": "", "RaxML on PF partitions": "", "RaxML on MtPAN partitions": "",
            "TreeCollection": "", "PhyloBayes": ""}
choose_which_to_run(programs)
#  file alignment
list_of_files = collect_list()  # list of files that will be grouped into samples
# for each sample, combine the alignments into one concatenated alignment, then generate
# a configuration file for it so it may be used by PartitionFinder
# also create list of OGs whose alignment is in the sample. This is needed by Tree collection

nex_dict = {}   # each element will have the 'sample number' as the key and the nexus alignment file-path as the value.
alignment_dict = {}  # each element will have the 'sample number' as the key and the
                     # PHYLIP alignment file-path as the value.

for i in range(1, no_of_samples + 1):
    nex_alignment = script_dir + "OGs/Sample_%i/sample_alignment_%i.nex" % (i, i)
    nex_dict.update({i: nex_alignment})
    phylip_alignment_path = script_dir + "OGs/Sample_%i/sample_alignment_%i.phy" % (i, i)
    alignment_dict.update({i: phylip_alignment_path})

job_filename = script_dir + "job_files/alignment_to_conversion.sh"
move_tasks = move_files_into_samples()
concat_tasks = "python concatenate.py\n"
# convert to PHYLIP-sequential format
conversion_tasks = write_conversion_tasks(nex_dict)
og_list_task = ""
if programs["TreeCollection"] == "y":
    og_list_task = "python make_OG_lists.py\n"

with open(r'%s' % job_filename, 'w') as newJob:
    job_resources = res_text("1G", 100, 1)  # reserved memory, number of hours and number of threads.
    newJob.write(job_resources + move_tasks + concat_tasks + conversion_tasks + og_list_task)

run_job(job_filename, True)


# RaxML w/o partitions
if programs["RaxML vanilla"] == "y":
    vanilla_raxml_job(alignment_dict, 10)  # dictionary of alignments, number of threads

# Make partitions using PartitionFinder, MtPAN, and then run RaxML on the resulting partitions.
if programs["PartitionFinder"] == "y":
    create_configuration_files(alignment_dict)


for j in range(1, no_of_samples + 1):
    pf_task = ""
    pf_raxml_task = ""
    mtpan_raxml_tasks = ""
    mtpan_task = ""
    partition_file_task = ""

    if programs["PartitionFinder"] == "y":
        pf_task = write_pf_task(j)
    if programs["RaxML on PF partitions"] == "y":
        raxml_part_file_dir = script_dir + "OGs/Sample_%i/sample_partitions.txt" % j
        # make RaxML-style partition file from PF schemes
        partition_file_task = "python make_partition_file.py %s %i\n" % (raxml_part_file_dir, j)
        pf_raxml_task = write_pf_raxml_task(alignment_dict, j, 12)  # dict of PHYLIP alignments, sample number, no of threads
    if programs["MtPAN"] == "y":
        mtpan_task = write_mtpan_task(nex_dict, j)

    if programs["RaxML on MtPAN partitions"] == "y":
        mtpan_raxml_tasks = mtpan_raxml(j, alignment_dict, 12)
        # make RaxML-style partition file from MtPAN clusters
        clusters_to_partitions_task = "python clusters_to_partitions.py %i\n" % j

    job_filename = script_dir + "job_files/partitioned_analyses_S%i.sh" % j
    with open(r'%s' % job_filename, 'w') as job:
        resources = res_text("1G", 1000, 12)
        job.write(resources + pf_task + partition_file_task + pf_raxml_task + mtpan_task + mtpan_raxml_tasks)

    run_job(job_filename, False)

if programs["PhyloBayes"] == "y":
    # Run PhyloBayes on samples
    phylobayes_job(alignment_dict)

if programs["TreeCollection"] == "y":
    # Run TreeCollection
    tree_coll_job()
