__author__ = 'Abraham'
# This script requires as input a folder of aligned OGs, and a number of folders in which the results
# for each sample will be stored (the number is specified in 'constants.py'
# The main outputs are the trees produced by the programs you wish to run, along with the monophyly results for each


from job_file_methods import *


def choose_which_to_run(progs):
    for prog_name in progs.keys():
        if progs[prog_name] == "":
            run_prog = ""
            run_prog = raw_input('Run %s? (y/n): ' % prog_name).lower()

            while run_prog != 'y' and run_prog != 'n':
                run_prog = raw_input('%s? (y/n): ' % prog_name).lower()

            progs[prog_name] = run_prog

            if prog_name == "PartitionFinder" and run_prog == 'n':
                print "You won't be able to run RaxML on PF partitions since PartitionFinder won't run."
                progs["RaxML on PF partitions"] = 'n'
            elif prog_name == "MtPAN" and run_prog == 'n':
                print "You won't be able to run RaxML on MtPAN partitions since MtPAN won't run."
                progs["RaxML on MtPAN partitions"] = 'n'

    return progs


def copy_files_into_samples():
    tasks = ""
    for k in range(1, no_of_samples + 1):
        for og in list_of_files[k-1]:
            tasks += "cp %s %s\n" % (og, alignments_dir + "Sample_%i/" % k)

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


def create_configuration_file(phy_alignment, sample_no, folder, option):
    print("Creating configuration files for %s..." % phy_alignment)
    with open(r'%spartition_finder.cfg' % folder, 'w') as new_config:
        filename = phy_alignment.replace(alignments_dir + "Sample_%i/" % sample_no, "")
        model = "WAG+I+G"
        scheme = option
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
                         "%s\n" % define_data_blocks(phy_alignment.replace(".phy", ".nex")) +
                         "# SCHEMES #\n"
                         "[schemes]\n"
                         "search = %s;\n" % scheme)

programs = {"PartitionFinder": "", "MtPAN": "", "RaxML vanilla": "", "RaxML on PF partitions": "", "RaxML on MtPAN partitions": "",
            "TreeCollection": "", "PhyloBayes": "", "partitioned-by-gene RaxML": "", "concatenation": "",
            "NEXUS to PHYLIP conversion": "", "create new samples": "", "Raxml on PF partitions (relaxed-clustering)": ""}


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
    sample_dirs = alignments_dir + "Sample_%i/" % i
    make_dir(sample_dirs)
    nex_alignment = sample_dirs + "sample_alignment.nex"
    nex_dict.update({i: nex_alignment})
    phylip_alignment_path = sample_dirs + "sample_alignment.phy"
    alignment_dict.update({i: phylip_alignment_path})

move_tasks = ""
if programs["create new samples"] == "y":
    move_tasks = copy_files_into_samples()

concat_tasks = ""
# convert to PHYLIP-sequential format
conversion_tasks = ""
og_list_task = ""
job_filename = script_dir + "job_files/alignment_to_conversion.sh"
if programs["TreeCollection"] == "y":
    og_list_task = "python make_OG_lists.py\n"

if programs["concatenation"] == "y":
    for m in range(1, no_of_samples + 1):
        concat_tasks += "python concatenate.py %i\n" % m
    # convert to PHYLIP-sequential format

if programs["NEXUS to PHYLIP conversion"] == "y":
    conversion_tasks = write_conversion_tasks(nex_dict)

if programs["concatenation"] == "y" or programs["create new samples"] == "y" or programs["NEXUS to PHYLIP conversion"] == "y":
    with open(r'%s' % job_filename, 'w') as newJob:
        job_resources = res_text_sge("1G", 100, 1)  # reserved memory, number of hours and number of threads.
        newJob.write(job_resources + move_tasks + concat_tasks + conversion_tasks + og_list_task)

    run_job(job_filename, True)


# RaxML w/o partitions
if programs["RaxML vanilla"] == "y":
    vanilla_raxml_job(alignment_dict, 12)  # dictionary of alignments, number of threads

if programs["partitioned-by-gene RaxML"] == "y":
    partitioned_by_gene_raxml(alignment_dict, nex_dict, 12)


if programs["Raxml on PF partitions (relaxed-clustering)"] == "y":
    for j in range(1, no_of_samples + 1):
        analysis_folder = alignments_dir + "Sample_%i/PF_analysis_RC/" % j  # RC = relaxed-clustering method for PF
        make_dir(analysis_folder)
        search_op_rc = "rcluster"
        pf_rc_tasks = write_pf_task(j, analysis_folder, search_op_rc)
        create_configuration_file(alignment_dict[j], j, analysis_folder, search_op_rc)
        part_dir = analysis_folder + "partitions.txt"
        pf_rc_raxml = write_pf_rc_raxml_task(alignment_dict, j, 12)
        job_filename = script_dir + "job_files/PF_RC_raxml_S%i.sh"% j
        with open(r'%s' % job_filename, 'w') as job:
            resources = res_text_sge("2G", 1000, 12)
            job.write(resources + pf_rc_tasks + pf_rc_raxml)
        run_job(job_filename, False)


if programs["RaxML on PF partitions"] == "y":
    for j in range(1, no_of_samples + 1):
        analysis_folder = alignments_dir + "Sample_%i/" % j
        search_op = "greedy"
        pf_task = write_pf_task(j, analysis_folder, search_op)
        create_configuration_file(alignment_dict[j], j, analysis_folder, search_op)
        raxml_part_file_dir = analysis_folder + "partitions.txt"
        # make RaxML-style partition file from PF schemes
        pf_raxml_task = write_pf_raxml_task(alignment_dict, j, 12)  # dict of PHYLIP alignments, sample number, no of threads
        job_filename = script_dir + "job_files/PF_raxml_S%i.sh" % j
        with open(r'%s' % job_filename, 'w') as job:
            resources = res_text_sge("2G", 1000, 12)
            job.write(resources + pf_task + pf_raxml_task)

        run_job(job_filename, False)

elif programs["PartitionFinder"] == "y":
    for j in range(1, no_of_samples + 1):
        pf_task = write_pf_task(j)
        job_filename = script_dir + "job_files/PF_S%i.sh" % j
        with open(r'%s' % job_filename, 'w') as job:
            resources = res_text_sge("2G", 1000, 12)
            job.write(resources + pf_task)

        run_job(job_filename, False)

if programs["RaxML on MtPAN partitions"] == "y":
    for j in range(1, no_of_samples + 1):
            mtpan_task = write_mtpan_task(nex_dict, j)
            mtpan_raxml_tasks = mtpan_raxml(j, alignment_dict, 4, 8, 12)
            # make RaxML-style partition file from MtPAN clusters
            clusters_to_partitions_task = "python clusters_to_partitions.py %i\n" % j
            job_filename = script_dir + "job_files/mtpan_raxml_S%i.sh" % j
            with open(r'%s' % job_filename, 'w') as job:
                resources = res_text_sge("2G", 1000, 12)
                job.write(resources + mtpan_task + clusters_to_partitions_task + mtpan_raxml_tasks)

            run_job(job_filename, False)

elif programs["MtPAN"] == "y":
        for j in range(1, no_of_samples + 1):
            mtpan_task = write_mtpan_task(nex_dict, j)
            job_filename = script_dir + "job_files/MtPAN_clustering_S%i.sh" % j
            with open(r'%s' % job_filename, 'w') as job:
                resources = res_text_sge("2G", 1000, 12)
                job.write(resources + mtpan_task)

            run_job(job_filename, False)

if programs["PhyloBayes"] == "y":
    # Run PhyloBayes on samples
    phylobayes_job(alignment_dict, 12)

if programs["TreeCollection"] == "y":
    # Run TreeCollection
    tree_coll_job(4, nex_dict)
