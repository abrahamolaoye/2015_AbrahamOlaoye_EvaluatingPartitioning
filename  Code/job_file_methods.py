__author__ = 'Abraham'
# a module containing methods relating to the creation and submission of jobs telling the cluster
# to execute a chosen program.
from shlex import *
from subprocess import *
from constants import *
from glob import *
from Bio.Phylo import *
from os import mkdir
from os import path

times_folder = script_dir + "times/"


def monophyly_task(tree_path, folder):
    clades_path = script_dir + "clade_list.txt"
    monophyly_script_path = script_dir + "external_tools/monophyly/max_monophyly.pl"
    task = "perl %s %s %s >> %s_monophyly.txt\n" % (monophyly_script_path, clades_path, tree_path, tree_path)
    task += 'echo "%s" >> %smon_results.txt\n' % (tree_path + ".tre_monophyly.txt", folder)

    return task


def write_pf_rc_raxml_task(dict_of_alignments, sample_no, no_threads):
    folder = alignments_dir + "Sample_%i/PF_analysis_RC/" % sample_no
    task = "cp %s %s\n" % (dict_of_alignments[sample_no], folder)  # need PHYLIP alignment to run PF
    partitions_path = folder + "sample_partitions.txt"  # name of partitions to be used by RaxML
    random_seed = choice(range(10000, 20000))
    pf_scheme_path = folder + "analysis/best_scheme.txt"
    task += "python make_partition_file.py %s %s\n" % (partitions_path, pf_scheme_path)
    substitution_model = "PROTGAMMAIWAG"
    task += time_code_beginning()
    params = (no_threads, substitution_model, random_seed, partitions_path, dict_of_alignments[sample_no], folder, "Tree")
    task += "raxmlHPC-PTHREADS -T %i -M -m %s -p %i -q %s -s %s -w %s -n %s\n" % params
    task += time_code_end(times_folder + "PF_RC_Raxml_time_S%i.txt" % sample_no)
    final_tree = folder + "RAxML_bestTree.Tree"
    task += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
    task += 'echo "Generated tree (RaxML with PF partitions (using relaxed-clustering)." >> S%i_progress.txt\n' % sample_no
    task += "python map_names_to_species.py %s %s\n" % (final_tree, dict_of_alignments[sample_no].replace(".phy",
                                                                                                          ".nex")
                                                        + "SpeciesMapping.txt")
    task += to_nexus_format_task(final_tree)
    task += monophyly_task(final_tree + ".tre", folder)
    task += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree, sample_no)

    return task


def partitioned_by_gene_raxml(dict_of_alignments, nex_alignments, no_threads):
    for sample_no in dict_of_alignments.keys():
        substitution_model = "PROTGAMMAIWAG"
        task = ""
        output_folder = alignments_dir + "Sample_%i/raxml_results/with_partitioning/PartitionByGene/" % sample_no
        make_dir(output_folder)
        random_seed = choice(range(10000, 20000))
        partitions_path = alignments_dir + "Sample_%i/by_gene_partition.txt" % sample_no  # name of partitions to be used by RaxML
        task += "python partition_by_genes.py %s %s\n" % (nex_alignments[sample_no], partitions_path)
        params = (no_threads, substitution_model, random_seed, partitions_path, dict_of_alignments[sample_no],
                  output_folder, "Tree")
        task += time_code_beginning()
        task += "raxmlHPC-PTHREADS -T %i -M -m %s -p %i -q %s -s %s -w %s -n %s\n" % params
        task += time_code_end(times_folder + "partitioned_by_gene_raxml_time_S%i.txt" % sample_no)
        final_tree = output_folder + "RAxML_bestTree.Tree"
        task += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
        task += 'echo "Generated tree (RaxML with partitions-by-gene)." >> S%i_progress.txt\n' % sample_no
        task += to_nexus_format_task(final_tree)
        task += monophyly_task(final_tree + ".tre", output_folder)
        task += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree, sample_no)
        job_filename = script_dir + "job_files/part_by_gene_raxml_S%i.sh" % sample_no
        with open(r'%s' % job_filename, 'w') as job:
            resources = res_text_sge("1G", 400, no_threads)
            job.write(resources + task)

        run_job(job_filename, False)


def tree_coll_job(no_threads, nex_alignments):
    for k in range(1, no_of_samples + 1):
        # making inputs for create_treecoll_input
        tree_coll_tasks = tc_tasks(k, nex_alignments)
        filename = script_dir + "job_files/tree_collection_S%i.sh" % k
        with open(r'%s' % filename, 'w') as newJob:
            tree_coll_resources = res_text_sge("1G", 500, no_threads)
            newJob.write(tree_coll_resources + tree_coll_tasks)

        run_job(filename, False)


def tc_tasks(sample_no, nex_alignments):
    tree_coll_tasks = ""
    output_directory = alignments_dir + "Sample_%i/tree_collection_results/" % sample_no
    make_dir(output_directory)
    og_list_address = output_directory + "OG_list.txt"
    tree_coll_tasks += "python %sexternal_tools/tree_collection/create_treecoll_input_aligned.py %s %s\n" % (script_dir, og_list_address, output_directory)
    tree_coll_tasks += 'echo "Created dist-var matrix, maps and labels." >> S%i_progress.txt\n' % sample_no
    # specify name of inputs for TreeCollection
    dist_var_name = output_directory + "DistVar.txt"
    map_name = output_directory + "Map.txt"
    label_name = output_directory + "Labels.txt"
    # make dist-var matrix made of weighted averages of pairwise distances
    weighted_distvar_name = output_directory + "tree_distvar"
    params = (script_dir, dist_var_name, map_name, label_name, weighted_distvar_name)
    tree_coll_tasks += "python %sexternal_tools/average_matrix_distvar_WEIGHTED.py %s %s %s %s\n" % params
    # make initial tree_topology from dist-var matrix using BIONJ
    output_tree_name = output_directory + "initial_tree.txt"
    tree_coll_tasks += "BIONJ %s.PHYLIP %s\n" % (weighted_distvar_name, output_tree_name)
    # Make bifurcating tree w/o negatives
    tree_coll_tasks += "python root_tree.py %s\n" % output_tree_name
    tree_coll_tasks += 'echo "Generated initial tree topology." >> S%i_progress.txt\n' % sample_no
    # Run TreeCollection
    tree_coll_tasks += "python map_names_to_species.py %s %s\n" % (output_tree_name, nex_alignments[sample_no] +
                                                                   "SpeciesMapping.txt")
    params = (script_dir, output_directory)
    tree_coll_tasks += time_code_beginning()
    tree_coll_tasks += "python %sexternal_tools/run_tree_collection.py %s\n" % params
    final_tree_name = output_directory + "tree.txt"  # remember to update the filename in run_tree_collection.py too
    tree_coll_tasks += time_code_end(times_folder + "TreeCollection_time_S%i.txt" % sample_no)
    tree_coll_tasks += 'echo "%s" >> %s\n' % (final_tree_name, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
    tree_coll_tasks += 'echo "Generated tree (TreeCollection)." >> S%i_progress.txt\n' % sample_no
    tree_coll_tasks += to_nexus_format_task(final_tree_name)
    tree_coll_tasks += monophyly_task(final_tree_name + ".tre", output_directory)
    tree_coll_tasks += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree_name, sample_no)

    return tree_coll_tasks


def write_pf_raxml_task(dict_of_alignments, sample_no, no_threads):
    substitution_model = "PROTGAMMAIWAG"
    task = ""
    input_folder = alignments_dir + "Sample_%i/" % sample_no
    partitions_path = input_folder + "sample_partitions.txt"  # name of partitions to be used by RaxML
    random_seed = choice(range(10000, 20000))
    pf_scheme_path = input_folder + "analysis/best_scheme.txt"
    task += "python make_partition_file.py %s %s\n" % (partitions_path, pf_scheme_path)
    output_folder = input_folder + "raxml_results/with_partitioning/PF/"
    make_dir(output_folder)
    task += time_code_beginning()
    params = (no_threads, substitution_model, random_seed, partitions_path, dict_of_alignments[sample_no], output_folder, "Tree")
    task += "raxmlHPC-PTHREADS -T %i -M -m %s -p %i -q %s -s %s -w %s -n %s\n" % params
    task += time_code_end(times_folder + "Raxml_W_PF_time_S%i.txt" % sample_no)
    final_tree = output_folder + "RAxML_bestTree.Tree"
    task += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
    task += 'echo "Generated tree (RaxML with PF partitions)." >> S%i_progress.txt\n' % sample_no
    task += "python map_names_to_species.py %s %s\n" % (final_tree, dict_of_alignments[sample_no].replace(".phy",
                                                                                                          ".nex")
                                                        + "SpeciesMapping.txt")
    task += to_nexus_format_task(final_tree)
    task += monophyly_task(final_tree + ".tre", output_folder)
    task += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree, sample_no)

    return task


def vanilla_raxml_job(dict_of_alignments, no_threads):
    substitution_model = "PROTGAMMAIWAG"
    for sample_no in dict_of_alignments.keys():
        filename = script_dir + "job_files/raxmlVanilla_S%i.sh" % sample_no
        output_folder = alignments_dir + "Sample_%i/raxml_results/without_partitioning/" % sample_no
        make_dir(output_folder)
        random_seed = choice(range(10000, 20000))
        params = (no_threads, substitution_model, random_seed, dict_of_alignments[sample_no], output_folder, "Tree")
        tasks = time_code_beginning()
        tasks += "raxmlHPC-PTHREADS -T %i -m %s -p %i -s %s -w %s -n %s\n" % params
        tasks += time_code_end(times_folder + "VanillaRaxML_time_S%i.txt" % sample_no)
        tasks += 'echo "Completed Vanilla-RaxML tree." >> S%i_progress.txt\n' % sample_no
        final_tree = output_folder + "RAxML_bestTree.Tree"
        tasks += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
        tasks += 'echo "Generated tree (Vanilla RaxML)." >> S%i_progress.txt\n' % sample_no
        tasks += "python map_names_to_species.py %s %s\n" % (final_tree, dict_of_alignments[sample_no].replace(".phy",
                                                                                                               ".nex")
                                                             + "SpeciesMapping.txt")
        tasks += to_nexus_format_task(final_tree)
        tasks += monophyly_task(final_tree + ".tre", output_folder)
        tasks += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree, sample_no)
        with open(r'%s' % filename, 'w') as raxmlJob:
            raxml_resources = res_text_sge("2G", 400, no_threads)
            raxmlJob.write(raxml_resources + tasks)

        run_job(filename, False)


def write_pf_task(sample_no, folder, search_option):
    if search_option == "greedy":
        job = 'echo "Starting partitioning (PartitionFinder, greedy)." >> S%i_progress.txt\n' % sample_no
        job += "python %sexternal_tools/PartitionFinderProtein.py %s\n" % (script_dir, folder)
    elif search_option == "rcluster":
        job = 'echo "Starting partitioning (PartitionFinder, relaxed-clustering)." >> S%i_progress.txt\n' % sample_no
        job += "python %sexternal_tools/PartitionFinderProtein.py %s --raxml --rcluster-percent 20\n" % (script_dir, folder)

    job += 'echo "Partitioning done (PartitionFinder)." >> S%i_progress.txt\n' % sample_no

    return job


def write_mtpan_task(nexus_files, sample_no):
    job = ""
    mtpan_script_dir = script_dir + "external_tools/MtPAN/Gap_script.R"
    job += 'echo "Starting partitioning (MtPAN)." >> S%i_progress.txt\n' % sample_no
    job += "Rscript %s %s long\n" % (mtpan_script_dir, nexus_files[sample_no])
    job += 'echo "Partitioning done (MtPAN)." >> S%i_progress.txt\n' % sample_no

    return job


def mtpan_raxml(sample_no, dict_of_alignments, min_cluster_size, max_cluster_size, no_threads):
    # perform RaxML analysis on partition schemes created by MtPAN. These partition schemes have a specified
    # range of sizes. Remember to set lines 57 and 58 to the min and max cluster sizes, respectively, in the
    # Gap_script.R script
    jobs = ""
    substitution_model = "PROTGAMMAIWAG"  # MtPAN 2013 matrix

    for k in range(min_cluster_size, max_cluster_size + 1):
        partition = alignments_dir + "Sample_%i/mtpan_partitions_%i.txt" % (sample_no, k)  # RaxML-style partition file.
        output_folder = alignments_dir + "Sample_%i/raxml_results/with_partitioning/MtPAN/size_%i/" % (sample_no, k)
        make_dir(output_folder)
        random_seed = choice(range(10000, 20000))
        params = (no_threads, substitution_model, random_seed, partition, dict_of_alignments[sample_no], output_folder,
                  "Tree")
        jobs += time_code_beginning()
        jobs += "raxmlHPC-PTHREADS -T %i -M -m %s -p %i -q %s -s %s -w %s -n %s\n" % params
        jobs += time_code_end(times_folder + "Raxml_W_MtPAN_time_S%i.txt" % sample_no)
        final_tree = output_folder + "RAxML_bestTree.Tree"
        jobs += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
        jobs += 'echo "Generated tree (RaxML with MtPAN partitions)." >> S%i_progress.txt\n' % sample_no
        jobs += "python map_names_to_species.py %s %s\n" % (final_tree, dict_of_alignments[sample_no].replace(".phy",
                                                                                                              ".nex")
                                                            + "SpeciesMapping.txt")
        jobs += to_nexus_format_task(final_tree)
        jobs += monophyly_task(final_tree + ".tre", output_folder)
        jobs += 'echo "Completed monophyly test for %s." >> S%i_progress.txt\n' % (final_tree, sample_no)

    return jobs


def phylobayes_job(dict_of_alignments, no_threads):

    for sample_no in dict_of_alignments.keys():
        filename = script_dir + "job_files/phylo_job_S%i.sh" % sample_no
        jobs = ""
        no_chains = 2
        burnin = 100
        max_discrepancy = 0.3
        min_sizes = 50
        site_variation_model = "cgam"
        aa_substitution_model = "wag"
        output_folder = alignments_dir + "Sample_%i/phylobayes_results/" % sample_no
        make_dir(output_folder)
        jobs += time_code_beginning()
        chain_address = output_folder + "chain"
        params = (dict_of_alignments[sample_no], no_chains, burnin, max_discrepancy, min_sizes, site_variation_model, aa_substitution_model, chain_address)
        jobs += 'echo "Generating tree (PhyloBayes)." >> S%i_progress.txt\n' % sample_no
        jobs += "pb -d %s -nchain %i %i %i %i -%s -%s %s\n" % params
        take_tree_every = 10  # trees
        jobs += "readpb -x %i %i %s\n" % (burnin, take_tree_every, chain_address)
        jobs += time_code_end(times_folder + "PhyloBayes_time_S%i.txt" % sample_no)
        jobs += 'echo "Generated tree (PhyloBayes)." >> S%i_progress.txt\n' % sample_no
        final_tree = chain_address + ".con.tre"
        jobs += "python map_names_to_species.py %s %s\n" % (final_tree, dict_of_alignments[sample_no].replace(".phy",
                                                                                                              ".nex")
                                                            + "SpeciesMapping.txt")
        jobs += 'echo "%s" >> %s\n' % (final_tree, alignments_dir + "Sample_%i/tree_list.txt" % sample_no)
        jobs += to_nexus_format_task(final_tree)
        jobs += monophyly_task(final_tree + ".tre", output_folder)

        with open(r'%s' % filename, 'w') as phylo_job:
            phylo_resources = res_text_sge("2G", 500, no_threads)
            phylo_job.write(phylo_resources + jobs)

        run_job(filename, False)


def write_conversion_tasks(nexus_files):
    jobs = ""
    for sample_no in nexus_files.keys():
        jobs += "python ./external_tools/nexus_to_phy_sequential.py %s\n" % nexus_files[sample_no]
        jobs += 'echo "Converted alignments from NEXUS to PHYLIP format." >> S%i_progress.txt\n' % sample_no

    return jobs


def to_nexus_format_task(tree):
    task = ""
    task += "python %s %s\n" % (script_dir + "to_nexus.py", tree)

    return task


def res_text_sge(memory, no_hours, no_threads):
    text = "#$ -l h_vmem=%s\n" % memory +\
           "#$ -l tmem=%s\n" % memory +\
           "#$ -S /bin/bash\n" + \
           "#$ -j y\n" +\
           "#$ -l h_rt=%i:00:0\n" % no_hours +\
           "#$ -pe smp %i\n" % no_threads +\
           "#$ -R y\n" +\
           "#$ -cwd\n" +\
           "\n"

    return text


def res_text_legion(memory, no_hours, no_threads):
    text = "#!/bin/bash -l\n" +\
           "#$ -l mem=%s\n" % memory +\
           "#$ -S /bin/bash\n" + \
           "#$ -j y\n" +\
           "#$ -l h_rt=%i:00:0\n" % no_hours +\
           "#$ -l thr=%i\n" % no_threads +\
           "#$ -wd /home/zcabaol/Scratch/\n" +\
           "cd $TMPDIR" +\
           "\n"

    return text


def make_dir(directory):
    if not path.exists(directory):
        mkdir(directory)


def run_job(filename, wait):
    args = split("qsub " + filename)
    p = Popen(args).wait()
    if wait:  # pipeline will wait for the job that does conversion to finish.
        s = Popen("qstat", stdin=PIPE, stdout=PIPE)
        status_msg = s.communicate("qstat")[0]
        while status_msg:
            s = Popen("qstat", stdin=PIPE, stdout=PIPE)
            status_msg = s.communicate("qstat")[0]

# methods for measuring time of program execution


def time_code_beginning():
    return '''_utime="$( TIMEFORMAT='%lU';time -p (\n'''


def time_code_end(fname):
    return ''') 2>&1 1>/dev/null )"\necho "$_utime" >> %s\n''' % fname


