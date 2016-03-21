__author__ = 'Lovecraft'
# break up a set of unaligned OGs into separate groups

from glob import *


def make_tasks(group_no):
    tasks = ""
    for og_path in list_of_og_groups[group_no - 1]:
        tasks += "prank -d=$%s -F -uselogs -o=$%s.prank" % (og_path, og_path)

    return tasks

all_ogs = [og for og in glob("./OGs/*.fa.out")]
aligned_ogs = [og.replace(".prank.best.fas", "") for og in glob("./OGs/*.fas")]

for og in aligned_ogs:  # skip OGs which have been aligned
    index = all_ogs.index(og)
    del all_ogs[index]

group_sizes = 150
list_of_og_groups = []
if len(all_ogs) / group_sizes == 0:
    no_jobs = len(all_ogs) / group_sizes
else:
    no_jobs = (len(all_ogs) / group_sizes) + 1

start_index = 0
end_index = group_sizes

for i in range(1, no_jobs + 1):
    if (len(all_ogs) - start_index) > group_sizes - 1:
        list_of_og_groups.append(all_ogs[start_index: end_index])
        start_index += group_sizes
        end_index += group_sizes
    else:
        list_of_og_groups.append(all_ogs[start_index:])


def make_tasks(group_no):
    tasks = ""
    print( "%i has %i groups" % (group_no, len(list_of_og_groups[group_no - 1])))
    for og_file in list_of_og_groups[group_no - 1]:
        tasks += "prank -d=%s -F -uselogs -o=%s.prank\n" % (og_file, og_file)

    return tasks


for k in range(1, no_jobs + 1):
    filename = "Run_prank_%i.sh" % k
    with open(r'%s' % filename, 'w') as prank_job:
        prank_job.write("#$ -l h_vmem=1G\n"
                        "#$ -l tmem=1G\n"
                        "#$ -S /bin/bash\n"
                        "#$ -j y\n"
                        "#$ -l h_rt=50:00:0\n"
                        "#$ -R y\n"
                        "#$ -cwd\n\n"

                        "%s" % make_tasks(k))


