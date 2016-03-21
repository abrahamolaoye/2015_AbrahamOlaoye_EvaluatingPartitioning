__author__ = 'Abraham'
from constants import *
from glob import *
# goes through a range of sample folders and, for each one, makes a list of OGs with their alignments in that folder.
# this is required by the script that makes distance-variance matrices, maps and labels for use by tree-collection
for k in range(1, no_of_samples + 1):
    orthologous_groups = [f for f in glob(alignments_dir + "Sample_%i/*.fas" % k)]
    filename = alignments_dir + "Sample_%i/tree_collection_results/OG_list.txt" % k
    with open(r'%s' % filename, 'w') as newFile:
        for aligned_og in orthologous_groups:
            newFile.write(aligned_og + '\n')

