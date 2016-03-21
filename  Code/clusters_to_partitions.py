__author__ = 'Abraham'
from constants import *
from glob import *
import sys

sample_no = int(sys.argv[1])


def clusters_to_partitions():
    clusters_in_sample = [f for f in glob(alignments_dir + 'Sample_%i/sample_alignment.nex.kclusters.*' % sample_no)]
    for cluster in clusters_in_sample:
        line = ""
        cluster_no = int(cluster[len(cluster) - 1])
        with open(r'%s' % cluster, 'r') as newCluster:
            line = newCluster.readline()
            line = newCluster.readline()  # read second line, number of numbers = number of blocks.
            no_of_blocks = len(line.split())
            site_posns = []
            for n1 in range(1, no_of_blocks + 1):
                site_posns.append([])  # each new list will hold the site positions for those in a given cluster.

            list_of_chars = []
            while line:
                line = newCluster.readline()
                list_of_chars.extend(line.split())
            posn = 1
            for char in list_of_chars:
                if char is not '':
                    site_posns[int(char) - 1].append(posn)
                    posn += 1
            write_file(site_posns, cluster_no)


def write_file(site_posns, file_no):
    block_no = 1
    filename = alignments_dir + "Sample_%i/mtpan_partitions_%i.txt" % (sample_no, file_no)
    with open(r'%s' % filename, 'w') as PartitionFile:
        blocks = ""
        for posns_for_cluster in site_posns:
            new_block = ""
            for posn in posns_for_cluster:
                new_block += str(posn) + ","
            new_block = new_block.rstrip(",")
            blocks += "WAG, p%i =%s\n" % (block_no, new_block)
            block_no += 1
        PartitionFile.write(blocks)

clusters_to_partitions()
