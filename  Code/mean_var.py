from __future__ import division
__author__ = 'Lovecraft'
# script that gets the mean and variance of the times taken for some method of tree_inference.

import sys
from glob import *

directory = sys.argv[1]
codename = sys.argv[2]  # codenames are the shortened names given to a method of tree inference
# Raxml with partition finder - Raxml_W_PF
# Tree collection - TreeCollection
# Raxml analysis using by-gene partitions - partitioned_by_gene_raxml
# Phylobayes - PhyloBayes
# Raxml without any partitioning - VanillaRaxML

time_files = [f for f in glob(directory + codename + "_time_S*")]
times = []
for txt_file in time_files:
    with open(r'%s' % txt_file, 'r') as time_file:
        line = time_file.readline()
        real_time = float((line.split()[1])) / float(3600)  # take number on first line, record it in terms of hours rather than
        #  seconds
        times.append(real_time)

mean = sum(times) / len(times)
variance = float(sum([float(pow(n - mean, 2)) / float(len(times) - 1) for n in times]))

print("Mean: %f, variance: %f" % (mean, variance))





