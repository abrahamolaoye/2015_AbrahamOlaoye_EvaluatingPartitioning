__author__ = 'Abraham'
from constants import *
import sys

# uses a PartitionFinder partition scheme to make a RaxML-style partition file.
partitions_path = sys.argv[1]
scheme_path = sys.argv[2]


def make_partition_file():
    # First, need to create the partition file for RaxML
    best_scheme_path = scheme_path
    text = ""
    with open(r'%s' % best_scheme_path, 'r') as best_scheme:
        line = best_scheme.readline()
        while line != "RaxML-style partition definitions\n":
            line = best_scheme.readline()
        while line:
            line = best_scheme.readline()
            text += line
    with open(r'%s' % partitions_path, 'w') as partition_file:
        partition_file.write(text)

make_partition_file()
