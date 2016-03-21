__author__ = 'Lovecraft'
import sys
from constants import *

alignment = sys.argv[1]
output_fname = sys.argv[2]
# defining blocks of sites for use by RaxML
# each OG is a partition.
# this method works for alignments in the Nexus format.
with open(r'%s' % alignment, 'r') as alignment_file:
        line = ""
        while line != "begin sets;\n":
            line = alignment_file.readline()
        partitions = ""
        for j in range(1, sample_size + 1):
            line = alignment_file.readline()
            partitions += "WAG, " + line.replace("charset ", "").replace(".fa.out.nex", "").replace(".", "")\
                .replace(";", '')

with open(r'%s' % output_fname, 'w') as newFile:
    newFile.write(partitions)


