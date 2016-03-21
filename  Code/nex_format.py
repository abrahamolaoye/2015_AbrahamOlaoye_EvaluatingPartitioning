__author__ = 'Lovecraft'
from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from glob import *
import sys

dir = sys.argv[1]

for new_file in [f for f in glob(dir + ".fas")]:
    alignment = AlignIO.read(open(new_file), "fasta", alphabet=Gapped(IUPAC.protein))
    nexus_file = open(new_file+".nex", "w")
    nexus_file.write(alignment.format("nexus"))
