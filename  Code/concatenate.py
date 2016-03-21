__author__ = 'Abraham'
# converts a set of fasta files to nexus format, concatenates each and combines them all into one large alignment.
from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from glob import *
from constants import *
import sys

sample_no = int(sys.argv[1])


def concat_alignment(files, output):
    nexi = [(fname.replace(alignments_dir, '').replace(".", "").replace("-", ""), Nexus.Nexus(fname)) for fname in
            files]
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open(output, 'w'))


nex_files = [f.replace("Sample_%i/" % sample_no, "") + ".nex" for f in glob(alignments_dir + "Sample_%i/*.fas" %
                                                                            sample_no)]
sample_alignment_path = alignments_dir + 'Sample_%i/sample_alignment.nex' % sample_no
concat_alignment(nex_files, sample_alignment_path)  # create concatenated alignment

