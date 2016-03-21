__author__ = 'Lovecraft'
# visits samples folders and, for each one, produces a text file containing a list of the trees made.
from constants import *

for j in range(1, no_of_samples + 1):
    filename = "tree_list.txt"
    with open(r'%s' % alignments_dir + "Sample_%i/" + filename, 'w') as list_file:


