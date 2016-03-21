__author__ = 'Lovecraft'
#!/usr/bin/python

import sys
filename = sys.argv[1]

from Bio import Phylo

Phylo.convert(filename,'newick',filename+".tre",'nexus')
