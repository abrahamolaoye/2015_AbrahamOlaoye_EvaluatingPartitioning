__author__ = 'Abraham'
# roots a tree and eliminates negatives

import sys
from ete2 import *
from constants import *

filename = sys.argv[1]
tree = Tree(filename)
midpoint_node = tree.get_midpoint_outgroup()
rooted_tree = tree.set_outgroup(midpoint_node)

with open(r'%s' % filename, 'w') as new_tree:
    new_tree.write(tree.write(format=5).replace(":-", ":"))  # get rid of negative signs too



