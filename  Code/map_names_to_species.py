__author__ = 'Abraham'
# Takes a tree and replaces it with one containing the actual species names rather than the mapped names.

import sys

filename = sys.argv[1]
mapping = sys.argv[2]

with open(r'%s' % filename, 'r') as treeFile:
    tree_txt = treeFile.read()


with open(r'%s' % mapping, 'r') as Mapping:
    line = Mapping.readline()
    new_tree_txt = tree_txt
    while line:
        name_and_species = line.split()
        new_tree_txt = new_tree_txt.replace("%s:" % name_and_species[0], "%s:" % name_and_species[1])
        line = Mapping.readline()

with open(r"%s" % filename, 'w') as newTreeFile:
    newTreeFile.write(new_tree_txt)







