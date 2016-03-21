#!/usr/bin/env python

import tree_collection
import sys

folder_name = sys.argv[1]

mat = open(folder_name + 'DistVar.txt').read()
map_ = open(folder_name + 'Map.txt').read()
lab = open(folder_name + "Labels.txt").read()
tre = open(folder_name + 'initial_tree.txt').read()

result = tree_collection.compute(mat, map_, lab, tre, 8, False)
with open(r'%s' % folder_name + "tree.txt", 'w') as tree_coll_tree:
    tree_coll_tree.write(result[0])


