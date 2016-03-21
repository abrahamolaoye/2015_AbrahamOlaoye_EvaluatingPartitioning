__author__ = 'Lovecraft'
# given a list of text files, each containing the monophyly results for their respective trees
# create a text file containing the tree names and their scores, listed in descending order
import sys
from constants import *

filename = sys.argv[1]  # text file containing a list of the names of the monophyly results
output_filename = sys.argv[2]  # output name of the rankings text file
clade_path = script_dir + "clade_list.txt"
list_of_tree_score_pairs = []


def get_score(pair):
    return pair[1]

with open(r'%s' % filename, 'r') as list_file:
    line = list_file.readline()

    while line:
        res_file = line.replace('\n', '')
        with open(r'%s' % res_file, 'r') as results:
            line = results.readline()
            line = results.readline()
            line = results.readline()  # monophyly results are on line 3
            split_line = line.replace('-', '0').split()
            tree_name = split_line[0]
            split_line = map(int, split_line[1:])
            numerator = sum(split_line[::2])
            denominator = sum(split_line[1::2])
            if denominator != 0:
                overall_monophly_perc = (float(numerator) / float(denominator)) * 100
            list_of_tree_score_pairs.append((tree_name, overall_monophly_perc))
        line = list_file.readline()

with open(r'%s' % output_filename, 'w') as monophyly_scores:
    list_of_tree_score_pairs.sort(key=get_score, reverse=True)
    for tree_and_score in list_of_tree_score_pairs:
        monophyly_scores.write("%s - %i%%\n" % tree_and_score)




