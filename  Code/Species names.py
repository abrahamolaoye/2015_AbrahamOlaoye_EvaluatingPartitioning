__author__ = 'Abraham'
# Returns the names of all species in the group of OGs
from glob import *

list_of_ogs = [f for f in glob("/Users/Lovecraft/Documents/Internship project/Dataset/OGs/*.fa.out")]
species_names = []

for og_file in list_of_ogs:
    with open(r'%s' % og_file, 'r') as og:
        line = og.readline()
        while line:
            if line[0] == '>':
                species = line[1:]
                if species not in species_names:
                    species_names.append(species)

            line = og.readline()

species_names.sort()
with open(r'Species_list.txt', 'w') as species_file:
    species_file.write(str(len(species_names)) + '\n')
    for name in species_names:
        species_file.write(name)

print(len(species_names))
