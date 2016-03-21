#!/usr/bin/python

import pyopa
import os
import threading
import sys

file_list = sys.argv[1]
output = sys.argv[2]


#threading.stack_size(1310720000000000000)
#apparently this is needed for pyopa to align/calc distances
threading.stack_size(1310720000)

###########################################################################################
#GENERATE ENVIRONMENTS
#taken from pyopa example

defaults = pyopa.load_default_environments()
env_list = defaults['environments']
log_pam1_env = defaults['log_pam1']
#generates a signle AlignmentEnvironment
# with a pam distance of 250
generated_env = pyopa.generate_env(log_pam1_env, 250)

#generates 1000 environments for different pam distances
gen_env_list = pyopa.generate_all_env(log_pam1_env, 1000)


#find pairwise distance and variance between two sequences
#needed to calc dist-var between two aligned sequences
dms = pyopa.MutipleAlEnv(gen_env_list, log_pam1_env)

###########################################################################################
#INITIALISE MATRIX TO 0's AND GENMAP to -1

#initialise matrix
def init_matrix(species_gene):
	matrix = [x[:] for x in [[0]*len(species_gene)]*len(species_gene)]
	return matrix


#initialise genmap
def init_genmap(label):
	map = [x[:] for x in [[-1]*len(label)]]
	return map
###########################################################################################
#read sequence file and put species name and sequence in dictionary
#this is done for each orthologous group

def read_sequence(filename):
	with open(filename) as f:
		data = f.readlines()
	list=[]
	pairs=dict()
	for i in range(len(data)):
		list.append(data[i].split('\n'))
		list[i]=list[i][:-1]
	for i in range(len(list)):	
		if list[i][0].startswith(">"):
			species_name=list[i][0][1:]
			pairs[species_name]=list[i+1][0]
		else:
			continue
	return pairs
###########################################################################################
#align two sequences and find the dist-var

#from pyopa: to align sequence
def nt_align(s1, s2, env, is_global, aligned_strs):
    tmp_aligned_strings = pyopa.align_strings(s1, s2, env, is_global)
    aligned_strs.extend(tmp_aligned_strings)

#align a pair of sequences
def align_sequence(s1,s2):
	aligned_strings = []
	t = threading.Thread(None, nt_align,
                     'Aligning Thread', (s1, s2, generated_env, False, aligned_strings))
	t.start()
	t.join()
	return aligned_strings


#find pam_distance and variance
#use for main function calc_distvar
def find_dist_var(alignment_pair):
	tmp = dms.estimate_pam(alignment_pair[0],alignment_pair[1])
	list=[]
	list.append(tmp[1])
	list.append(tmp[2])
	return list



#CALCULATE DISTVAR

def calc_distvar(sequence1,sequence2):
	sp1=pyopa.Sequence(sequence1)
	sp2=pyopa.Sequence(sequence2)		
	aligned_pair = align_sequence(sp1,sp2)
	distvar = find_dist_var(aligned_pair)
	return distvar


################################################################################################
#first need to check through every orthologous group and record every species
#then sort them alphabetically
#gives us "labels" output and ordering for "genmap"
#USE THESE FUNCTIONS FOR FUNCTION IN NEXT SECTION
#TO CREATE FILES FOR GENMAP
 
def find_species(filename,species_dict):
	with open(filename) as f:
		data = f.readlines()
	list=[]
	for i in range(len(data)):
		list.append(data[i].split('\n'))
		list[i]=list[i][:-1]
	for i in range(len(list)):	
		if list[i][0].startswith(">"):
			species_name=list[i][0][1:]
			species_dict[species_name]=0
	return species_dict


#list of filenames from file list
def read_filelist(file_list):
	with open(file_list) as f:
		content = f.readlines()
	for i in range(len(content)):
		content[i]=content[i][:-1]
	return content


def get_filename(file_list,i):
	return file_list[i]

##############################################
#THIS CREATES ORDERED LABELS FOR SPECIES.... NEEDED FOR GENMAP

def create_complete_species_list(file_list):
	labels_list=[]
	files = read_filelist(file_list)
	species_global=dict()
	for i in range(len(files)):
		labels = find_species(files[i],species_global)
	for name in species_global:
		labels_list.append(name)
	labels=sorted(labels_list)
	return labels

################################################################################################
#sort dictionary_gene list

def sort_gene_species(species_dict):
	species_gene=[]
	for name in species_dict:
		species_gene.append(name)
	species_gene = sorted(species_gene)
	return species_gene
################################################################################################
#find the genmap index... place species into correct position in genmap

def find_genmap_index(species_gene,labels):
	for i in range(len(labels)):
		if species_gene==labels[i]:
			return i
################################################################################################
#retrieve dist and var from pyopa output
def get_dist(dist_var):
	return dist_var[0]

def get_var(dist_var):
	return dist_var[1]

########################################################################################################
#creates dist var matrix and genmap

def create_matrix(sorted_species_gene,gene_sequence_index,genmap,matrix,global_species_dict,k):
	for i in range(len(sorted_species_gene)):
		tmp_index = find_genmap_index(sorted_species_gene[i],global_species_dict)
		genmap[k][0][tmp_index] = i+1 
		for j in range(i+1,len(sorted_species_gene),1):
			sp1_sequence = gene_sequence_index[sorted_species_gene[i]]
			sp2_sequence = gene_sequence_index[sorted_species_gene[j]]
			dist_var = calc_distvar(sp1_sequence,sp2_sequence)
			dist = get_dist(dist_var)
			var = get_var(dist_var)
			matrix[k][i][j] = dist
			matrix[k][j][i] = var

##################################################################






#reads filenames from filelist
#creates labels (ordered list of species)
#creates empty dist-var matrix
#creates empty map

#for each orthologous group:
#creates sorted index (so they are in same order as labels) of species:sequences (eg. sp1: seq1, sp2:seq2 etc.)
#add empty matrix of size (no.species in group x no.species in group)
#add empty map of size (no. species in whole set)
#fills matrix and map 
#for matrix: distances in upper triangle, variances in lower triangle.
#for map: -1 if species isn't present in group. integer value for position species is in dist-var matrix.
#eg. 4 in 8th column of map means that Sp8 is the species represented in 4th column and row of matrix


filenames = read_filelist(file_list)
treecollection_species = create_complete_species_list(file_list)
treecollection_matrix=[]
treecollection_genmap=[]
for k in range(len(filenames)):
	filename = get_filename(filenames,k)
	gene_sequence_index = read_sequence(filename)
	sorted_gene_list = sort_gene_species(gene_sequence_index)
	treecollection_matrix.append(init_matrix(sorted_gene_list))
	treecollection_genmap.append(init_genmap(treecollection_species))
	create_matrix(sorted_gene_list,gene_sequence_index,treecollection_genmap,treecollection_matrix,treecollection_species,k)
	




##################################################################
#write to files..

#write genmap

with open(output+"Map" +".txt", "w") as my_file:
	my_file.write(str(len(filenames)) + " " + str(len(treecollection_species)))
	for i in range(len(treecollection_genmap)):
		my_file.write("\n")
		for j in range(len(treecollection_genmap[i][0])):
			my_file.write(str(treecollection_genmap[i][0][j]) + " ")

#write labels

with open(output+"Labels" +".txt", "w") as my_file:
	my_file.write(str(len(treecollection_species)) + "\n")
	for i in range(len(treecollection_species)):
		my_file.write(treecollection_species[i] + " ")

#write matrices

with open(output+"DistVar" +".txt", "w") as my_file:
	my_file.write(str(len(treecollection_matrix)))
	for i in range(len(treecollection_matrix)):
		my_file.write("\n")
		my_file.write(str(len(treecollection_matrix[i])) + " " + str(len(treecollection_matrix[i])) + " " +  str(i+1) + "\n")
		for j in range(len(treecollection_matrix[i])):
			for k in range(len(treecollection_matrix[i][j])):
				if k==len(treecollection_matrix[i][j])-1:
					my_file.write(str(treecollection_matrix[i][j][k]) + "\n")
				else:
					my_file.write(str(treecollection_matrix[i][j][k]) + " ")






