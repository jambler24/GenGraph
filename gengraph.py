# Dependencies

import networkx as nx


# Built in

from subprocess import call

import argparse

import sys

import csv

import copy

import time

import pickle


path_to_muscle = 'muscle3.8.31_i86darwin64'

path_to_clustal = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/clustalo'

path_to_mafft = 'mafft'

path_to_kalign = '/Users/panix/Dropbox/Programs/tools/kalign2/kalign'

path_to_progressiveMauve = '/Applications/Mauve.app/Contents/MacOS/progressiveMauve'


global_aligner = 'progressiveMauve'

local_aligner = 'mafft'

global_thread_use = '16'

'''

Conventions:

Orientation is relative to the sequence in the node, which is default the reference organism sequence.

- For normal orientation 				Left(+) < Right(+)		Left 10  Right 20		ATGC	(Reference)
- For reverse-compliment orientation 	Left(-) < Right(-)		Left -10 Right -20		GCAT

'''
# ---------------------------------------------------- New functions under testing

def get_neighbours_context(graph, source_node, label, dir='out'):
	'''Retrieves neighbours of a node using only edges containing a certain label'''
	'''Created by shandu.mulaudzi@gmail.com'''
	list_of_neighbours = []
	in_edges = graph.in_edges(source_node, data='sequence')
	out_edges = graph.out_edges(source_node, data='sequence')
	for (a,b,c) in in_edges:
		if label in c['sequence'].strip(','):
			list_of_neighbours.append(a)
	for (a,b,c) in out_edges:
		if label in c.strip(','):
			list_of_neighbours.append(b)
		
	list_of_neighbours = list(set(list_of_neighbours)) # to remove any duplicates (like bi_correct_node)	
		
	return list_of_neighbours



# ---------------------------------------------------- General functions 

def input_parser(file_path, parse_as='default'):
	if file_path[-3:] == ".fa" or file_path[-6:] == ".fasta":
		input_file = open(file_path, "r")
		output_list = []
		# set variables
		sequence_details = ""
		sequence = ""

		for line in input_file:
			if line[0] == ">":
				if len(sequence_details) > 0:
					sequence_details = sequence_details.strip()
					sequence = sequence.strip()
					sequence = sequence.replace("\n","")
					gene_ID_dict = {"gene_details" : sequence_details[1:], "DNA_seq" : sequence}
					output_list.append(gene_ID_dict)
					sequence_details = ""
					sequence = ""
				sequence_details = line
			else:
				sequence += line
		
		sequence_details = sequence_details.strip()
		
		sequence = sequence.strip()
		sequence = sequence.replace("\n","")
		sequence = sequence.replace(" ","")
		
		gene_ID_dict = {"gene_details" : sequence_details[1:], "DNA_seq" : sequence}

		output_list.append(gene_ID_dict)

		return output_list

	if file_path[-4:] == ".bim":
		input_file = open(file_path, "r")
		return input_file
		
	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result=numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".ssv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=';')
		data_matrix = list(data_table)
		result=numpy.array(data_matrix)
		return result
		
	if file_path[-4:] == ".txt":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts
		
	if file_path[-4:] == ".vcf":
		list_of_dicts  = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if line[0:2] != "##" and line[0] == "#":
				vcf_headder_line = line.split('\t')
				vcf_headder_line[0] = vcf_headder_line[0][1:]
				vcf_headder_line[-1] = vcf_headder_line[-1].strip()

			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {vcf_headder_line[0]:entries[0], vcf_headder_line[1]:entries[1], vcf_headder_line[2]:entries[2], vcf_headder_line[3]:entries[3], vcf_headder_line[4]:entries[4], vcf_headder_line[5]:entries[5], vcf_headder_line[6]:entries[6], vcf_headder_line[7]:entries[7], vcf_headder_line[8]:entries[8], vcf_headder_line[9]:entries[9].strip(), 'ORIGIN':entry_label} 

				list_of_dicts.append(entry_dict)
		return list_of_dicts
	
	if file_path[-5:] == ".diff":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts

	if file_path[-5:] == "kbone":
		backb_listOlists = []
		in_file = open(file_path, 'r')
		for line in in_file:
			curr_list = []
			line = line.strip()
			curr_list = line.split('\t')
			backb_listOlists.append(curr_list)
		return backb_listOlists

	if file_path[-4:] == ".gtf":
		list_of_lists = []
		in_file = open(file_path, 'r')
		for line in in_file:
			entries = line.split('\t')
			entries[8] = entries[8].strip('\n')
			#entries_extra_info = entries[8].split(';')

			if '; ' in entries[8]:
				entries_extra_info = entries[8].split('; ')
			else:
				entries_extra_info = entries[8].split(';')

			extra_info_dict = {}

			for info_byte in entries_extra_info:

				if info_byte != '#' and len(info_byte) > 1:

					info_byte = info_byte.split(' ')

					extra_info_dict[info_byte[0]] = info_byte[1]

					entries[8] = extra_info_dict

			list_of_lists.append(entries)


		return list_of_lists

	if file_path[-5:] == ".gff3" or file_path[-4:] == ".gff":
		list_of_lists  = []
		in_file = open(file_path, 'r')
		entry_label = file_path
		
		if parse_as == 'default':
			for line in in_file:
				if not line.startswith('#'):
					entries = line.split('\t')
					if len(entries) > 5:
						entries[8] = entries[8].strip('\n')

						entries_extra_info = entries[8].split(';')

						extra_info_dict = {}

						for info_byte in entries_extra_info:
							info_byte = info_byte.split('=')
							extra_info_dict[info_byte[0]] = info_byte[1]

						entries[8] = extra_info_dict
				
						list_of_lists.append(entries)
		if parse_as == 'gtf':
			# Reformatting as gtf
			for line in in_file:
				if not line.startswith('#'):
					entries = line.split('\t')
					if len(entries) > 5:
						entries[8] = entries[8].strip('\n')

						entries_extra_info = entries[8].split(';')

						extra_info_dict = {}
						gtf_info_dict = {}

						for info_byte in entries_extra_info:
							info_byte = info_byte.split('=')
							extra_info_dict[info_byte[0]] = info_byte[1]

						if 'locus_tag' in extra_info_dict.keys():

							gtf_info_dict['gene_id'] = extra_info_dict['locus_tag']

							entries[8] = gtf_info_dict
					
							list_of_lists.append(entries)


		return list_of_lists

	if file_path[-4:] == ".gff OLD":
		list_of_dicts  = []
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if not line.startswith('#'):
				entries = line.split('\t')
				entries[8] = entries[8].strip('\n')
				entries_extra_info = entries[8].split(';')

				if entries[2] == 'gene':
					
					NOTE = ''
					for add_info in entries_extra_info:
						if 'locus_tag' in add_info:
							LOCUS = add_info[10:]
						if 'Name'in add_info:
							SYMBOL = add_info[5:]
						if 'Note'in add_info:
							NOTE = add_info[5:]



						#row['LOCUS'] row['SYMBOL'] row['SYNOYM'] row['LENGTH']  row['START'] row['STOP']  row['STRAND']  row['NAME']  row['CHROMOSOME']  row['GENOME ONTOLOGY']  row['ENZYME CODE']  row['KEGG']  row['PATHWAY']  row['REACTION'] row['COG'] row['PFAM']  row['OPERON'] 

					entry_dict = {'CHROM':entries[0], 'LOCUS':LOCUS, 'START':entries[3], 'STOP':entries[4], 'STRAND':entries[6], 'SYMBOL':SYMBOL, 'INFO':entries_extra_info, 'NOTE':NOTE}
					list_of_dicts.append(entry_dict)
		return list_of_dicts

def parse_seq_file(path_to_seq_file):

	seq_file_dict = input_parser(path_to_seq_file)

	A_seq_label_dict = {}
	A_input_path_dict = {}
	ordered_paths_list = []
	anno_path_dict = {}

	for a_seq_file in seq_file_dict:
		print a_seq_file
		A_seq_label_dict[a_seq_file['aln_name']] = a_seq_file['seq_name']
		A_input_path_dict[a_seq_file['seq_name']] = a_seq_file['seq_path']
		ordered_paths_list.append(a_seq_file['seq_path'])
		anno_path_dict[a_seq_file['seq_name']] = a_seq_file['annotation_path']

	return A_seq_label_dict, A_input_path_dict, ordered_paths_list, anno_path_dict

def compliment_DNA(seq_in):
	# Returns the string with the complimentry base pairs
	seq_out = ""
	for char in seq_in:
		seq_out += compliment_Base(char)
	return seq_out

def reverse(text):
	# returns a reversed string using recursion
	return text[::-1]

def reverse_compliment(seq_in):
	"Returns the reverse compliment of a sequence"
	comp_seq = compliment_DNA(seq_in)
	rev_comp_seq = reverse(comp_seq)
	return rev_comp_seq

def compliment_Base(base_in):
	# Returns the compliment of the given base. Assuning DNA
	# Non-DNA bases are returned as a empty string

	base_in = base_in.upper()
	if base_in == "A":
		return "T"
	if base_in == "T":
		return "A"
	if base_in == "C":
		return "G"
	if base_in == "G":
		return "C"
	if base_in == "N":
		return "N"
	else:
		return ""

def is_even(a_num):
    x = int(a_num)
    if x % 2 == 0:
        return True
    else:
        return False

def bp_distance(pos_A, pos_B):
	# Assuming that you cannot get a negitive and a positive base pair value set
	abs_pos_A = abs(int(pos_A))
	abs_pos_B = abs(int(pos_B))

	dist = abs(abs_pos_A - abs_pos_B)

	dist += 1

	return dist

def bp_length_node(node_name):
	''' Returns the length of a node in base pairs '''
	iso_string = node_name['present_in'].split(',')[0]
	left_pos = node_name[iso_string + '_leftend']
	right_pos = node_name[iso_string + '_rightend']

	return bp_distance(left_pos, right_pos)


# ---------------------------------------------------- Graph generating functions 


def link_nodes_2(graph_obj, sequence_name, node_prefix='gn'):

	print 'link nodes 2'
	
	# To make sure all new nodes have unique names, a count is kept

	if 'add_count' in graph_obj.graph.keys():
		#print 'found'
		add_count_val = graph_obj.graph['add_count']
	else:
		#print 'none found'
		graph_obj.graph['add_count'] = 0
		add_count_val = 0

	# This list is a list of all the existing nodes start and stop positions with the node name formatted at start, stop, node_name

	pos_lol = []

	for node,data in graph_obj.nodes_iter(data=True):
		#print data
		#print data[left_end_name], data[right_end_name], node
		#print 'gggggg'
		#print data

		left_end_name = sequence_name + '_leftend'
		right_end_name = sequence_name + '_rightend'		
		#print '------------'
		#print data
		#print node
		#print sequence_name
		if sequence_name in data['present_in'].split(','):
			#print 'right node'
			#print data

			new_left_pos = abs(int(data[left_end_name]))
			new_right_pos = abs(int(data[right_end_name]))

			if new_left_pos > new_right_pos:
				temp_left = new_right_pos
				temp_right = new_left_pos
				new_left_pos = temp_left
				new_right_pos = temp_right


			if int(new_left_pos) != 0 and int(new_right_pos) != 0 or node == sequence_name + '_start':
				pos_lol.append([new_left_pos, new_right_pos, node])

	count = 0

	gap_lol = []

	pos_lol = sorted(pos_lol)

	#for a_list in pos_lol:
	#	print a_list

	sorted_pos_lol = pos_lol

	#print sorted_pos_lol

	# -------------------------------------------------------- Here we determine if there are any gaps in the sequence 

	while count < (len(sorted_pos_lol) - 1):

		if sorted_pos_lol[count][1] != (sorted_pos_lol[count + 1][0] - 1) and sorted_pos_lol[count][1] != (sorted_pos_lol[count + 1][0]):
			gap_list = [(sorted_pos_lol[count][1] + 1), ((sorted_pos_lol[count + 1][0] - 1)), sequence_name]
			gap_lol.append(gap_list)

		count +=1


	for entry in gap_lol:
		if entry[0] > entry[1]:
			#print 'wtf'
			#print entry
			new_entry_0 = entry[1]
			new_entry_1 = entry[0]
			entry[1] = new_entry_1
			entry[0] = new_entry_0

	# -------------------------------------------------------- Add new nodes

	#print gap_lol

	new_node_list = []

	for new_seq_node in gap_lol:
		#print new_seq_node
		gap_len = new_seq_node[1] - new_seq_node[0] + 1
		node_dict = {new_seq_node[2] + '_leftend':str(new_seq_node[0]), new_seq_node[2] + '_rightend':str(new_seq_node[1]), 'present_in':sequence_name}

		#print node_dict
		node_name = node_prefix + '_' + str(add_count_val)

		new_seq_node.append(node_name)

		graph_obj.add_node(node_name, node_dict)


		new_node_list.append([new_seq_node[0], new_seq_node[1], node_name])

		add_count_val += 1





	#print pos_lol
	#print gap_lol
	#print new_node_list

	# -------------------------------------------------------- Combine the origional node list with the list of gaps generated earlier

	all_node_list = pos_lol + new_node_list
	#print '\n'

	all_node_list = sorted(all_node_list)


	#print 'This is the list*************************'
	#print all_node_list


	# -------------------------------------------------------- Here we add edges to the graph, weaving in the new nodes

	count = 0

	#graph_obj.add_edges_from([('H37Rv_seq_4', 'Aln_region_7', dict(sequence='testing'))])

	while count < (len(all_node_list) - 1):

		if (all_node_list[count][2], all_node_list[count+1][2]) in graph_obj.edges():

			#print 'has edge'
			#print count

			#print all_node_list[count][2]
			#print all_node_list[count+1][2]

			#print '------'
			#print graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]]


			if sequence_name not in graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'].split(','):

				new_seq_list = graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'] + ',' + sequence_name

			else:
				new_seq_list = graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence']
			#print new_seq_list

			graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'] = new_seq_list

			#print 'appending'
		else:

			graph_obj.add_edges_from([(all_node_list[count][2], all_node_list[count+1][2], dict(sequence=sequence_name))])

		count += 1

	graph_obj.graph['add_count'] = add_count_val

	print 'returning'

	return graph_obj

def link_nodes(graph_obj, sequence_name, node_prefix='gn'):
	print 'link nodes 3'
	
	# To make sure all new nodes have unique names, a count is kept

	if 'add_count' in graph_obj.graph.keys():
		#print 'found'
		add_count_val = graph_obj.graph['add_count']
	else:
		#print 'none found'
		graph_obj.graph['add_count'] = 0
		add_count_val = 0

	# This list is a list of all the existing nodes start and stop positions with the node name formatted at start, stop, node_name

	pos_lol = []

	print 'getting positions'

	for node,data in graph_obj.nodes_iter(data=True):
		#print data
		#print data[left_end_name], data[right_end_name], node
		#print 'gggggg'
		#print data

		left_end_name = sequence_name + '_leftend'
		right_end_name = sequence_name + '_rightend'		
		#print '------------'
		#print data
		#print node
		#print sequence_name
		if sequence_name in data['present_in'].split(','):
			#print 'right node'
			#print data

			new_left_pos = abs(int(data[left_end_name]))
			new_right_pos = abs(int(data[right_end_name]))

			if new_left_pos > new_right_pos:
				print 'Doing the switch!'
				temp_left = new_right_pos
				temp_right = new_left_pos
				new_left_pos = temp_left
				new_right_pos = temp_right


			if int(new_left_pos) != 0 and int(new_right_pos) != 0 or node == sequence_name + '_start':
				pos_lol.append([new_left_pos, new_right_pos, node])

	count = 0

	gap_lol = []

	print 'Sorting'

	pos_lol = sorted(pos_lol, key=lambda start_val: start_val[0])

	#for a_list in pos_lol:
	#	print a_list

	sorted_pos_lol = pos_lol

	for a_pos in sorted_pos_lol:
		if a_pos[0] > a_pos[1]:
			print 'ERROR ON LINE 734'

	#print sorted_pos_lol

	# -------------------------------------------------------- Here we determine if there are any gaps in the sequence 
	print 'Look for gaps'

	while count < (len(sorted_pos_lol) - 1):

		# if the end of the current node is not equal to the end of the next node -1 and it is not equal to the start of the next node
			
		current_node_end = sorted_pos_lol[count][1]
		current_node_start = sorted_pos_lol[count][0]
		next_node_start = sorted_pos_lol[count + 1][0]
		next_node_stop = sorted_pos_lol[count + 1][1]

		if current_node_end > next_node_start or current_node_end > next_node_stop or current_node_start > next_node_start or current_node_start > next_node_stop:
			print 'Node order wrong'
			print current_node_start, current_node_end
			print next_node_start, next_node_stop


		if current_node_end != (next_node_start - 1):

			if current_node_end != next_node_start:

				if current_node_start > current_node_end or next_node_start > next_node_stop or current_node_end > next_node_start:

					print "Strange thing --------------------------<"
					print current_node_start, current_node_end
					print next_node_start, next_node_stop

				gap_list = [(current_node_end + 1), ((next_node_start - 1)), sequence_name]
				gap_lol.append(gap_list)

		count +=1


	print 're-orientate'
	for entry in gap_lol:
		if entry[0] > entry[1]:
			print 'wtf'
			print entry
			new_entry_0 = entry[1]
			new_entry_1 = entry[0]
			entry[1] = new_entry_1
			entry[0] = new_entry_0

	# -------------------------------------------------------- Add new nodes

	#print gap_lol

	new_node_list = []

	print 'Add new nodes'

	for new_seq_node in gap_lol:
		#print new_seq_node
		#gap_len = new_seq_node[1] - new_seq_node[0] + 1
		node_dict = {new_seq_node[2] + '_leftend':str(new_seq_node[0]), new_seq_node[2] + '_rightend':str(new_seq_node[1]), 'present_in':sequence_name}

		#print node_dict
		node_name = node_prefix + '_' + str(add_count_val)

		new_seq_node.append(node_name)

		graph_obj.add_node(node_name, node_dict)


		new_node_list.append([new_seq_node[0], new_seq_node[1], node_name])

		add_count_val += 1





	#print pos_lol
	#print gap_lol
	#print new_node_list

	# -------------------------------------------------------- Combine the original node list with the list of gaps generated earlier

	all_node_list = pos_lol + new_node_list
	#print '\n'

	all_node_list = sorted(all_node_list)


	#print 'This is the list*************************'
	#print all_node_list


	# -------------------------------------------------------- Here we add edges to the graph, weaving in the new nodes

	count = 0

	#graph_obj.add_edges_from([('H37Rv_seq_4', 'Aln_region_7', dict(sequence='testing'))])

	print 'Add new edges'

	edges_obj = graph_obj.edges()

	while count < (len(all_node_list) - 1):
	

		if (all_node_list[count][2], all_node_list[count+1][2]) in edges_obj:

			#print 'has edge'
			#print count

			#print all_node_list[count][2]
			#print all_node_list[count+1][2]

			#print '------'
			#print graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]]

	
			if sequence_name not in graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'].split(','):

				new_seq_list = graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'] + ',' + sequence_name

			else:
				new_seq_list = graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence']
			#print new_seq_list

			graph_obj.edge[all_node_list[count][2]][all_node_list[count+1][2]][0]['sequence'] = new_seq_list

			#print 'appending'
		else:

			graph_obj.add_edges_from([(all_node_list[count][2], all_node_list[count+1][2], dict(sequence=sequence_name))])

		count += 1

	print 'Edges added'

	graph_obj.graph['add_count'] = add_count_val


	return graph_obj



def create_genome_alignment_graph(backbone_file, seq_label_dict, input_path_dict, node_prefix=''):

	#print seq_label_dict
	#print input_path_dict
	# Getting graph info

	all_iso_in_graph = ''

	for isolate_name in seq_label_dict.keys():
		all_iso_in_graph = all_iso_in_graph + seq_label_dict[isolate_name] + ','

	all_iso_in_graph = all_iso_in_graph[:-1]

	#print all_iso_in_graph


	# ------------------------------------------------------- Determining node names 

	if len(node_prefix) > 0:
		list_node_prefix = node_prefix.split('_')
		node_new_mid_fix = ''.join(list_node_prefix[2:])


	backbone_lol = input_parser(backbone_file)

	genome_network = nx.MultiDiGraph()

	genome_network.graph['isolates'] = all_iso_in_graph

	als_seq_number = len(seq_label_dict)
	#print '================================'
	#print seq_label_dict
	aln_number = 0

	headder_line = backbone_lol[0]

	#print headder_line

	backbone_lol_headless = backbone_lol[1:]


	# Replacing generic names based on dict
	new_headder_line = []

	for headder_seq_name in headder_line:
		seq_name_list = headder_seq_name.split('_')
		#print seq_name_list
		for seq_key in seq_label_dict.keys():
			#print seq_key
			if seq_name_list[0] == seq_key:
				new_headder_name = seq_label_dict[seq_key] + '_' + seq_name_list[1]
				new_headder_line.append(new_headder_name)

	#print new_headder_line
	#print '---------'

	# Adding the aligned region nodes in this segment of code
	for aligned_region in backbone_lol_headless:
		aln_number += 1
		#print aligned_region
		node_dict = {}

		# creating length dict
		length_dict = {}


		row_number = 0

		# Adding isolate start/stop aka left/right to the node dict
		while row_number < len(headder_line):
			#pos_label = headder_line[row_number]
			#print new_headder_line
			#print node_dict[new_headder_line[row_number]]
			#print aligned_region[row_number]
			#print new_headder_line
			node_dict[new_headder_line[row_number]] = aligned_region[row_number]

			row_number += 1
		#print node_dict

		# creating length dict
		length_dict = {}

		for seq_name in seq_label_dict:
			left_val = node_dict[seq_label_dict[seq_name] + '_leftend']
			right_val = node_dict[seq_label_dict[seq_name] + '_rightend']

			#print left_val, right_val

			seq_aln_length = abs(int(left_val) - int(right_val))


			# if int(left_val) != 0 and int(right_val) != 0 and int(left_val) != 1 and int(right_val) != 1:   (Changed 15 july 2016 to allow nodes starting at 1)
			if int(left_val) != 0 and int(right_val) != 0:
				#print 'hmmmm'
				#print aln_number
				#print left_val, right_val
				seq_aln_length = seq_aln_length + 1
				

			length_dict[seq_label_dict[seq_name]] = seq_aln_length



		present_in_str = ''


		for seq_name in length_dict.keys():
			if int(length_dict[seq_name]) != 0:
				present_in_str = present_in_str + seq_name
				present_in_str = present_in_str + ','

		
		present_in_str = present_in_str[:-1]

		node_dict['present_in'] = present_in_str
		#print present_in_str


		# Assign the length of the node

		iso_for_length = present_in_str.split(',')[0]

		#node_dict['length'] = length_dict[iso_for_length]


		# ----------- ----------- ----------- ----------- ----------- Adding prefix if needed to node name -----------

		if len(node_prefix) > 0:
			node_name = "Aln_" + node_new_mid_fix + "_" +  str(aln_number)
		else:
			node_name = "Aln_" + str(aln_number)

		add_node = True

		if len(node_dict['present_in'].split(',')) == 0:
			add_node = False

		test_node_dict = {}

		if add_node == True:
			#print '=========='
			#print node_name
			#print node_dict


			genome_network.add_node(node_name, node_dict)

	# ---------- ----------- ----------- ----------- -----------  Add seq start and stop ends ----------- -----------

	for seq_name in input_path_dict.keys():
		seq_fasta = input_parser(input_path_dict[seq_name])
		#print seq_fasta[0].keys()
		#print len(seq_fasta[0]['DNA_seq'])
		seq_total_length = len(seq_fasta[0]['DNA_seq'])

		start_node_name = seq_name + '_start'
		stop_node_name = seq_name + '_stop'

		genome_network.add_node(start_node_name)
		genome_network.node[start_node_name][seq_name + '_leftend'] = '0'
		genome_network.node[start_node_name][seq_name + '_rightend'] = '0'
		genome_network.node[start_node_name]['present_in'] = seq_name

		genome_network.add_node(stop_node_name)
		genome_network.node[stop_node_name][seq_name + '_leftend'] = str(seq_total_length)
		genome_network.node[stop_node_name][seq_name + '_rightend'] = str(seq_total_length)
		genome_network.node[stop_node_name]['present_in'] = seq_name


	# And now to add the edges...

	for seq_name in seq_label_dict.keys():
		#print seq_label_dict
		genome_network = link_nodes(genome_network, seq_label_dict[seq_name])




	return genome_network


def add_sequences_to_graph(graph_obj, paths_dict):
	
	print 'Adding sequences'

	for node,data in graph_obj.nodes_iter(data=True):

		seq_source = data['present_in'].split(',')[0]
		is_reversed = False
		is_comp = False

		if len(seq_source) < 1:
			print 'No present_in current node'
			print node

		else:
			ref_seq = input_parser(paths_dict[seq_source])[0]['DNA_seq']

			# Check orientation
			if int(data[seq_source + '_leftend']) < 0:
				is_reversed = True

			seq_start = abs(int(data[seq_source + '_leftend']))
			seq_end = abs(int(data[seq_source + '_rightend']))

			if seq_start > seq_end:
				new_seq_start = seq_end
				new_seq_end = seq_start
				seq_end = new_seq_end
				seq_start = new_seq_start


			seq_start = seq_start - 1

			node_seq = ref_seq[seq_start:seq_end].upper()

			if is_reversed == True:
				print 'seq was rev comp' + node
				node_seq = reverse_compliment(node_seq)

			graph_obj.node[node]['sequence'] = node_seq

	return graph_obj

def make_circular(graph_obj, seq_name):

	start_node_name = seq_name + '_start'
	end_node_name = seq_name + '_stop'

	graph_obj.add_edges_from([(start_node_name, end_node_name, dict(sequence=seq_name))])

	return graph_obj

def check_isolates_in_region(graph_obj, start_pos, stop_pos, reference_name, threshold=1.0, return_dict=False, simmilarity_measure='percentage'):
	'''Retrieve the nodes from a graph spanning a region'''


	print start_pos, stop_pos

	#print '\nchecking iso in region'

	graph_isolate_list = graph_obj.graph['isolates'].split(',')

	expected_ref_length = int(stop_pos) - int(start_pos)
	expected_ref_length = expected_ref_length + 1

	if expected_ref_length < 0:
		print 'POSSIBLE ERROR, start larger than stop'
		print int(stop_pos), int(start_pos)
		print expected_ref_length

	#print expected_ref_length
	#print reference_name

	# Extract subgraph
	node_leftend_label = reference_name + '_leftend'
	node_rightend_label = reference_name + '_rightend'

	# Identify the nodes that contain the start and stop positions

	#print int(start_pos)
	#print int(stop_pos)
	print reference_name

	for node,data in graph_obj.nodes_iter(data=True):
		if reference_name in data['present_in']:

			if int(data[node_leftend_label]) > 0:

				if abs(int(data[node_leftend_label])) <= int(start_pos) <= abs(int(data[node_rightend_label])):
					start_node = node
					print 'Ping'
					#print data[node_leftend_label]
					print start_pos
					#print data[node_rightend_label]
					start_overlap =  bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) <= int(stop_pos) <= abs(int(data[node_rightend_label])) or abs(int(data[node_leftend_label])) >= int(stop_pos) >= abs(int(data[node_rightend_label])):
					print 'ping'
					stop_node = node
					print stop_pos
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])
					#stop_overlap = int(stop_pos) - int(data[node_leftend_label]) + 1

			if int(data[node_leftend_label]) < 0:
				
				if abs(int(data[node_leftend_label])) >= int(start_pos) >= abs(int(data[node_rightend_label])):
					start_node = node
					print 'Ping'
					print 'neg node'
					#print data[node_leftend_label]
					print start_pos
					#print data[node_rightend_label]
					start_overlap =  bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) >= int(stop_pos) >= abs(int(data[node_rightend_label])) or abs(int(data[node_leftend_label])) >= int(stop_pos) >= abs(int(data[node_rightend_label])):
					print 'ping'
					print 'neg node'
					stop_node = node
					print stop_pos
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])
					#stop_overlap = int(stop_pos) - int(data[node_leftend_label]) + 1



	#print start_node
	#print start_overlap
	#print stop_node
	#print stop_overlap
	#print '\n'

	# Dealing with genes that occur at the start and stop nodes of the graph... needs a propper solution
	# Caused by the lack of a 'length' attribute in the start and stop node
	# Temp fix = just return a nope.

	if start_node.split('_')[-1] == 'start' or stop_node.split('_')[-1] == 'stop':
		if return_dict == True:
			# THIS WILL FAIL maybe...
			return {reference_name:1}
		else:
			return [reference_name]

	# Calculating overlap

	start_node_pb_length = bp_length_node(graph_obj.node[start_node])
	stop_node_pb_length = bp_length_node(graph_obj.node[stop_node])

	#print start_node_pb_length


	start_node_nonoverlap = start_node_pb_length - start_overlap
	stop_node_nonoverlap = stop_node_pb_length - stop_overlap

	total_overlap = start_overlap + stop_overlap

	# Finding shortest path between nodes (IF start and stop nodes are not the same)
	nodes_in_path = []
	ref_nodes_in_path = []
	alt_nodes_in_path = []

	'''
	if start_node != stop_node:
		for path in nx.shortest_simple_paths(graph_obj, source=start_node, target=stop_node):
			#print path
			for path_node in path:
				if path_node not in nodes_in_path:
					nodes_in_path.append(path_node)
	else:
		nodes_in_path.append(start_node)
		nodes_in_path.append(stop_node)

	'''

	if start_node != stop_node:
		for path in nx.shortest_path(graph_obj, source=start_node, target=stop_node):
			#print path
			nodes_in_path.append(path)
	else:
		nodes_in_path.append(start_node)
		nodes_in_path.append(stop_node)

	#print 'nodes in path'
	#print nodes_in_path


	for path_node in nodes_in_path:
		if reference_name in graph_obj.node[path_node]['present_in']:
			ref_nodes_in_path.append(path_node)
		else:
			alt_nodes_in_path.append(path_node)

	#print ref_nodes_in_path
	#print alt_nodes_in_path

	total_alt_length = 0
	total_ref_length = 0

	for alt_node in alt_nodes_in_path:
		total_alt_length = total_alt_length + bp_length_node(graph_obj.node[alt_node])

	

	for ref_node in ref_nodes_in_path:
		if ref_node != start_node and ref_node != stop_node:
			total_ref_length = total_ref_length + bp_length_node(graph_obj.node[ref_node])

	total_ref_length = total_ref_length + start_overlap + stop_overlap

	#print total_alt_length
	#print total_ref_length

	iso_diff_dict = {}
	iso_sim_dict = {}

	for node_iso in graph_isolate_list:
		#print node_iso
		iso_diff_dict[node_iso] = 0
		iso_sim_dict[node_iso] = 0


	for a_node in nodes_in_path:
		for a_isolate in graph_isolate_list:
			if a_isolate in graph_obj.node[a_node]['present_in'] and reference_name in graph_obj.node[a_node]['present_in']:
				iso_sim_dict[a_isolate] = iso_sim_dict[a_isolate] + bp_length_node(graph_obj.node[a_node])

			elif a_isolate in graph_obj.node[a_node]['present_in'] and reference_name not in graph_obj.node[a_node]['present_in']:
				iso_diff_dict[a_isolate] = iso_diff_dict[a_isolate] + bp_length_node(graph_obj.node[a_node])

	for a_iso in graph_isolate_list:
		if a_iso in graph_obj.node[start_node]['present_in']:
			iso_sim_dict[a_iso] = iso_sim_dict[a_iso] - start_node_nonoverlap
		if a_iso in graph_obj.node[stop_node]['present_in']:
			iso_sim_dict[a_iso] = iso_sim_dict[a_iso] - stop_node_nonoverlap


	#print iso_sim_dict
	#print iso_diff_dict

	iso_sim_score_dict = {}

	# Calculating percentage similarity 

	if simmilarity_measure == 'percentage':

		for a_isolate in graph_isolate_list:
			iso_sim_score_dict[a_isolate] = float(iso_sim_dict[a_isolate] - iso_diff_dict[a_isolate]) / float(iso_sim_dict[reference_name])

	if simmilarity_measure == 'levenshtein':

		levenshtein(seq_1, seq_2)

	# Applying filter
	ret_list = []

	for a_isolate in graph_isolate_list:
		if iso_sim_score_dict[a_isolate] >= threshold:
			ret_list.append(a_isolate)

	#print iso_sim_score_dict

	if return_dict == True:
		return iso_sim_score_dict
	else:
		return ret_list

def convert_coordinates(graph_obj, q_start, q_stop, ref_iso, query_iso):

	# Find the right node

	node_leftend_label = ref_iso + '_leftend'
	node_rightend_label = ref_iso + '_rightend'

	query_leftend_label = query_iso + '_leftend'
	query_rightend_label = query_iso + '_rightend'

	ref_node_left = ''
	ref_node_right = ''

	query_node_left = ''
	query_node_right = ''

	for node,data in graph_obj.nodes_iter(data=True):
		if ref_iso in data['present_in'] and query_iso in data['present_in']:
			if int(data[node_leftend_label]) < int(q_start) and int(q_stop) < int(data[node_rightend_label]):
				#print 'local found!!!!'
				ref_node_left = data[node_leftend_label]
				ref_node_right = data[node_rightend_label]
				query_node_left = data[query_leftend_label]
				query_node_right = data[query_rightend_label]

				#print int(ref_node_left) - int(ref_node_right)
				#print int(query_node_left) - int(query_node_right)

				ref_low_num = abs(int(ref_node_left))
				if abs(int(ref_node_left)) > abs(int(ref_node_right)):
					ref_low_num = abs(int(ref_node_right))

				#print ref_low_num

				query_low_num = abs(int(query_node_left))
				if abs(int(query_node_left)) > abs(int(query_node_right)):
					query_low_num = abs(int(query_node_right))

				#print query_low_num

				conversion_factor = ref_low_num - query_low_num

				#print conversion_factor

				new_r_start = int(q_start) - conversion_factor
				new_r_stop = int(q_stop) - conversion_factor

				#print new_r_start
				#print new_r_stop

				return {query_iso + '_leftend':new_r_start, query_iso + '_rightend':new_r_stop}


def convert_coordinate(graph_obj, q_position, ref_iso, query_iso):

	# Find the right node

	node_leftend_label = ref_iso + '_leftend'
	node_rightend_label = ref_iso + '_rightend'

	query_leftend_label = query_iso + '_leftend'
	query_rightend_label = query_iso + '_rightend'


	for node,data in graph_obj.nodes_iter(data=True):

		# Look at all nodes that contain both the isolates

		if ref_iso in data['present_in'] and query_iso in data['present_in']:

			# find the node

			if int(data[node_leftend_label]) > 0:
				#print 'positive node search'
				if abs(int(data[node_leftend_label])) <= abs(int(q_position)) <= abs(int(data[node_rightend_label])):
					print 'positive local found!!!!'
					#print data[node_leftend_label]
					#print data[node_rightend_label]

					#print data

					ref_node_left = data[node_leftend_label]
					ref_node_right = data[node_rightend_label]
					query_node_left = data[query_leftend_label]
					query_node_right = data[query_rightend_label]

					#print int(ref_node_left) - int(ref_node_right)
					#print int(query_node_left) - int(query_node_right)

					ref_low_num = abs(int(ref_node_left))
					
					if abs(int(ref_node_left)) > abs(int(ref_node_right)):
						ref_low_num = abs(int(ref_node_right))

					#print ref_low_num

					query_low_num = abs(int(query_node_left))
					if abs(int(query_node_left)) > abs(int(query_node_right)):
						query_low_num = abs(int(query_node_right))

					#print query_low_num

					conversion_factor = ref_low_num - query_low_num

					#print conversion_factor

					new_r_start = int(q_position) - conversion_factor

					#print new_r_start
					#print new_r_stop

					return {query_iso:new_r_start}

			if int(data[node_leftend_label]) < 0:
				#print 'negative ref search'
				if abs(int(data[node_leftend_label])) >= abs(int(q_position)) >= abs(int(data[node_rightend_label])):
					print 'Negative local found'
					print q_position
					print data[node_leftend_label], data[node_rightend_label]

					#print data

					ref_node_left = data[node_leftend_label]
					ref_node_right = data[node_rightend_label]
					query_node_left = data[query_leftend_label]
					query_node_right = data[query_rightend_label]

					print query_node_left, query_node_right

					#print int(ref_node_left) - int(ref_node_right)
					#print int(query_node_left) - int(query_node_right)

					ref_node_smallest_val = ref_node_right

					if abs(int(ref_node_left)) < abs(int(ref_node_right)):
						ref_node_smallest_val = ref_node_left

					query_high_num = abs(int(query_node_right))

					if abs(int(query_node_left)) > abs(int(query_node_right)):
						query_high_num = abs(int(query_node_right))

					#print query_low_num

					conversion_factor = abs(int(q_position)) - abs(int(ref_node_smallest_val))

					print 'conversion_factor'
					print conversion_factor
					print query_high_num
					new_r_start = query_high_num - conversion_factor

					print '-----'
					print new_r_start
					#print new_r_stop

					if int(query_node_left) < 0:
						new_r_start = new_r_start * (-1)

					return {query_iso:new_r_start}



	else:
		return 'pos not found'


def import_gtf_dict_to_massive_dict(gtf_dict):

	all_dict = {}
	#print gtf_dict

	for isolate in gtf_dict.keys():

		curr_gtf = input_parser(gtf_dict[isolate], parse_as='gtf')
		
		for entry in curr_gtf:
			#print entry
			pos_deets = isolate + ',' + entry[3] + ',' + entry[4] + ',' + entry[5]
			if 'gene_id' not in entry[8].keys():
				print 'key error'
				print entry
				print isolate
				quit()
			else:
				gene_name = entry[8]['gene_id']
				all_dict[gene_name] = pos_deets

	#print all_dict
	return all_dict


def fasta_alignment_to_subnet(fasta_aln_file, true_start={}, node_prefix='X', orientation={}, re_link_nodes=True):
	'''New and improved conversion function. Needs to be modified to work with inverted seq still'''
	# Created 11/01/2017
	#print true_start


	aln_lol = input_parser(fasta_aln_file)

	offset_dict = {}
	gap_dict = {}
	all_isolate_list = []
	# Get the length of each of the sequences without gaps
	seq_len_dict = {}

	for fasta_entry in aln_lol:

		#print '$$$$$$'
		#print fasta_entry
		#print fasta_entry['gene_details']
		strip_seq = fasta_entry['DNA_seq'].replace('-','')

		seq_len_dict[fasta_entry['gene_details']] = len(strip_seq)

		all_isolate_list.append(fasta_entry['gene_details'])
		fasta_entry['DNA_seq'] = fasta_entry['DNA_seq'].upper()

	#print seq_len_dict

	# Setting the orientation dict for later
	if len(orientation) == 0:
		#print 'Orientation dict empty'
		for a_isolate in all_isolate_list:
			orientation[a_isolate] = '+'

	# Creating graph. May remove later 
	local_node_network = nx.MultiDiGraph()
	local_node_network.graph['isolates'] = ','.join(all_isolate_list) 

	# Getting needed vals
	total_len = len(aln_lol[0]['DNA_seq'])
	count = 0


	# Creating block list
	block_list = []
	block_ends_list = []

	while count < total_len:

		block_dict = {}
		#block_dict['pos'] = count

		for a_seq in aln_lol:

			#offset_dict[a_seq['gene_details']] = 0
			#gap_dict[a_seq['gene_details']] = []
			


			current_keys = block_dict.keys()
			if a_seq['DNA_seq'][count] in current_keys:
				block_dict[a_seq['DNA_seq'][count]].append(a_seq['gene_details'])
			else:
				block_dict[a_seq['DNA_seq'][count]] = [a_seq['gene_details']]

		block_list.append(block_dict)

		count += 1

	#print '------------'
	#print block_list

	relative_pos_dict = {}

	#print orientation

	# Making sure the start pos is correct
	
	for an_isolate in all_isolate_list:
		true_start[an_isolate] = true_start[an_isolate] - 1

	# Adding additional info and indexing

	last_bline = {'blocklist':[]}

	bpos_count = 1

	# Working along the black from here ----------------------
	for bline in block_list:

		for an_isolate in all_isolate_list:

			true_start[an_isolate] += 1


			for base in bline.keys():
				if base == '-' and an_isolate in bline[base]:
					#print an_isolate
					#print '***********'
					true_start[an_isolate] = true_start[an_isolate] - 1



		bline['global_pos'] = bpos_count
		bline['relative_pos'] = copy.deepcopy(true_start)

		# See if this bline is a new block / node set
		block_list = []
		for base in bline.keys():
			if base != '-' and base != 'global_pos' and base != 'relative_pos':

				block_list.append(bline[base])

				# make [ordered list [ ] ]

		block_list = sorted(block_list)

		bline['blocklist'] = block_list

		


		# Here we call blocks
		if bline['blocklist'] != last_bline['blocklist']:
			'''
			print '\n'
			print 'New block' 

			print 'last block end'
			print last_bline
			print 'new lock start'
			print bline
			print '\n'
			'''
			block_ends_list.append(last_bline)
			block_ends_list.append(bline)
			



		#print 'last line'
		#print last_bline
		#print 'this line'
		#print bline

		last_bline = bline

		bpos_count += 1

	# Last block
	block_ends_list.append(last_bline)

	#print '\n'

	# Trimming off that forst start block
	block_ends_list = block_ends_list[1:]


	# Here we start getting the nodes paired and into the correct format
	# Start and stop are in pairs

	start_block_list = []
	end_block_list = []

	count = 0
	for end_block in block_ends_list:
		#print end_block, '\t', end_block['blocklist']
		if is_even(count) == True:
			start_block_list.append(end_block)
		else:
			end_block_list.append(end_block)

		count += 1


	# Now we pair the stop and start 
	# Also begin adding nodes to the local_node_network
	#print '\n'
	count = 0
	node_count = 1
	while count < len(start_block_list):
		#print start_block_list[count], end_block_list[count]

		if start_block_list[count]['blocklist'] != end_block_list[count]['blocklist']:
			print "ERROR IN BLOCK LIST PAIRS"

		for block_group in start_block_list[count]['blocklist']:
			# Each of these becomes a node
			#print block_group
			new_node_name = node_prefix + '_' + str(node_count)
			block_group_string = ",".join(block_group)

			local_node_network.add_node(new_node_name, present_in = block_group_string)

			# Adding start / stop for node
			# this is where we do the +/- thing. 
			for block_isolate in block_group:
				if orientation[block_isolate] == '+':
				
					local_node_network.node[new_node_name][block_isolate + '_leftend'] = str(start_block_list[count]['relative_pos'][block_isolate])
					local_node_network.node[new_node_name][block_isolate + '_rightend'] = str(end_block_list[count]['relative_pos'][block_isolate])
				
				elif orientation[block_isolate] == '-':
					
					#So,
					# IF the positions are 1-10, 11-15, 15-40, what would the reverse positions be?
					# len node - pos - 1
					# So, we need the total length of the nodes, found in seq_len_dict

					local_node_network.node[new_node_name][block_isolate + '_leftend'] = '-' + str(int(seq_len_dict[block_isolate]) - start_block_list[count]['relative_pos'][block_isolate] - 1)
					local_node_network.node[new_node_name][block_isolate + '_rightend'] = '-' + str(int(seq_len_dict[block_isolate]) - end_block_list[count]['relative_pos'][block_isolate] - 1)
			
				else:
					print "ORIENTATION MISSING"
					print orientation
					quit()

			node_count += 1

		count += 1

	if re_link_nodes == True:
		for a_isolate in all_isolate_list:
			local_node_network = link_nodes(local_node_network, a_isolate, node_prefix='gn')

	#print '----------------'
	return local_node_network

	#print 'ORIENTATION'
	#print orientation

	print '---------------------------- block list created ----------------------------'


def create_ungapped_graph(in_graph, seq_fasta_paths_dict):

	print 'creating ungapped graph'

	realign_node_list = []

	iso_list = in_graph.graph['isolates'].split(',')

	for node,data in in_graph.nodes_iter(data=True):

		if len(data['present_in'].split(',')) > 1:

			#print node

			realign_node_list.append(node)

			#local_node_realign(in_graph, seq_fasta_paths_dict)

	count = 0

	for a_node in realign_node_list:
		print count
		print a_node
		in_graph = local_node_realign_fast(in_graph, a_node, seq_fasta_paths_dict)
		count += 1

	# Saving unlinked graph
	nx.write_graphml(in_graph, 'intermediate_split_unlinked.xml')

	for iso in iso_list:
		in_graph = link_nodes(in_graph, iso)

	return in_graph

def local_node_realign_fast(in_graph, node_ID, seq_fasta_paths_dict):

	print 'Fast local node realign: ' + node_ID

	in_graph = nx.MultiDiGraph(in_graph)

	node_data_dict = in_graph.node[node_ID]

	# Make temp fasta file and record the start positions into a dict

	node_seq_start_pos = {}

	# Store the orientation of the sequences in the node to pass to the fasta to graph conversion function
	orientation_dict = {}

	temp_fasta_file = open('temp_unaligned.fasta', 'w')

	for node_isolate in node_data_dict['present_in'].split(','):
		iso_full_seq = input_parser(seq_fasta_paths_dict[node_isolate])[0]['DNA_seq'].upper()

		#print node_isolate
		#print '----------------+------------------'

		''' Currenty only rev comp sequences are seen in the BBone file, represented by a - but not reversed start / stop '''
		#print node_data_dict

		if int(node_data_dict[node_isolate + '_leftend']) >= -1:
			orientation_dict[node_isolate] = '+'
		else:
			orientation_dict[node_isolate] = '-'
		

		if int(node_data_dict[node_isolate + '_leftend']) <= int(node_data_dict[node_isolate + '_rightend']):
			#print 'not reversed'
			#print 'Not compliment'
			iso_node_seq = iso_full_seq[int(node_data_dict[node_isolate + '_leftend']) - 1:int(node_data_dict[node_isolate + '_rightend']) ]
			node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])


		else:
			#print 'reversed'
			#print abs(int(node_data_dict[node_isolate + '_rightend']))
			#print abs(int(node_data_dict[node_isolate + '_leftend']))

			iso_node_seq = iso_full_seq[abs(int(node_data_dict[node_isolate + '_leftend'])) - 1:abs(int(node_data_dict[node_isolate + '_rightend'])) ]
			node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])
			iso_node_seq = reverse_compliment(iso_node_seq)

		#print 'seq'
		#print iso_node_seq

		temp_fasta_file.write('>' + node_isolate + '\n')
		temp_fasta_file.write(iso_node_seq + '\n')
		node_seq_len_est = len(iso_node_seq)
		#print 'passed start positions'
		#print node_seq_start_pos

	temp_fasta_file.close()

	if local_aligner == 'muscle':
		print 'conducting muscle alignment'
		muscle_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'mafft':
		print 'conducting mafft alignment'
		print node_seq_len_est
		mafft_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'clustalo':
		print 'conducting clustal alignment'
		clustalo_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'kalign':
		print node_seq_len_est
		if node_seq_len_est < 7:
			print 'conducting mafft alignment'
			mafft_alignment('temp_unaligned.fasta', 'temp_aln.fasta')
		else:
			print 'conducting kalign alignment'
			kalign_alignment('temp_unaligned.fasta', 'temp_aln.fasta')


	#fasta_alignment_to_bbone('temp_aln.fasta', 'temp_aln', true_start=node_seq_start_pos)

	new_subgraph = fasta_alignment_to_subnet('temp_aln.fasta', true_start=node_seq_start_pos, node_prefix=node_ID, orientation=orientation_dict, re_link_nodes=False)

	#nx.write_graphml(new_subgraph, 'test_new_subg_327.xml')



	# -------------------------------- Here we take the new graph for the aligned node, and replace the origional node with it



	iso_list = new_subgraph.graph['isolates'].split(',')


	# remove the old node

	in_graph.remove_node(node_ID)

	# add new subgraph to the old graph

	in_graph.add_nodes_from(new_subgraph.nodes(data=True))

	in_graph.add_edges_from(new_subgraph.edges(data=True))

	new_merged_graph = in_graph
	
	#nx.write_graphml(new_merged_graph, 'test_merged_graphs_linked.xml')

	return new_merged_graph


def local_node_realign(in_graph, node_ID, seq_fasta_paths_dict):

	print 'local node realign: ' + node_ID

	in_graph = nx.MultiDiGraph(in_graph)

	node_data_dict = in_graph.node[node_ID]

	# Make temp fasta file and record the start positions into a dict

	node_seq_start_pos = {}

	# Store the orientation of the sequences in the node to pass to the fasta to graph conversion function
	orientation_dict = {}

	temp_fasta_file = open('temp_unaligned.fasta', 'w')

	for node_isolate in node_data_dict['present_in'].split(','):
		iso_full_seq = input_parser(seq_fasta_paths_dict[node_isolate])[0]['DNA_seq'].upper()

		#print node_isolate
		#print '----------------+------------------'

		''' Currenty only rev comp sequences are seen in the BBone file, represented by a - but not reversed start / stop '''
		#print node_data_dict

		if int(node_data_dict[node_isolate + '_leftend']) > 0:
			orientation_dict[node_isolate] = '+'
		else:
			orientation_dict[node_isolate] = '-'
		
		if int(node_data_dict[node_isolate + '_leftend']) < int(node_data_dict[node_isolate + '_rightend']):
			#print 'not reversed'
			if int(node_data_dict[node_isolate + '_leftend']) > 0:
				#print 'Not compliment'
				iso_node_seq = iso_full_seq[int(node_data_dict[node_isolate + '_leftend']) - 1:int(node_data_dict[node_isolate + '_rightend']) ]
				node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])

			else:
				#print 'Complimented'
				iso_node_seq = iso_full_seq[int(node_data_dict[node_isolate + '_leftend']) - 1:int(node_data_dict[node_isolate + '_rightend']) ]
				node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])
				iso_node_seq = compliment_DNA(iso_node_seq)

		else:
			#print 'reversed'
			#print abs(int(node_data_dict[node_isolate + '_rightend']))
			#print abs(int(node_data_dict[node_isolate + '_leftend']))

			if int(node_data_dict[node_isolate + '_leftend']) > 0:
				#print 'Not complimented'
				iso_node_seq = iso_full_seq[abs(int(node_data_dict[node_isolate + '_rightend'])) - 1:abs(int(node_data_dict[node_isolate + '_leftend'])) ]
				node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_rightend'])
			else:
				#print 'Compliment'
				iso_node_seq = iso_full_seq[abs(int(node_data_dict[node_isolate + '_leftend'])) - 1:abs(int(node_data_dict[node_isolate + '_rightend'])) ]
				node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])
				iso_node_seq = reverse(iso_node_seq)
				iso_node_seq = compliment_DNA(iso_node_seq)

		#print 'seq'
		#print iso_node_seq

		temp_fasta_file.write('>' + node_isolate + '\n')
		temp_fasta_file.write(iso_node_seq + '\n')
		#print 'passed start positions'
		#print node_seq_start_pos

	temp_fasta_file.close()

	if local_aligner == 'muscle':
		print 'conducting muscle alignment'
		muscle_alignment('temp_unaligned.fasta', 'temp_aln.fasta')
	if local_aligner == 'clustalo':
		print 'conducting clustal alignment'
		clustalo_alignment('temp_unaligned.fasta', 'temp_aln.fasta')


	#fasta_alignment_to_bbone('temp_aln.fasta', 'temp_aln', true_start=node_seq_start_pos)

	new_subgraph = fasta_alignment_to_subnet('temp_aln.fasta', true_start=node_seq_start_pos, node_prefix=node_ID, orientation=orientation_dict)

	nx.write_graphml(new_subgraph, 'test_new_subg_327.xml')



	# -------------------------------- Here we take the new graph for the aligned node, and replace the origional node with it



	iso_list = new_subgraph.graph['isolates'].split(',')


	# remove the old node

	in_graph.remove_node(node_ID)

	# add new subgraph to the old graph

	in_graph.add_nodes_from(new_subgraph.nodes(data=True))

	in_graph.add_edges_from(new_subgraph.edges(data=True))

	for iso in iso_list:
		print 'tag'
		print iso
		in_graph = link_nodes(in_graph, iso)

	new_merged_graph = in_graph
	
	#nx.write_graphml(new_merged_graph, 'test_merged_graphs_linked.xml')

	return new_merged_graph

# ---------------------------------------------------- Alignment functions 


def kalign_alignment(fasta_unaln_file, out_aln_name):

	kalign_command_call = [path_to_kalign, '-i', fasta_unaln_file, '-o', out_aln_name, '-f', 'fasta', '-q']

	return call(kalign_command_call)

def muscle_alignment(fasta_unaln_file, out_aln_name):


	muscle_command_call = [path_to_muscle, '-in', fasta_unaln_file, '-out', out_aln_name]

	return call(muscle_command_call)

def clustalo_alignment(fasta_unaln_file, out_aln_name):

	clustalo_command_call = [path_to_clustal, '-i', fasta_unaln_file, '-o', out_aln_name, '--force']

	return call(clustalo_command_call)

def mafft_alignment(fasta_unaln_file, out_aln_name):

	out_temp_fa = open(out_aln_name, 'w')

	call([path_to_mafft, '--retree', '2', '--maxiterate', '2', '--quiet', '--thread', '-1', fasta_unaln_file], stdout=out_temp_fa)

def progressiveMauve_alignment(fasta_path_list, out_aln_name):

	progressiveMauve_call = [path_to_progressiveMauve, '--output=globalAlignment_'+out_aln_name] + fasta_path_list

	return call(progressiveMauve_call)

# ---------------------------------------------------- Utility functions 

def nodes_connected(u, v, graph_obj):
	return u in graph_obj.neighbors(v)

def convert_transcriptome_to_aln_input(trans_fasta_path, out_name):

	chrom_name = 'MTB_Pan1'

	break_length = 120

	fasta_lol = input_parser(trans_fasta_path)

	#print fasta_lol[1]

	# Make masked break string

	break_str = 'N' * break_length

	#print break_str


	out_file_fasta = open(out_name + '.fasta', 'w')
	out_file_gtf = open(out_name + '.gtf', 'w')

	out_file_fasta.write('>' + chrom_name + '\n')

	# Writing the fasta and gft files

	cur_start = 0

	cur_stop = 0

	for trans_seq in fasta_lol:	
		
		seq_len = len(trans_seq['DNA_seq'])

		cur_start = cur_start + break_length

		cur_stop = cur_start + seq_len


		out_file_fasta.write(break_str)

		out_file_fasta.write(trans_seq['DNA_seq'])


		gtf_line = chrom_name + '\t' + 'gg_pan_genome' + '\t' + 'exon' + '\t' + str(cur_start) + '\t' + str(cur_stop) + '\t' + '.' + '\t' + '+' + '\t' + '0' + '\t' + 'gene_id "' + trans_seq['gene_details'] +'"; transcript_id "' + trans_seq['gene_details'] + '";\n'

		cur_start = cur_stop

		out_file_gtf.write(gtf_line)


# --------------------- PanGenome related

def extract_pan_genome(graph_obj, gtf_dict, out_file_name):
	'''Currently not functioning'''
	#isolate_list = ['H37Rv', 'H37Ra', 'CDC1551'] <-- needs to come from the graph['isolates'] list
	added_list = []

	outfile_obj = open(out_file_name + '.fa', 'w')

	for isolate in isolate_list:
		gtf_lol = input_parser(gtf_dict[isolate])
		for entry in gtf_lol:

			if entry[2] == 'exon':	
				found_in_string = check_isolates_in_region(graph_obj, entry[3], entry[4], isolate)
				found_in_list = found_in_string.split(',')
				if len(list(set(found_in_list) & set(added_list))) < 1:
					curr_seq = extract_seq_region(graph_obj, entry[3], entry[4], isolate)
					curr_seq_headder = '>' + entry[8] + found_in_string + '\n'
					print curr_seq_headder
					outfile_obj.write(curr_seq_headder)
					outfile_obj.write(curr_seq)
					outfile_obj.write('\n')

		added_list.append(isolate)

def extract_pan_genome_csv(graph_obj, gtf_dict, out_file_name, hom_threshold=1.0, refseq=''):

	isolate_list = []

	for a_key in gtf_dict.keys():
		isolate_list.append(a_key)

	if len(refseq) > 1:
		isolate_list.remove(refseq)
		isolate_list = [refseq] + isolate_list

	added_list = []

	outfile_obj = open(out_file_name + '.csv', 'w')
	
	csv_headder = 'gene'

	for iso in isolate_list:
		csv_headder = csv_headder + ',' + iso

	csv_headder = csv_headder + '\n'

	outfile_obj.write(csv_headder)

	for isolate in isolate_list:
		gtf_lol = input_parser(gtf_dict[isolate], parse_as='gtf')
		for entry in gtf_lol:

			if entry[2] == 'exon' or entry[2] == 'gene':
				#print '_____________________'	
				#print entry[3], entry[4]
				#print entry
				if int(entry[3]) < int(entry[4]):
					ent_start_pos = entry[3]
					ent_stop_pos = entry[4]

				if int(entry[3]) > int(entry[4]):
					ent_start_pos = entry[4]
					ent_stop_pos = entry[3]
				#print ent_start_pos, ent_stop_pos

				found_in_list = check_isolates_in_region(graph_obj, ent_start_pos, ent_stop_pos, isolate, threshold=hom_threshold)
				
				if len(list(set(found_in_list) & set(added_list))) < 1:
					line_str = csv_headder

					if 'ID' in entry[8].keys():
						curr_gene = entry[8]['ID']

					if 'gene_id' in entry[8].keys():
						curr_gene = entry[8]['gene_id']

					#print found_in_list
					#print line_str
					line_str = line_str.replace('gene',curr_gene)
					#print line_str
					line_str = line_str.replace(isolate,'1')
					#print line_str


					for iso in isolate_list:
						if iso in found_in_list:
							line_str = line_str.replace(iso,'1')
						else:
							line_str = line_str.replace(iso,'0')



					outfile_obj.write(line_str)


		added_list.append(isolate)

def create_fasta_from_pangenome_csv(pg_csv, gtf_dict, seq_file_dict, out_name):
	in_pg_obj = open(pg_csv, "r")
	out_transcriptome = open(out_name + '.fasta', "w")
	reader = csv.reader(in_pg_obj)
	next(reader, None)

	genome_dict = {}

	for iso_seq_dict in seq_file_dict[1].keys():
		
		fasta_seq_dict = input_parser(seq_file_dict[1][iso_seq_dict])
		genome_dict[iso_seq_dict] = fasta_seq_dict[0]['DNA_seq']


	all_anno_dict = import_gtf_dict_to_massive_dict(gtf_dict)

	for line in reader:
		#print '----'
		#print line

		out_headder = ">" + line[0] + '\n'
		
		seq_details = all_anno_dict[line[0]].split(',')

		#print seq_details

		if seq_details[3] == '-':
			out_seq = genome_dict[seq_details[0]][int(seq_details[1])-1:int(seq_details[2])]
			#out_seq = reverse_compliment(out_seq)
		else:
			out_seq = genome_dict[seq_details[0]][int(seq_details[1])-1:int(seq_details[2])]


		out_transcriptome.write(out_headder)
		out_transcriptome.write(out_seq)
		out_transcriptome.write('\n')
		out_transcriptome.write('\n')


def extract_anno_pan_genome_csv(graph_obj, gtf_dict, out_file_name, refseq='', sim_threshold=1.0):

	isolate_list = gtf_dict.keys()



	added_list = []

	outfile_obj = open(out_file_name + 'anno.csv', 'w')
	
	csv_headder = 'gene'

	# Create headder for csv file

	for iso in isolate_list:
		csv_headder = csv_headder + ',' + iso

	csv_headder = csv_headder + '\n'

	outfile_obj.write(csv_headder)

	for isolate in isolate_list:
		gtf_lol = input_parser(gtf_dict[isolate], parse_as='gtf')
		timer = 0

		# For each gene for this isolate, see which other isolates have the same sequence

		for entry in gtf_lol:
			
			# For this gene for this isolate
			#print entry

			if entry[2] == 'exon' or entry[2] == 'gene':

				# Entries that are genes	
				
				found_in_list = check_isolates_in_region(graph_obj, entry[3], entry[4], isolate, threshold=sim_threshold, return_dict=False)
				
				if abs(int(entry[4])) < abs(int(entry[3])):
					print entry
				print '****'
				# this gene, is also found in these isolates
				print found_in_list

				if len(list(set(found_in_list) & set(added_list))) < 1:
					line_str = csv_headder

					#print entry[8].keys()

					if 'ID' in entry[8].keys():
						curr_gene = entry[8]['ID']

					if 'gene_id' in entry[8].keys():
						curr_gene = entry[8]['gene_id']
					print curr_gene

					line_str = line_str.replace('gene',curr_gene)
					line_str = line_str.replace(isolate,curr_gene)
					
					#print 'here we get the other iso annotations'

					for found_iso in found_in_list:

						found_iso = str(found_iso)

						#print '------'
						#print entry
						#print isolate
						#print str(found_iso)

						#print convert_coordinate(graph_obj, 100029, isolate, 'CDC1551')

						#new_coord_dict = convert_coordinates(graph_obj, entry[3], entry[4], isolate, str(found_iso))

						left_pos = convert_coordinate(graph_obj, entry[3], isolate, str(found_iso))

						right_pos = convert_coordinate(graph_obj, entry[4], isolate, str(found_iso))

						#print 'the pos list'
						print 'new pos'
						print left_pos, right_pos
						#print 'old pos'
						#print entry[3], entry[4]

						iso_gtf_lol = input_parser(gtf_dict[found_iso])

						#print iso_gtf_lol

						#print found_iso

						if left_pos != 'pos not found' and right_pos != 'pos not found':
							homo_gene = get_anno_from_coordinates(iso_gtf_lol, left_pos[str(found_iso)], right_pos[str(found_iso)], 10)

							print 'gene found!!'
							print homo_gene

							line_str = line_str.replace(found_iso, homo_gene)

							print line_str
						else:
							line_str = line_str.replace(found_iso, 'partial')


					for remaining_iso in isolate_list:
						line_str = line_str.replace(remaining_iso, '0')



					timer += 1
					#print timer
					outfile_obj.write(line_str)


		added_list.append(isolate)

def extract_unique_sequences(pangenome_csv):
	'''Extract the genes unique to each isolate, return a transcript file'''

	in_file = open(pangenome_csv, 'r')

	out_file_name = pangenome_csv[:-4] + '_unique.csv'

	out_trans_file = open(out_file_name, 'w')

	reader = csv.reader(in_file)

	rownum = 0

	for row in reader:
		# Save header row.
		if rownum == 0:
			header = row
			new_head = ''
			for col in header:
				new_head = new_head + col + ','

			new_head = new_head[:-1] + '\n'
			out_trans_file.write(new_head)

		else:
			col_total = 0
			for col in row[1:]:
				col_total = col_total + int(col)


			if col_total == 1:
				#print 'unique'
				#print row
				out_line = ",".join(row)
				out_line = out_line + '\n'
				out_trans_file.write(out_line)

		rownum += 1

	in_file.close()


def get_anno_from_coordinates(in_gtf_lol, start_pos, stop_pos, tollerence):
	
	if int(start_pos) > int(stop_pos):
		temp_start = stop_pos
		temp_stop = start_pos
		start_pos = temp_start
		stop_pos = temp_stop


	for anno in in_gtf_lol:
		if anno[2] == 'exon':
			if abs(int(anno[3]) - int(start_pos)) <= int(tollerence) and abs(int(anno[4]) - int(stop_pos)) <= int(tollerence):
				return anno[8]['gene_id']

		# if using gff3
		if anno[2] == 'gene':
			if abs(int(anno[3]) - int(start_pos)) <= int(tollerence) and abs(int(anno[4]) - int(stop_pos)) <= int(tollerence):
				return anno[8]['ID']

	return '1'

def get_gene_homo_gff(graph_obj, gtf_file, reference_name):

	gff_lod = input_parser(gtf_file)

	for anno in gff_lod:
		start_pos = anno[3]
		stop_pos = anno[4]

		found_in = check_isolates_in_region(graph_obj, start_pos, stop_pos, reference_name)

		print anno[8], found_in


# --------------------- Sequence extraction related

def extract_seq(graph_obj, seq_name):
	
	extracted_seq = ''

	pos_lol = []

	test_count = 0

	for node,data in graph_obj.nodes_iter(data=True):


		left_end_name = seq_name + '_leftend'
		right_end_name = seq_name + '_rightend'

		if seq_name in data['present_in']:
			#print 'right node'

			new_left_pos = data[left_end_name]
			new_right_pos = data[right_end_name]

			if int(data[left_end_name]) < 0:
				new_left_pos = abs(int(data[left_end_name]))
				new_right_pos = abs(int(data[right_end_name]))
				


			if int(new_left_pos) != 0 and int(new_right_pos) != 0 or node == seq_name + '_start':
				pos_lol.append([int(new_left_pos), int(new_right_pos), node])


	#print pos_lol

	pos_lol = sorted(pos_lol)

	#for aposthing in pos_lol:
	#	print aposthing
	#print len(pos_lol)
	#print 'oooo'
	#print test_count

	# REMOVE AND PUT AS TEST

	last_node_line = [0,0,0]

	for segment in pos_lol:
		#print '\n'
		#print segment
		#print graph_obj.node[str(segment[2])]['sequence']

		if segment[0] == last_node_line[1]:
			print '-'
			print segment

		if segment[1] == last_node_line[0]:
			print '+'
			print segment

		if segment[2].split('_')[1] != 'start':

			extracted_seq = extracted_seq + graph_obj.node[str(segment[2])]['sequence']

		last_node_line = segment

	return extracted_seq

def extract_original_seq(graph_obj, seq_name):
	
	extracted_seq = ''
	pos_lol = []

	test_count = 0

	for node,data in graph_obj.nodes_iter(data=True):

		is_negative = False

		left_end_name = seq_name + '_leftend'
		right_end_name = seq_name + '_rightend'

		if seq_name in data['present_in']:
			#print 'right node'

			new_left_pos = data[left_end_name]
			new_right_pos = data[right_end_name]

			if int(data[left_end_name]) < 0:
				is_negative = True
				new_left_pos = abs(int(data[left_end_name]))
				new_right_pos = abs(int(data[right_end_name]))
				


			if int(new_left_pos) != 0 and int(new_right_pos) != 0 or node == seq_name + '_start':
				if node != seq_name + '_stop':
					pos_lol.append([int(new_left_pos), int(new_right_pos), node, is_negative])


	#print pos_lol

	pos_lol = sorted(pos_lol)

	#for aposthing in pos_lol:
	#	print aposthing
	#print len(pos_lol)
	#print 'oooo'
	#print test_count

	# REMOVE AND PUT AS TEST

	last_node_line = [0,0,0]

	for segment in pos_lol:
		#print '\n'
		#print segment
		#print graph_obj.node[str(segment[2])]['sequence']

		if segment[0] == last_node_line[1]:
			print '-'
			print segment

		if segment[1] == last_node_line[0]:
			print '+'
			print segment

		if segment[2].split('_')[-1] != 'start':

			if segment[3] == True:
				# Node is negative
				#print reverse_compliment(graph_obj.node[str(segment[2])]['sequence'])
				extracted_seq = extracted_seq + reverse_compliment(graph_obj.node[str(segment[2])]['sequence'])
			else:
				extracted_seq = extracted_seq + graph_obj.node[str(segment[2])]['sequence']

		last_node_line = segment

	# Deal with the end node:
	'''
	stop_seq = graph_obj.node[seq_name + '_stop']['sequence']
	print stop_seq
	if len(stop_seq) > 0:
		print 'stop seq added'
		extracted_seq = extracted_seq + stop_seq
	'''

	return extracted_seq

def extract_seq_region(graph_obj, region_start, region_stop, seq_name):
	seq_string = extract_seq(graph_obj, seq_name)

	region_start = int(region_start) - 1
	region_stop = int(region_stop)

	return seq_string[region_start:region_stop]

def extract_original_seq_region(graph_obj, region_start, region_stop, seq_name):
	
	seq_string = extract_original_seq(graph_obj, seq_name)

	region_start = int(region_start) - 1
	region_stop = int(region_stop)

	return seq_string[region_start:region_stop]

def extract_iso_subgraph(graph_obj, isolate):
	iso_graph = nx.MultiDiGraph()
	iso_graph.graph['isolate'] = 'isolate'

	for node,data in graph_obj.nodes_iter(data=True):
		if isolate in data['present_in']:
			#print node
			iso_graph.add_node(node,data)

	iso_graph = link_nodes(iso_graph, isolate)

	return iso_graph

# For heaviest path function

def get_neighbour_most_iso(list_of_nodes, graph_obj):

	longest_list_node = 'nope'
	longest_list_length = 0

	for neigh_node in list_of_nodes:

		if len(graph_obj.node[neigh_node]['present_in'].split(',')) > longest_list_length:
			longest_list_node = neigh_node
			longest_list_length = len(graph_obj.node[neigh_node]['present_in'].split(','))

	return longest_list_node

def extract_heaviest_path(graph_obj, phyloinfo):

	# setup
	out_graph = nx.MultiDiGraph()

	# Determine starting node

	start_node = phyloinfo + '_start'

	start_node = graph_obj.successors(start_node)[0]

	# Adding, and let's go!

	out_graph.graph['start_node'] = start_node

	node_list = [start_node]

	curr_node = start_node
	
	while curr_node != 'nope':

		#print graph_obj.node[curr_node]

		neighbors_out = graph_obj.successors(curr_node)

		new_node = get_neighbour_most_iso(neighbors_out, graph_obj)

		node_list.append(new_node)

		curr_node = new_node

	node_list = node_list[:-1]

	for heavy_node in node_list:
		out_graph.add_node(heavy_node, graph_obj.node[heavy_node])


	# linking nodes
	head = 0
	tail = 1
	while tail < len(node_list):
		out_graph.add_edge(node_list[head], node_list[tail])
		head += 1
		tail += 1

	
	return out_graph

def extract_seq_heavy(graph_obj):
	
	start_node = graph_obj.graph['start_node']

	extracted_seq = graph_obj.node[start_node]['sequence']

	curr_node = start_node

	neighbors_list = graph_obj.successors(start_node)

	while len(neighbors_list) > 0:

		neighbors_list = graph_obj.successors(neighbors_list[0])

		if len(neighbors_list) > 0:
			if 'sequence' in graph_obj.node[neighbors_list[0]].keys():
				extracted_seq = extracted_seq + graph_obj.node[neighbors_list[0]]['sequence']


	return extracted_seq


def levenshtein(s1, s2):
	# This code was obtained from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

# ---------------------------------------------------- # Testing functions


def seq_addition_test(in_graph, node_ID, seq_fasta_paths_dict):

	# Get seq from node
	pass_seq_test = True

	the_node_seq = in_graph.node[node_ID]['sequence']
	print '\n'
	print 'node seq'
	print the_node_seq
	print '\n'

	node_seq_iso_list = in_graph.node[node_ID]['present_in'].split(',')

	for seq_isolate in node_seq_iso_list:
		
		reverseSeq = False
		complimentSeq = False

		iso_genome_seq = input_parser(seq_fasta_paths_dict[seq_isolate])

		start_pos = int(in_graph.node[node_ID][seq_isolate + '_leftend'])
		end_pos = int(in_graph.node[node_ID][seq_isolate + '_rightend'])


		if int(end_pos) < 0:
			end_pos = abs(end_pos)
			start_pos = abs(start_pos)
			reverseSeq = True
			#print "REVERSE"

		if start_pos > end_pos:
			tempstart = end_pos
			tempstop = start_pos
			end_pos = tempstop
			start_pos = tempstart
			complimentSeq = True
			#print 'COMP'


		sub_iso_seq = iso_genome_seq[0]['DNA_seq'][start_pos - 1:end_pos]

		if reverseSeq == True:
			sub_iso_seq = reverse_compliment(sub_iso_seq)

		if sub_iso_seq.upper() == the_node_seq.upper():
			1 == 1
			print seq_isolate + ' - Passed'
			print sub_iso_seq
		else:
			pass_seq_test = False
			print seq_isolate + ' - Failed ' + str(start_pos) + ' - ' + str(end_pos)
			print sub_iso_seq


	return pass_seq_test

def generate_graph_report(in_graph, out_file_name):
	nx_summary = nx.info(in_graph)

	report_file = open(out_file_name + '_report.txt', 'w')

	for line in nx_summary:
		report_file.write(line)

	report_file.write("\nIsolates: " + str(in_graph.graph['isolates']))

	# Length of sequence in the graph

	len_dict = {}

	for an_iso in in_graph.graph['isolates'].split(','):
		len_dict[an_iso] = 0

	comp_len = 0

	for n,d in in_graph.nodes_iter(data=True):
		if 'sequence' in d.keys():
			comp_len = comp_len + len(d['sequence'])

			for node_iso in d['present_in'].split(','):
				len_dict[node_iso] = len_dict[node_iso] + len(d['sequence'])


	#print comp_len

	report_file.write("\nCompressed sequence length: " + str(comp_len))

	#print len_dict

	for iso_name in len_dict.keys():
		report_file.write('\n' + iso_name + " length: " + str(len_dict[iso_name]))

	#print "Density: " + str(nx.density(in_graph))

	report_file.write("\nDensity: " + str(nx.density(in_graph)))

	return nx_summary

# ---------------------------------------------------- # Main loop of the program

if __name__ == '__main__':

	parser=argparse.ArgumentParser(
	description='''Welcome to GenGraph v0.1''',
	
	epilog="""Insanity is trying the same thing over and over and expecting different results""")

	parser.add_argument('toolkit', type=str,  default='test_mode', help='Select what you would like to do')

	parser.add_argument('--out_file_name', type=str, help='Prefix of the created file')

	parser.add_argument('--alignment_file', nargs=1,  help='The path to the alignment file')

	parser.add_argument('--backbone_file', nargs=1, default='default', help='The path to the backbone file')

	parser.add_argument('--out_format', nargs=1, default='default', help='Format for the output')

	parser.add_argument('--seq_file', type=str, help='Tab delimited text file with paths to the aligned sequences')

	parser.add_argument('--no_seq', dest='should_add_seq', action='store_false',  help='Create a graph genome with no sequence stored within')

	parser.add_argument('--make_circular', type=str, default='No', help='To circularise the graph for a sequence, give the name of that sequence')

	parser.add_argument('--recreate_check', dest='rec_check', action='store_true',  help='Set to True to attempt to recreate the input sequences from the graph and comapre to the origionals')

	parser.add_argument('--extract_sequence', type=str,  default='some_isolate', help='To circularise the graph for a sequence, give the name of that sequence')

	parser.add_argument('--isolate', type=str,  default='some_isolate', help='pass the isolate variable')

	parser.add_argument('--extract_sequence_range', nargs=2,  default=['all','all'], help='Extract sequence between two positions')

	parser.add_argument('--graph_file', type=str, help='Give the path to the graph file')

	parser.add_argument('--max_node_length', type=int,  default=-1, help='Max sequence length that can be aligned')

	parser.add_argument('--input_file', type=str, help='Generic input')

	parser.set_defaults(should_add_seq=True)

	parser.set_defaults(rec_check=False)


	args=parser.parse_args()


	print 'Current settings'
	print args


	if args.toolkit == 'test_mode':

		print "Place code here for quick testing purposes"


	if args.toolkit == 'test':
		print "Testing cuttent installation"
		curr_graph = './test_files/test_g.gml'
		#test_gtf_dict = {'H37Rv':h37rv_gtf, 'H37Ra':hra_gtf, 'CDC1551':cdc1551_gtf, 'F11':f11_gtf, 'c':c_gtf, 'W-148':w_gtf}
		test_input_file = './test_files/input_test.txt'
		test_backbone = './test_files/aln_test.backbone'


	if args.toolkit == 'make_genome_graph':

		# Requires:
		# --out_file_name
		# --seq_file

		# optional:
		# --recreate_check
		# --add_seq


		start_time = time.time()

		parsed_input_dict = parse_seq_file(args.seq_file)



		# --------------------------------------------------------------------------------- Initial global alignment

		if global_aligner == 'progressiveMauve' and args.backbone_file == 'default':
			
			print 'Conducting progressiveMauve'

			progressiveMauve_alignment(parsed_input_dict[2], args.out_file_name)
		
		# --------------------------------------------------------------------------------- Conversion to graph


		if args.backbone_file == 'default':
			bbone_file = 'globalAlignment_' + args.out_file_name + '.backbone'
		else:
			bbone_file = args.backbone_file[0]


		genome_aln_graph = create_genome_alignment_graph(bbone_file, parsed_input_dict[0], parsed_input_dict[1])

		nx.write_graphml(genome_aln_graph, 'intermediate_Graph.xml')

		# --------------------------------------------------------------------------------- node splitting

		if args.max_node_length != -1:

			genome_aln_graph = split_all_long_nodes(genome_aln_graph, args.max_node_length)

			print 'Nodes split'

			nx.write_graphml(genome_aln_graph, 'intermediate_split_Graph.xml')

	

		# --------------------------------------------------------------------------------- Local node realignment

		print 'Conducting local node realignment'

		genome_aln_graph = create_ungapped_graph(genome_aln_graph, parsed_input_dict[1])
		
		print 'Genome graph created'


		# --------------------------------------------------------------------------------- adding sequence to the graph

		if args.should_add_seq == True:

			if args.isolate[0] == 'some_isolate':
				ref_isolate = graph_obj.graph['isolates'].split(',')[0]

			print 'Sequence dict used:'
			print parsed_input_dict[1]

			genome_aln_graph = add_sequences_to_graph(genome_aln_graph, parsed_input_dict[1])
			print 'Sequence added'

		if args.make_circular != 'No':
			make_circular(genome_aln_graph, args.make_circular)
			print 'Graph circularised'


		if args.rec_check == True:
			for isolate in parsed_input_dict[1].keys():
				extracted_seq = extract_original_seq(genome_aln_graph, isolate)
				original_seq_from_fasta = input_parser(parsed_input_dict[1][isolate])

				if extracted_seq.upper() == original_seq_from_fasta[0]['DNA_seq'].upper():
					print 'Sequence recreate pass'
					recreate_check_result = 'Pass'
				else:
					print 'Sequence recreate fail'
					print len(extracted_seq)
					print len(original_seq_from_fasta[0]['DNA_seq'])
					print extracted_seq[-10:]
					print original_seq_from_fasta[0]['DNA_seq'][-10:]
					print '\n'
					print extracted_seq[:10]
					print original_seq_from_fasta[0]['DNA_seq'][:10]
					recreate_check_result = 'Fail'


		# Saving output

		if args.out_format == 'default':
			print 'Writing to GraphML'
			out_filename_created = args.out_file_name + '.xml'
			nx.write_graphml(genome_aln_graph, out_filename_created)

		if args.out_format[0] == 'serialize':
			print 'Writing to serialized file'
			pickle.dump(genome_aln_graph, open(args.out_file_name + '.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

		end_time = (time.time() - start_time)

		print "run time: " + str(end_time)

		generate_graph_report(genome_aln_graph, args.out_file_name)



	if args.toolkit == 'make_graph_from_fasta':

		# Requires:
		# --input_file
		# --out_file_name


		print 'In testing'

		fasta_object = input_parser(args.input_file)

		print 'Adding', len(fasta_object), 'sequences'

		seqStartDict = {}

		for seq_entry in fasta_object:
			seqStartDict[seq_entry['gene_details']] = 1

		

		new_graph = fasta_alignment_to_subnet(args.input_file, true_start=seqStartDict)

		nx.write_graphml(new_graph, 'intermediate_virus_Graph.xml')

		new_graph = add_sequences_to_graph(new_graph, fasta_object)

		nx.write_graphml(new_graph, args.out_file_name)
		


	if args.toolkit == 'region_alignment_score':

		print args.graph_file
		graph_obj = nx.read_graphml(args.graph_file)

		print 'this is the region alignment section'
		result = check_isolates_in_region(graph_obj, args.extract_sequence_range[0], args.extract_sequence_range[1], args.isolate, threshold=1.0, return_dict=True)

		print 'result is'
		print result



	if args.toolkit == 'extract_fasta_file':

		out_fasta = open(args.out_file_name, 'w')

		graph_obj = nx.read_graphml(args.graph_file)

		extracted_seq = extract_original_seq(graph_obj, args.isolate)

		fasta_headder = '>' + args.isolate + '\n'

		out_fasta.write(fasta_headder)

		n = 70
		for seq_line in [extracted_seq[i:i+n] for i in range(0, len(extracted_seq), n)]:

			out_fasta.write(seq_line + '\n')

		out_fasta.close()

	if args.toolkit == 'extract_region':
		imported_genome = nx.read_graphml(args.graph_file)

		print "Extracting from " + str(args.extract_sequence_range[0]) + " to " + str(args.extract_sequence_range[1])

		print extract_original_seq_region(imported_genome, args.extract_sequence_range[0], args.extract_sequence_range[1], args.isolate)




	if args.toolkit == 'extract_pan_transcriptome':

		# This needs to deal with the required GTF dict and seq_file_dict

		# Needs
		# --graph_file
		# --seq_file
		# --out_file_name

		print 'Extracting...'

		parsed_input_dict = parse_seq_file(args.seq_file)

		graph_obj = nx.read_graphml(args.graph_file)

		test_gtf_dict = parsed_input_dict[3]


		if args.isolate == 'some_isolate':
			ref_isolate = ''
		else:
			ref_isolate = args.isolate

		# Extracting the csv from the graph
		print 'Extracting annotated pan genome csv'
		extract_anno_pan_genome_csv(graph_obj, test_gtf_dict, args.out_file_name, sim_threshold=0.95)
		print 'Extracting pan genome csv'
		extract_pan_genome_csv(graph_obj, test_gtf_dict, args.out_file_name, hom_threshold=0.95, refseq=ref_isolate)
		print 'Extracting pan genome transcriptome'
		create_fasta_from_pangenome_csv(args.out_file_name + '.csv', test_gtf_dict, parsed_input_dict, args.out_file_name)

		# Converting the csv to a fasta file of transcripts
		#create_fasta_from_pangenome_csv(args.input_file, test_gtf_dict, parsed_input_dict, args.out_file_name)

	if args.toolkit == 'generate_report':
		imported_genome = nx.read_graphml(args.graph_file)

		print generate_graph_report(imported_genome, args.out_file_name)














