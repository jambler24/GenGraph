# Dependencies

import networkx as nx


# Built in

from subprocess import call

from operator import itemgetter

import argparse

import sys, os

import csv

import copy

import time

try:
	import cPickle as pickle
except:
	print('cPickle not found')
	import pickle

import logging



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

# ---------------------------------------------------- New gg object class


class GgDiGraph(nx.DiGraph):

	def get_sequence(self, region_start, region_stop, seq_name):
		'''
		Returns the sequence found between the given coordinates for a strain.
		:param region_start: The start coordinate
		:param region_stop: The stop coordinate
		:param seq_name: The ID for the sequence that the coordinates refer to. Often the isolate.
		:return: A sequence string
		'''

		#seq_string = extract_original_seq_region_fast(self, region_start, region_stop, seq_name)
		# Change once fast version is fixed. Currently it doesn't handle inv sequence
		seq_string = extract_original_seq_region(self, region_start, region_stop, seq_name)

		return seq_string

	def transfer_annotations(self, annofile, source_seq_id, target_seq_id, out_file_name, chr_name=''):
		'''
		Given a annotation file (GTF, GFF) for a sequence, convert the coordinates relative to another sequence.
		:param graph: Genome graph object
		:param annofile: The filepath to the annotation file
		:param source_seq_id:
		:param target_seq_id:
		:param out_file_name:
		:return:
		'''

		# parse the annotation file

		anno_dict = input_parser(annofile)

		out_file = open(out_file_name, 'w')
		out_file.write('##gff-version 3\n')

		if len(chr_name) == 0:
			chr_name = target_seq_id

		for line in anno_dict:

			conv_coords = convert_coordinates(self, line[3], line[4], source_seq_id, target_seq_id)

			if conv_coords is not None:

				line[0] = chr_name
				line[3] = str(conv_coords[target_seq_id + '_leftend'])
				line[4] = str(conv_coords[target_seq_id + '_rightend'])
				info_dict = line[8]
				line[8] = 'ID=' + info_dict['ID'] + ';' + 'Name=' + info_dict['Name']

				line_string = '\t'.join(line)
				out_file.write(line_string)

		return 'Done'


# ---------------------------------------------------- New functions under testing


def import_gg_graph(path):
	'''
	Import a GG graph
	:param path: file path
	:return: a GG graph genome object
	'''
	graph_obj = nx.read_graphml(path)

	out_graph = GgDiGraph(graph_obj)

	return out_graph


def get_neighbours_context(graph, source_node, label, dir='out'):
	'''Retrieves neighbours of a node using only edges containing a certain label'''
	'''Created by shandu.mulaudzi@gmail.com'''

	# Changes seequence to ids
	list_of_neighbours = []
	in_edges = graph.in_edges(source_node, data='ids')
	out_edges = graph.out_edges(source_node, data='ids')
	for (a,b,c) in in_edges:
		if label in c['ids'].strip(','):
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
				vcf_header_line = line.split('\t')
				vcf_header_line[0] = vcf_header_line[0][1:]
				vcf_header_line[-1] = vcf_header_line[-1].strip()

			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {vcf_header_line[0]:entries[0], vcf_header_line[1]:entries[1], vcf_header_line[2]:entries[2], vcf_header_line[3]:entries[3], vcf_header_line[4]:entries[4], vcf_header_line[5]:entries[5], vcf_header_line[6]:entries[6], vcf_header_line[7]:entries[7], vcf_header_line[8]:entries[8], vcf_header_line[9]:entries[9].strip(), 'ORIGIN':entry_label} 

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
							#print info_byte
							info_byte = info_byte.split('=')
							extra_info_dict[info_byte[0]] = info_byte[1]


						if 'locus_tag' in extra_info_dict.keys():

							gtf_info_dict['gene_id'] = extra_info_dict['locus_tag']
							gtf_info_dict['locus_tag'] = extra_info_dict['locus_tag']

							entries[8] = gtf_info_dict
					
							list_of_lists.append(entries)

							#quit()



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


def reshape_fastaObj(in_obj):
	out_fastaObj = {}

	for seq_ent in in_obj:
		out_fastaObj[seq_ent['gene_details']] = seq_ent['DNA_seq']

	return out_fastaObj


def parse_seq_file(path_to_seq_file):

	seq_file_dict = input_parser(path_to_seq_file)

	A_seq_label_dict = {}
	A_input_path_dict = {}
	ordered_paths_list = []
	anno_path_dict = {}

	for a_seq_file in seq_file_dict:
		logging.info(a_seq_file)
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
	'''
	Returns the reverse compliment of a sequence
	:param seq_in: A string representing a DNA sequence
	:return: The reverse-compliment of that sequence
	'''
	comp_seq = compliment_DNA(seq_in)
	rev_comp_seq = reverse(comp_seq)
	return rev_comp_seq


def compliment_Base(base_in):
	'''
	Returns the compliment of the given base. Assuming DNA. Non-DNA bases are returned as a empty string
	:param base_in: A single letter string.
	:return: The compliment of that base. So A -> T
	'''
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
	'''
	Returns the length of a node in base pairs
	:param node_name:
	:return:
	'''
	iso_string = node_name['ids'].split(',')[0]
	left_pos = node_name[iso_string + '_leftend']
	right_pos = node_name[iso_string + '_rightend']

	return bp_distance(left_pos, right_pos)


def export_to_fasta(sequence, header, filename):
	file_obj = open(filename + '.fa', 'w')

	header_line = '>' + header + '\n'
	file_obj.write(header_line)

	n=80
	seq_split = [sequence[i:i+n] for i in range(0, len(sequence), n)]

	for line in seq_split:
		line = line + '\n'
		file_obj.write(line)
	file_obj.close()

# ---------------------------------------------------- Graph generating functions


def add_missing_nodes(a_graph, input_dict):

	from operator import itemgetter

	logging.info(input_dict[1].keys())

	iso_list = a_graph.graph['isolates'].split(',')

	for isolate in iso_list:

		isolate_Seq = input_parser(input_dict[1][isolate])
		isolate_Seq = isolate_Seq[0]['DNA_seq']

		isolate_node_list = []
		for node, data in a_graph.nodes_iter(data=True):
			if isolate in data['ids'].split(','):
				isolate_node_list.append(node)

		presorted_list = []
		for a_node in isolate_node_list:
			presorted_list.append((a_node, abs(a_graph.node[a_node][isolate + '_leftend']), abs(a_graph.node[a_node][isolate + '_rightend'])))
			if abs(a_graph.node[a_node][isolate + '_leftend']) > abs(a_graph.node[a_node][isolate + '_rightend']):
				logging.warning('problem node' + str(a_node))
				logging.warning(a_graph.node[a_node][isolate + '_leftend'])

		sorted_list = sorted(presorted_list, key=itemgetter(1))

		count = 0

		while count < len(sorted_list) - 1:

			if sorted_list[count][2] != sorted_list[count + 1][1] - 1:
				#print 'gap', sorted_list[count], sorted_list[count + 1]
				new_node_dict = {isolate + '_leftend':sorted_list[count][2] + 1, isolate + '_rightend':sorted_list[count + 1][1] - 1, 'ids':isolate}
				#print new_node_dict
				new_node_dict['sequence'] = isolate_Seq[sorted_list[count][2]:sorted_list[count + 1][1] - 1]

				node_ID = isolate + "_" + str(count)
				#print node_ID
				a_graph.add_node(node_ID, new_node_dict)
				#a_graph.add_node(node_ID)
				#nx.set_node_attributes(a_graph, {node_ID: new_node_dict})

			count+=1


def node_check(a_graph):
	from operator import itemgetter
	logging.info('checking nodes')
	doespass = True

	iso_list = a_graph.graph['isolates'].split(',')

	for isolate in iso_list:
		isolate_node_list = []
		for node, data in a_graph.nodes_iter(data=True):
			if isolate in data['ids'].split(','):
				isolate_node_list.append(node)

		#print isolate_node_list
		presorted_list = []
		for a_node in isolate_node_list:
			presorted_list.append((a_node, abs(int(a_graph.node[a_node][isolate + '_leftend'])), abs(int(a_graph.node[a_node][isolate + '_rightend']))))
			if abs(int(a_graph.node[a_node][isolate + '_leftend'])) > abs(int(a_graph.node[a_node][isolate + '_rightend'])):
				logging.warning('problem node', a_node)
				logging.warning(a_graph.node[a_node][isolate + '_leftend'])

		sorted_list = sorted(presorted_list,key=itemgetter(1))

		count = 0

		while count < len(sorted_list) - 1:

			if sorted_list[count][2] != sorted_list[count + 1][1] - 1:
				logging.error('error 1: gaps in graph')
				logging.error(isolate)
				logging.error(sorted_list[count])
				logging.error(sorted_list[count + 1])
				logging.error('Gap length:', sorted_list[count + 1][1] - sorted_list[count][2] - 1)
				doespass = False

			if sorted_list[count][2] >= sorted_list[count + 1][1]:
				logging.error('error 2')
				logging.error(isolate)
				logging.error(sorted_list[count])
				logging.error(sorted_list[count + 1])
				doespass = False

			if sorted_list[count][1] == sorted_list[count - 1][2] - 1:
				logging.error('error 3: last node end close to node start. If this is not a SNP, there is an error')
				logging.error(isolate)
				logging.error(sorted_list[count])
				logging.error(sorted_list[count - 1])
				doespass = False

			if sorted_list[count][1] > sorted_list[count][2]:
				logging.error('error 4: start greater than stop')
				logging.error(isolate)
				logging.error(sorted_list[count])
				logging.error(a_graph.node[sorted_list[count][0]])
				doespass = False


			count+=1
	return doespass


def refine_initGraph(a_graph):

	iso_list = a_graph.graph['isolates'].split(',')

	for isolate in iso_list:
		isolate_node_list = []
		for node, data in a_graph.nodes_iter(data=True):
			if isolate in data['ids'].split(','):
				isolate_node_list.append(node)

		presorted_list = []
		for a_node in isolate_node_list:
			presorted_list.append((a_node, abs(a_graph.node[a_node][isolate + '_leftend']), abs(a_graph.node[a_node][isolate + '_rightend'])))
			if abs(a_graph.node[a_node][isolate + '_leftend']) > abs(a_graph.node[a_node][isolate + '_rightend']):
				logging.warning('problem node' + str(a_node))
				logging.warning(a_graph.node[a_node][isolate + '_leftend'])

		sorted_list = sorted(presorted_list, key=itemgetter(1))

		count = 0

		while count < len(sorted_list) - 1:

			if sorted_list[count][2] == sorted_list[count + 1][1]:
				logging.warning('problem node')
				logging.warning(sorted_list[count])
				logging.warning(sorted_list[count + 1])
				logging.warning(a_graph.node[sorted_list[count][0]])
				for a_node_isolate in a_graph.node[sorted_list[count][0]]['ids'].split(','):
					if a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] < 0:
						a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] = a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] + 1
					if a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] > 0:
						a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] = a_graph.node[sorted_list[count][0]][a_node_isolate + '_rightend'] - 1

				logging.warning(a_graph.node[sorted_list[count][0]])

			count+=1


def bbone_to_initGraph(bbone_file, input_dict):

	"""
	This function takes a bbone file created by MAUVE that represents the colinear blocks of sequences
	and uses it to create the initial graph.
	:param bbone_file:
	:param input_dict:
	:return:
	"""

	# This is the network that will be returned in the end
	genome_network = nx.MultiDiGraph()

	all_iso_in_graph = ''
	all_iso_in_graph_list = []
	iso_length_dict = {}
	iso_largest_node = {}

	has_start_dict = {}
	has_stop_dict = {}

	# Get the identifiers of all the sequences in the bbone file from the input_dict.
	for isolate_name in input_dict[0].keys():
		all_iso_in_graph_list.append(input_dict[0][isolate_name])
		all_iso_in_graph = all_iso_in_graph + input_dict[0][isolate_name] + ','

	for iso in all_iso_in_graph_list:
		has_start_dict[iso + '_leftend'] = False
		has_stop_dict[iso] = False
		# The iso_largest_node keeps track of the current highest nucleotide position for an isolate. This is so when
		# the sequences are aligned and there is a bit at the end og the sequence with no match, it can be made into an
		# end node.
		iso_largest_node[iso] = 0

	for iso in all_iso_in_graph_list:
		iso_length = len(input_parser(input_dict[1][iso])[0]['DNA_seq'])
		iso_length_dict[iso] = iso_length

	all_iso_in_graph = all_iso_in_graph[:-1]

	genome_network.graph['isolates'] = all_iso_in_graph

	# Parse the BBone file
	backbone_lol = input_parser(bbone_file)

	# Get the first line, which identifies which column relates to what identifier.
	header_line = backbone_lol[0]

	new_header_line = []

	# Because MAUVE does not use teh sequence names in the bbone file, you need to replace the generic names like 'seq0'
	# with the identifiers like 'H37Rv'
	for header_item in header_line:
		for seqID in input_dict[0].keys():
			if seqID in header_item:
				new_header_line.append(header_item.replace(seqID, input_dict[0][seqID]))

	# This is now the list of lists without the header line.
	backbone_lol_headless = backbone_lol[1:]

	#print backbone_lol_headless[0]

	# When naming nodes, keep track of how many have been created so each name is unique.
	node_count = 1

	for line in backbone_lol_headless:
		# Organise the info into a dict
		node_dict = {}
		# This header_count links the sequence positions with the sequence identifiers.
		header_count = 0
		# For each isolate / sequence identifier in the list. The header_item is the sequence identifier (eg., H37Rv)
		for header_item in new_header_line:
			# if the value is 0 at this position, the sequence was not found for this block.
			if line[header_count] != '0':
				# This is the record of the sequences represented by this node and their relative nucleotide positions.
				node_dict[header_item] = int(line[header_count])

				# Check for start node
				if abs(int(line[header_count])) == 1:
					has_start_dict[header_item] = 'Aln_' + str(node_count)
				# Check for end node
				curr_iso = header_item.replace('_leftend', '')
				curr_iso = curr_iso.replace('_rightend', '')

				if abs(int(line[header_count])) == iso_length_dict[curr_iso]:
					has_stop_dict[curr_iso] = 'Aln_' + str(node_count)

				# Check if largest node, if it is larger that the last largest, make it the new one.

				if abs(int(line[header_count])) > iso_largest_node[curr_iso]:
					iso_largest_node[curr_iso] = abs(int(line[header_count]))

			header_count += 1

		found_in_list = []
		for item in node_dict.keys():
			item = item.split('_')[:-1]
			item = '_'.join(item)
			if item not in found_in_list:
				found_in_list.append(item)

		#print found_in_list
		found_in_string = ','.join(found_in_list)
		#print found_in_string
		node_dict['ids'] = found_in_string

		node_ID = 'Aln_' + str(node_count)

		node_dict['name'] = node_ID

		genome_network.add_node(node_ID, node_dict)
		#genome_network.add_node(node_ID)
		#nx.set_node_attributes(genome_network, {node_ID: node_dict})

		node_count += 1

	for an_iso in all_iso_in_graph_list:
		for a_largest_node in iso_largest_node.keys():
			an_LGN_iso = a_largest_node.replace('_leftend', '')
			an_LGN_iso = an_LGN_iso.replace('_leftend', '')
			if an_iso == an_LGN_iso:
				if iso_largest_node[a_largest_node] != iso_length_dict[an_iso]:
					#print an_iso
					node_ID = 'Aln_' + str(node_count)
					node_dict = {'ids': an_iso, an_iso + '_leftend': int(iso_largest_node[a_largest_node]) + 1, an_iso + '_rightend': int(iso_length_dict[an_iso])}
					#print node_dict
					genome_network.add_node(node_ID, node_dict)
					#genome_network.add_node(node_ID)
					#nx.set_node_attributes(genome_network, {node_ID: node_dict})

					has_stop_dict[an_iso] = node_ID

					node_count += 1

	# print(has_stop_dict)

	# Use coords to extract fasta files for alignment

	# Parse the alignment files and place in graph
	return genome_network


def realign_all_nodes(inGraph, input_dict):
	logging.info('Running realign_all_nodes')

	realign_node_list = []

	iso_list = inGraph.graph['isolates'].split(',')

	# Load genomes into memory

	# Only need to realign nodes with more than one isolate in them
	for node, data in inGraph.nodes_iter(data=True):
		# print(data)
		if len(data['ids'].split(',')) > 1:

			realign_node_list.append(node)

	# Realign the nodes. This is where multiprocessing will come in.
	for a_node in realign_node_list:

		inGraph = local_node_realign_new(inGraph, a_node, input_dict[1])

	nx.write_graphml(inGraph, 'intermediate_split_unlinked.xml')

	return inGraph


def link_nodes(graph_obj, sequence_name, node_prefix='gn'):
	"""
	Create edges between all the nodes from a given sequence. The sequence is generally a
	isolate, and is specified in the sequence file. If there is a gap in the sequence that is not
	represented by a node, a new node will be created.
	:param graph_obj: A networkx graph object containing nodes.
	:param sequence_name: The identifier used for the sequence.
	:param node_prefix: Depreciated, will be removed.
	:return: A networkx graph object.
	"""

	# This list is a list of all the existing nodes start and stop positions with the node name formatted at start, stop, node_name

	pos_lol = []

	for node, data in graph_obj.nodes_iter(data=True):

		left_end_name = sequence_name + '_leftend'
		right_end_name = sequence_name + '_rightend'

		if sequence_name in data['ids'].split(','):

			new_left_pos = abs(int(data[left_end_name]))
			new_right_pos = abs(int(data[right_end_name]))

			if new_left_pos > new_right_pos:
				print('Doing the switch!')
				temp_left = new_right_pos
				temp_right = new_left_pos
				new_left_pos = temp_left
				new_right_pos = temp_right

			if int(new_left_pos) != 0 and int(new_right_pos) != 0 or node == sequence_name + '_start':
				pos_lol.append([new_left_pos, new_right_pos, node])

	count = 0


	#print 'Sorting'

	all_node_list = sorted(pos_lol, key=lambda start_val: start_val[0])


	# -------------------------------------------------------- Here we add edges to the graph, weaving in the new nodes

	count = 0

	#print 'Add new edges'

	edges_obj = graph_obj.edges()

	while count < (len(all_node_list) - 1):

		node_1 = all_node_list[count][2]
		node_2 = all_node_list[count+1][2]

		if (node_1, node_2) in edges_obj:

			if sequence_name not in graph_obj.edge[node_1][node_2][0]['ids'].split(','):


				new_seq_list = graph_obj.edge[node_1][node_2][0]['ids'] + ',' + sequence_name
				#print new_seq_list

			else:
				print('Error found.')
				quit()
				new_seq_list = graph_obj.edge[node_1][node_2]['ids']
			#print new_seq_list

			graph_obj.edge[node_1][node_2][0]['ids'] = new_seq_list

			#print 'appending'
		else:
			#print 'new edge'
			#print [(node_1, node_2, dict(ids=sequence_name))]

			graph_obj.add_edge(node_1, node_2, ids=sequence_name)

		count += 1

	logging.info('Edges added')

	return graph_obj


def link_all_nodes(graph_obj):
	logging.info('Running link_all_nodes')
	for isolate in graph_obj.graph['isolates'].split(','):
		logging.info(isolate)
		graph_obj = link_nodes(graph_obj, isolate)
	return graph_obj


def add_sequences_to_graph(graph_obj, paths_dict):

	logging.info('Adding sequences')

	for node, data in graph_obj.nodes_iter(data=True):

		seq_source = data['ids'].split(',')[0]
		is_reversed = False
		is_comp = False

		if len(seq_source) < 1:
			logging.error('No ids current node')
			logging.error(node)

		else:
			ref_seq = input_parser(paths_dict[1][seq_source])[0]['DNA_seq']

			# Check orientation
			if int(data[seq_source + '_leftend']) < 0:
				is_reversed = True

			seq_start = abs(int(data[seq_source + '_leftend']))
			seq_end = abs(int(data[seq_source + '_rightend']))

			if seq_start > seq_end:
				#new_seq_start = seq_end
				#new_seq_end = seq_start
				#seq_end = new_seq_end
				#seq_start = new_seq_start
				logging.error('Something wrong with orientation')

			'''
			if is_reversed != True:
				seq_start = seq_start - 1

			if is_reversed == True:
				seq_end = seq_end + 1
			'''
			seq_start = seq_start - 1

			node_seq = ref_seq[seq_start:seq_end].upper()

			if is_reversed:
				#print 'seq was rev comp' + node
				node_seq = reverse_compliment(node_seq)

			graph_obj.node[node]['sequence'] = node_seq

	return graph_obj


def add_sequences_to_graph_fastaObj(graph_obj, imported_fasta_object):

	logging.info('Adding sequences')

	seqObj = reshape_fastaObj(imported_fasta_object)

	for node, data in graph_obj.nodes_iter(data=True):

		seq_source = data['ids'].split(',')[0]
		is_reversed = False
		is_comp = False

		if len(seq_source) < 1:
			logging.error('No ids current node')
			logging.error(node)

		else:
			ref_seq = seqObj[seq_source]

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

			if is_reversed is True:
				seq_start = seq_start - 1

			if is_reversed:
				seq_start = seq_start
				seq_end = seq_end + 1

				logging.info(str(seq_start) + ' ' + str(seq_end))
				logging.info(ref_seq[seq_start:seq_end].upper())

			node_seq = ref_seq[seq_start:seq_end].upper()

			graph_obj.node[node]['sequence'] = node_seq

	return graph_obj


def make_circular(graph_obj, seq_name):

	start_node_name = seq_name + '_start'
	end_node_name = seq_name + '_stop'

	graph_obj.add_edges_from([(start_node_name, end_node_name, dict(ids=seq_name))])

	return graph_obj


def check_isolates_in_region(graph_obj, start_pos, stop_pos, reference_name, threshold=1.0, return_dict=False, simmilarity_measure='percentage'):
	'''Retrieve the nodes from a graph spanning a region'''


	print(start_pos, stop_pos)

	#print '\nchecking iso in region'

	graph_isolate_list = graph_obj.graph['isolates'].split(',')

	expected_ref_length = int(stop_pos) - int(start_pos)
	expected_ref_length = expected_ref_length + 1

	if expected_ref_length < 0:
		print('POSSIBLE ERROR, start larger than stop')
		print(int(stop_pos), int(start_pos))
		print(expected_ref_length)

	#print expected_ref_length
	#print reference_name

	# Extract subgraph
	node_leftend_label = reference_name + '_leftend'
	node_rightend_label = reference_name + '_rightend'

	# Identify the nodes that contain the start and stop positions

	#print int(start_pos)
	#print int(stop_pos)
	print(reference_name)

	for node, data in graph_obj.nodes_iter(data=True):
		if reference_name in data['ids'].split(','):

			if int(data[node_leftend_label]) > 0:
				#print 'positive'
				if abs(int(data[node_leftend_label])) <= int(start_pos) <= abs(int(data[node_rightend_label])):
					start_node = node
					print('Ping')
					#print data[node_leftend_label]
					print(start_pos)
					#print data[node_rightend_label]
					start_overlap =  bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) <= int(stop_pos) <= abs(int(data[node_rightend_label])) or abs(int(data[node_leftend_label])) >= int(stop_pos) >= abs(int(data[node_rightend_label])):
					print('ping')
					stop_node = node
					print(stop_pos)
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])
					#stop_overlap = int(stop_pos) - int(data[node_leftend_label]) + 1

			if int(data[node_leftend_label]) < 0:
				#print 'negative'
				#print data[node_leftend_label]
				#print data[node_rightend_label]
				if abs(int(data[node_leftend_label])) <= int(start_pos) <= abs(int(data[node_rightend_label])):
					start_node = node
					print('Ping')
					print('neg node')
					#print data[node_leftend_label]
					print(start_pos)
					#print data[node_rightend_label]
					start_overlap =  bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) <= int(stop_pos) <= abs(int(data[node_rightend_label])):
					logging.info('neg node')
					stop_node = node
					logging.info(stop_pos)
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])
					#stop_overlap = int(stop_pos) - int(data[node_leftend_label]) + 1




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
		if reference_name in graph_obj.node[path_node]['ids']:
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
			if a_isolate in graph_obj.node[a_node]['ids'] and reference_name in graph_obj.node[a_node]['ids']:
				iso_sim_dict[a_isolate] = iso_sim_dict[a_isolate] + bp_length_node(graph_obj.node[a_node])

			elif a_isolate in graph_obj.node[a_node]['ids'] and reference_name not in graph_obj.node[a_node]['ids']:
				iso_diff_dict[a_isolate] = iso_diff_dict[a_isolate] + bp_length_node(graph_obj.node[a_node])

	for a_iso in graph_isolate_list:
		if a_iso in graph_obj.node[start_node]['ids']:
			iso_sim_dict[a_iso] = iso_sim_dict[a_iso] - start_node_nonoverlap
		if a_iso in graph_obj.node[stop_node]['ids']:
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
	'''
	Convert coordinates from one sequence to another
	:param graph_obj:
	:param q_start:
	:param q_stop:
	:param ref_iso:
	:param query_iso:
	:return:
	'''

	# Find the right node

	node_leftend_label = ref_iso + '_leftend'
	node_rightend_label = ref_iso + '_rightend'

	query_leftend_label = query_iso + '_leftend'
	query_rightend_label = query_iso + '_rightend'

	ref_node_left = ''
	ref_node_right = ''

	query_node_left = ''
	query_node_right = ''

	for node, data in graph_obj.nodes_iter(data=True):
		if ref_iso in data['ids'].split(',') and query_iso in data['ids'].split(','):
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


	for node, data in graph_obj.nodes_iter(data=True):

		# Look at all nodes that contain both the isolates

		if ref_iso in data['ids'].split(',') and query_iso in data['ids'].split(','):

			# find the node

			if int(data[node_leftend_label]) > 0:
				#print 'positive node search'
				if abs(int(data[node_leftend_label])) <= abs(int(q_position)) <= abs(int(data[node_rightend_label])):
					logging.info('positive local found!!!!')
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
					logging.info('Negative local found')
					logging.info(q_position)
					logging.info(data[node_leftend_label], data[node_rightend_label])

					#print data

					ref_node_left = data[node_leftend_label]
					ref_node_right = data[node_rightend_label]
					query_node_left = data[query_leftend_label]
					query_node_right = data[query_rightend_label]

					logging.info(query_node_left, query_node_right)

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

					logging.info('conversion_factor')
					logging.info(conversion_factor)
					logging.info(query_high_num)
					new_r_start = query_high_num - conversion_factor

					logging.info(new_r_start)
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
			pos_deets = isolate + ',' + entry[3] + ',' + entry[4] + ',' + entry[6]
			if 'gene_id' not in entry[8].keys():
				logging.error('key error')
				logging.error(entry)
				logging.error(isolate)
				quit()
			else:
				gene_name = entry[8]['gene_id']
				all_dict[gene_name] = pos_deets

	#print all_dict
	return all_dict


def fasta_alignment_to_subnet(fasta_aln_file, true_start={}, node_prefix='X', orientation={}, re_link_nodes=True, add_seq=False):
	"""
	Conversion of a fasta alignment file from a MSA to a subnetwork that can replace a block node in the main network.
	:param fasta_aln_file:
	:param true_start:
	:param node_prefix:
	:param orientation:
	:param re_link_nodes:
	:param add_seq:
	:return:
	"""

	aln_lol = input_parser(fasta_aln_file)

	offset_dict = {}
	gap_dict = {}
	all_isolate_list = []
	# Get the length of each of the sequences without gaps
	seq_len_dict = {}

	for fasta_entry in aln_lol:

		strip_seq = fasta_entry['DNA_seq'].replace('-','')

		seq_len_dict[fasta_entry['gene_details']] = len(strip_seq)

		all_isolate_list.append(fasta_entry['gene_details'])
		fasta_entry['DNA_seq'] = fasta_entry['DNA_seq'].upper()

	# Setting the orientation dict for later
	if len(orientation) == 0:
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

		for a_seq in aln_lol:

			current_keys = block_dict.keys()
			if a_seq['DNA_seq'][count] in current_keys:
				block_dict[a_seq['DNA_seq'][count]].append(a_seq['gene_details'])
			else:
				block_dict[a_seq['DNA_seq'][count]] = [a_seq['gene_details']]

		block_list.append(block_dict)

		count += 1

	# Making sure the start pos is correct

	if len(true_start.keys()) > 0:
		for an_isolate in all_isolate_list:
			true_start[an_isolate] = true_start[an_isolate] - 1

	else:
		for an_isolate in all_isolate_list:
			true_start[an_isolate] = 0

	# Adding additional info and indexing

	last_bline = {'blocklist':[]}

	bpos_count = 1

	# Working along the block from here ----------------------
	# bline is block line
	for bline in block_list:

		for an_isolate in all_isolate_list:

			true_start[an_isolate] += 1

			for base in bline.keys():
				if base == '-' and an_isolate in bline[base]:

					true_start[an_isolate] = true_start[an_isolate] - 1


		bline['global_pos'] = bpos_count
		#bline['relative_pos'] = copy.deepcopy(true_start)
		bline['relative_pos'] = pickle.loads(pickle.dumps(true_start, -1))

		# See if this bline is a new block / node set
		block_list = []
		for base in bline.keys():
			if base != '-' and base != 'global_pos' and base != 'relative_pos':

				block_list.append(bline[base])

		block_list = sorted(block_list)

		bline['blocklist'] = block_list




		# Here we call blocks
		if bline['blocklist'] != last_bline['blocklist']:

			block_ends_list.append(last_bline)
			block_ends_list.append(bline)


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
		if is_even(count):
			start_block_list.append(end_block)
		else:
			end_block_list.append(end_block)

		count += 1


	# Now we pair the stop and start
	# Also begin adding nodes to the local_node_network
	count = 0
	node_count = 1
	while count < len(start_block_list):

		if start_block_list[count]['blocklist'] != end_block_list[count]['blocklist']:
			logging.error("ERROR IN BLOCK LIST PAIRS")

		for block_group in start_block_list[count]['blocklist']:
			# Each of these becomes a node
			new_node_name = node_prefix + '_' + str(node_count)
			block_group_string = ",".join(block_group)

			local_node_network.add_node(new_node_name, ids=block_group_string)

			# Adding start / stop for node
			# this is where we do the +/- thing.
			for block_isolate in block_group:
				if orientation[block_isolate] == '+':

					local_node_network.node[new_node_name][block_isolate + '_leftend'] = int(start_block_list[count]['relative_pos'][block_isolate])
					local_node_network.node[new_node_name][block_isolate + '_rightend'] = int(end_block_list[count]['relative_pos'][block_isolate])

				elif orientation[block_isolate] == '-':

					#So,
					# IF the positions are 1-10, 11-15, 15-40, what would the reverse positions be?
					# len node - pos - 1
					# So, we need the total length of the nodes, found in seq_len_dict

					local_node_network.node[new_node_name][block_isolate + '_rightend'] = -1 * int(int(seq_len_dict[block_isolate]) - start_block_list[count]['relative_pos'][block_isolate] - 1)
					local_node_network.node[new_node_name][block_isolate + '_leftend'] = -1 * int(int(seq_len_dict[block_isolate]) - end_block_list[count]['relative_pos'][block_isolate] - 1)

				else:
					logging.error("ORIENTATION MISSING")
					logging.error(orientation)
					quit()

			node_count += 1

		count += 1

	if re_link_nodes == True:
		for a_isolate in all_isolate_list:
			local_node_network = link_nodes(local_node_network, a_isolate, node_prefix='gn')

	# Here we add the seq if required

	if add_seq == True:
		new_fasta_list = []

		for a_seq in aln_lol:
			new_fasta_list.append({'DNA_seq':a_seq['DNA_seq'].replace('-', ''), 'gene_details':a_seq['gene_details']})

		local_node_network = add_sequences_to_graph_fastaObj(local_node_network, new_fasta_list)

	node_check(local_node_network)

	return local_node_network


def local_node_realign_new(in_graph, node_ID, seq_fasta_paths_dict):

	logging.info('Fast local node realign: ' + node_ID)
	logging.info(in_graph.node[node_ID])

	in_graph = nx.MultiDiGraph(in_graph)

	node_data_dict = in_graph.node[node_ID]

	# Make temp fasta file and record the start positions into a dict

	node_seq_start_pos = {}

	# Store the orientation of the sequences in the node to pass to the fasta to graph conversion function
	orientation_dict = {}

	temp_fasta_file = open('temp_unaligned.fasta', 'w')

	for node_isolate in node_data_dict['ids'].split(','):
		iso_full_seq = input_parser(seq_fasta_paths_dict[node_isolate])[0]['DNA_seq'].upper()

		#print node_isolate
		#print '----------------+------------------'

		''' Currenty only rev comp sequences are seen in the BBone file, represented by a - but not reversed start / stop '''
		#print node_data_dict

		if int(node_data_dict[node_isolate + '_leftend']) > 0:
			orientation_dict[node_isolate] = '+'
		else:
			orientation_dict[node_isolate] = '-'


		if orientation_dict[node_isolate] == '+':
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
		logging.info('conducting muscle alignment')
		muscle_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'mafft':
		logging.info('conducting mafft alignment')
		logging.info(node_seq_len_est)
		mafft_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'clustalo':
		logging.info('conducting clustal alignment')
		clustalo_alignment('temp_unaligned.fasta', 'temp_aln.fasta')

	if local_aligner == 'kalign':
		logging.info(node_seq_len_est)
		if node_seq_len_est < 7:
			logging.info('conducting mafft alignment')
			mafft_alignment('temp_unaligned.fasta', 'temp_aln.fasta')
		else:
			logging.info('conducting kalign alignment')
			kalign_alignment('temp_unaligned.fasta', 'temp_aln.fasta')


	#fasta_alignment_to_bbone('temp_aln.fasta', 'temp_aln', true_start=node_seq_start_pos)

	new_subgraph = fasta_alignment_to_subnet('temp_aln.fasta', true_start=node_seq_start_pos, node_prefix=node_ID, orientation=orientation_dict, re_link_nodes=False, add_seq=True)

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


def seq_recreate_check(graph_obj, input_dict):
	for isolate in input_dict[1].keys():
		extracted_seq = extract_original_seq(graph_obj, isolate)
		original_seq_from_fasta = input_parser(input_dict[1][isolate])

		count = 0

		while count < len(extracted_seq):
			if extracted_seq[count] != original_seq_from_fasta[0]['DNA_seq'][count]:
				logging.warning(count)
				logging.warning(extracted_seq[count])
				logging.warning(original_seq_from_fasta[0]['DNA_seq'][count])
				logging.warning(extracted_seq[count-10:count + 10])
				logging.warning(original_seq_from_fasta[0]['DNA_seq'][count-10:count + 10])
			count += 1


		if extracted_seq.upper() == original_seq_from_fasta[0]['DNA_seq'].upper():
			logging.info('Sequence recreate pass')
			print('Sequence recreate pass')
			recreate_check_result = 'Pass'
		else:
			logging.error('Sequence recreate fail')
			logging.error(len(extracted_seq))
			logging.error(len(original_seq_from_fasta[0]['DNA_seq']))
			logging.error(extracted_seq[-10:])
			logging.error(original_seq_from_fasta[0]['DNA_seq'][-10:])
			logging.error(extracted_seq[:10])
			logging.error(original_seq_from_fasta[0]['DNA_seq'][:10])
			recreate_check_result = 'Fail'


def add_graph_data(graph_obj):

	count_dict = {}

	# Add start nodes
	for node, data in graph_obj.nodes_iter(data=True):

		logging.info(node)
		logging.info(data)

		for an_isolate in data['ids'].split(','):
			if abs(int(data[an_isolate + '_leftend'])) == 1:
				graph_obj.graph[an_isolate + '_startnode'] = node
				if node not in count_dict.keys():
					count_dict[node] = 1
				else:
					count_dict[node] = count_dict[node] + 1

	logging.info(count_dict)
	most_start_node = ''
	most_start_node_number = 0

	for a_node in count_dict.keys():
		if count_dict[a_node] > most_start_node_number:
			most_start_node = a_node
			most_start_node_number = count_dict[a_node]

	graph_obj.graph['start_node'] = most_start_node


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

	logging.info(path_to_progressiveMauve)
	progressiveMauve_call = [path_to_progressiveMauve, '--output=globalAlignment_' + out_aln_name, '--scratch-path-1=/Volumes/HDD/mauveTemp', '--scratch-path-2=/Volumes/HDD/mauveTemp'] + fasta_path_list

	return call(progressiveMauve_call, stdout=open(os.devnull, 'wb'))

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
					curr_seq_header = '>' + entry[8] + found_in_string + '\n'
					logging.info(curr_seq_header)
					outfile_obj.write(curr_seq_header)
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

	csv_header = 'gene'

	for iso in isolate_list:
		csv_header = csv_header + ',' + iso

	csv_header = csv_header + '\n'

	outfile_obj.write(csv_header)

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
					line_str = csv_header

					if 'ID' in entry[8].keys():
						curr_gene = entry[8]['ID']

					if 'gene_id' in entry[8].keys():
						curr_gene = entry[8]['gene_id']

					if 'locus_tag' in entry[8].keys():
						curr_gene = entry[8]['locus_tag']

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


def create_fasta_from_pangenome_csv(pg_csv, seq_file_dict, out_name):
	in_pg_obj = open(pg_csv, "r")
	out_transcriptome = open(out_name + '.fasta', "w")
	reader = csv.reader(in_pg_obj)
	next(reader, None)

	genome_dict = {}

	for iso_seq_dict in seq_file_dict[1].keys():

		fasta_seq_dict = input_parser(seq_file_dict[1][iso_seq_dict])
		genome_dict[iso_seq_dict] = fasta_seq_dict[0]['DNA_seq']


	all_anno_dict = import_gtf_dict_to_massive_dict(seq_file_dict[3])


	for line in reader:


		if len(line) > 0:
			out_header = ">" + line[0] + '\n'

			seq_details = all_anno_dict[line[0]].split(',')

			#print seq_details

			if seq_details[3] == '-':
				out_seq = genome_dict[seq_details[0]][int(seq_details[1])-1:int(seq_details[2])]
				#out_seq = reverse_compliment(out_seq)
			else:
				out_seq = genome_dict[seq_details[0]][int(seq_details[1])-1:int(seq_details[2])]


			out_transcriptome.write(out_header)
			out_transcriptome.write(out_seq)
			out_transcriptome.write('\n')
			out_transcriptome.write('\n')


def pangenome_csv_to_virtual_genome_fasta(pg_csv, seq_file_dict, out_name, chrom_name='virtChromosome'):
	in_pg_obj = open(pg_csv, "r")
	out_fasta = open(out_name + '.fasta', "w")
	out_gff = open(out_name + '.gff3', "w")
	reader = csv.reader(in_pg_obj)
	next(reader, None)

	genome_dict = {}

	for iso_seq_dict in seq_file_dict[1].keys():

		fasta_seq_dict = input_parser(seq_file_dict[1][iso_seq_dict])
		genome_dict[iso_seq_dict] = fasta_seq_dict[0]['DNA_seq']


	all_anno_dict = import_gtf_dict_to_massive_dict(seq_file_dict[3])

	filler = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

	out_header = ">" + chrom_name + '\n'
	out_fasta.write(out_header)

	out_gff.write('##gff-version 3\n')

	current_pos = len(filler) + 1
	temp_seq = ''

	count = 0

	for line in reader:

		if len(line) > 0:

			seq_details = all_anno_dict[line[0]].split(',')

			out_seq = genome_dict[seq_details[0]][int(seq_details[1])-1:int(seq_details[2])]

			feature_name = line[0]
			feature_start = current_pos
			feature_stop = feature_start + len(out_seq)

			out_gff.write(chrom_name + '\t' + 'gengraph' + '\t' + 'gene' + '\t' + str(feature_start) + '\t' + str(feature_stop) + '\t' + '.' + '\t' + seq_details[3] + '\t' + '.' + '\t' + 'ID=' + feature_name + '\n')
			current_pos = feature_stop + len(filler)

			out_fasta.write(filler)
			temp_seq = temp_seq + filler
			out_fasta.write(out_seq)
			temp_seq = temp_seq + out_seq


def extract_anno_pan_genome_csv(graph_obj, gtf_dict, out_file_name, refseq='', sim_threshold=1.0):

	isolate_list = gtf_dict.keys()



	added_list = []

	outfile_obj = open(out_file_name + 'anno.csv', 'w')

	csv_header = 'gene'

	# Create header for csv file

	for iso in isolate_list:
		csv_header = csv_header + ',' + iso

	csv_header = csv_header + '\n'

	outfile_obj.write(csv_header)

	for isolate in isolate_list:
		gtf_lol = input_parser(gtf_dict[isolate])
		timer = 0

		# For each gene for this isolate, see which other isolates have the same sequence

		for entry in gtf_lol:

			# For this gene for this isolate
			#print entry

			if entry[2] == 'gene':

				# Entries that are genes

				found_in_list = check_isolates_in_region(graph_obj, entry[3], entry[4], isolate, threshold=sim_threshold, return_dict=False)

				if abs(int(entry[4])) < abs(int(entry[3])):
					logging.info(entry)
				# this gene, is also found in these isolates
				logging.info(found_in_list)

				if len(list(set(found_in_list) & set(added_list))) < 1:
					line_str = csv_header

					logging.info(entry[8].keys())

					if 'locus_tag' in entry[8].keys():
						curr_gene = entry[8]['locus_tag']


					logging.info(curr_gene)

					logging.info(line_str)
					#line_str = line_str.replace('gene',curr_gene)
					line_str = line_str.replace(isolate,curr_gene)
					line_str = line_str.replace('gene',curr_gene)
					logging.info(line_str)

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
						logging.info('new pos')
						logging.info(left_pos, right_pos)
						#print 'old pos'
						#print entry[3], entry[4]

						iso_gtf_lol = input_parser(gtf_dict[found_iso])

						#print iso_gtf_lol

						#print found_iso

						if left_pos != 'pos not found' and right_pos != 'pos not found':
							homo_gene = get_anno_from_coordinates(iso_gtf_lol, left_pos[str(found_iso)], right_pos[str(found_iso)], 10)

							logging.info('gene found!!')
							logging.info(homo_gene)
							logging.info(found_iso)

							line_string_list = line_str.split(',')
							for n,i in enumerate(line_string_list):
								if str(i).replace('\n','') == found_iso:
									logging.info('yes')
									line_string_list[n] = homo_gene

							line_str = ','.join(line_string_list)
							logging.info('line string')
							logging.info(line_str)
						else:

							line_str = line_str.replace(found_iso, 'partial')


					for remaining_iso in isolate_list:
						line_string_list = line_str.split(',')
						for n,i in enumerate(line_string_list):
							if i.replace('\n','') == remaining_iso:
								line_string_list[n] = '0'
						line_str = ','.join(line_string_list)

					logging.info('NB OUT --------------------------------------------')
					logging.info(line_str)

					timer += 1

					if line_str[-2:] != '\n':
						line_str = line_str + '\n'

					logging.info('Writing line')
					logging.info(line_str)

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
		'''
		if anno[2] == 'exon':
			if abs(int(anno[3]) - int(start_pos)) <= int(tollerence) and abs(int(anno[4]) - int(stop_pos)) <= int(tollerence):
				return anno[8]['gene_id']
		'''
		# if using gff3
		if anno[2] == 'gene':
			if abs(int(anno[3]) - int(start_pos)) <= int(tollerence) and abs(int(anno[4]) - int(stop_pos)) <= int(tollerence):
				return anno[8]['locus_tag']
			if abs(int(anno[4]) - int(start_pos)) <= int(tollerence) and abs(int(anno[3]) - int(stop_pos)) <= int(tollerence):
				logging.warning('get_anno_from_coordinates problem')
				quit()

	else:
		return '1'


def get_gene_homo_gff(graph_obj, gtf_file, reference_name):

	gff_lod = input_parser(gtf_file)

	for anno in gff_lod:
		start_pos = anno[3]
		stop_pos = anno[4]

		found_in = check_isolates_in_region(graph_obj, start_pos, stop_pos, reference_name)

		logging.info(anno[8], found_in)


def calc_simmilarity_matrix(graph_obj, method='node'):

	import pandas as pd

	iso_list = graph_obj.graph['isolates'].split(',')

	distance_dict = {}

	for ref_iso in iso_list:

		distance_dict[ref_iso] = []

		for other_iso in iso_list:
			node_count = 0
			for a_node, data in graph_obj.nodes_iter(data=True):
				if ref_iso in data['ids'].split(',') and other_iso in data['ids'].split(','):
					node_count += 1
			distance_dict[ref_iso].append(float(node_count))

	distMatrix = pd.DataFrame(distance_dict, index=iso_list)

	#print distMatrix

	self_matrix = {}
	for a_iso in iso_list:
		self_matrix[a_iso] = distMatrix[a_iso][a_iso]

	#print pd.DataFrame(self_matrix, index=self_matrix.keys())


	#print distMatrix / pd.DataFrame(self_matrix, index=self_matrix.keys())

	return distMatrix / pd.DataFrame(self_matrix, index=self_matrix.keys())


# --------------------- Sequence extraction related

def extract_original_seq(graph_obj, seq_name):
	"""
	Given a networkx graph object created by GenGraph, return the sequence of an individual isolate
	as a string. This can be used to extract a genome for a particular isolate.
	:param graph_obj: A networkx graph object created by or formatted in a similar way to GenGraph.
	:param seq_name: The identifier of the sequence to be extracted, as was used in the sequence file and found in the
	node attributes.
	:return: A string.
	"""
	extracted_seq = ''

	# A position list of lists (lol)
	pos_lol = []

	test_count = 0

	for node, data in graph_obj.nodes_iter(data=True):

		is_negative = False

		left_end_name = seq_name + '_leftend'
		right_end_name = seq_name + '_rightend'

		if seq_name in data['ids'].split(','):
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
			logging.info('-')
			logging.info(segment)

		if segment[1] == last_node_line[0]:
			logging.info('+')
			logging.info(segment)

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


def extract_original_seq_region(graph_obj, region_start, region_stop, seq_name):
	"""
	Retrieve a subsequence from the graph for an isolate given the positions. This can be used for example to extract
	a gene given the gene start and stop positions using the normal positions found in an gff file.
	:param graph_obj: A networkx graph object created by or formatted in a similar way to GenGraph.
	:param region_start: The nucleotide position to start extracting from
	:param region_stop: The nucleotide position to stop extracting from
	:param seq_name: The identifier of the sequence that the positions are associated with (The isolate).
	:return: A string
	"""

	seq_string = extract_original_seq(graph_obj, seq_name)

	region_start = int(region_start) - 1
	region_stop = int(region_stop)

	return seq_string[region_start:region_stop]


def extract_original_seq_region_fast(graph_obj, region_start, region_stop, seq_name):
	"""
	A Newer version of this function, attempting to find a faster way.
	Requires testing for negative nodes.
	Retrieve a subsequence from the graph for an isolate given the positions. This can be used for example to extract
	a gene given the gene start and stop positions using the normal positions found in an gff file.
	:param graph_obj: A networkx graph object created by or formatted in a similar way to GenGraph.
	:param region_start: The nucleotide position to start extracting from
	:param region_stop: The nucleotide position to stop extracting from
	:param seq_name: The identifier of the sequence that the positions are associated with (The isolate).
	:return: A string
	"""

	region_start = int(region_start)
	region_stop = int(region_stop)

	multi_node_dict = {}

	for node, data in graph_obj.nodes_iter(data=True):

		if seq_name in data['ids'].split(','):

			node_start = abs(int(data[seq_name + '_leftend']))
			node_stop = abs(int(data[seq_name + '_rightend']))

			if int(data[seq_name + '_leftend']) < 0:
				orientation = '-'
			else:
				orientation = '+'

			if region_start > node_start and region_stop < node_stop:

				if orientation == '-':
					# Reverse comp seq stuff
					print('still to code')
					quit()

				else:
					slce_start = node_start - region_start
					slce_end = region_stop - node_stop

					return data['sequence'][slce_start:slce_end + 1]

			elif node_start <= region_start <= node_stop:

				#print('start node')

				if orientation == '-':
					# Reverse comp seq stuff
					print('still to code')
					quit()
				else:

					sub_seq = data['sequence'][region_start - node_start:]
					multi_node_dict[sub_seq] = node_start

			elif node_start <= region_stop <= node_stop:

				#print('end node')

				if orientation == '-':
					# Reverse comp seq stuff
					print('still to code')
					quit()
				else:

					sub_seq = data['sequence'][0: region_stop - node_start + 1]
					multi_node_dict[sub_seq] = node_start

			elif region_start <= node_start and node_stop <= region_stop:

				#print('mid node')

				if orientation == '-':
					# Reverse comp seq stuff
					print('still to code')
					quit()
				else:

					multi_node_dict[data['sequence']] = node_start

	result_string = ''
	for key, value in sorted(multi_node_dict.iteritems(), key=lambda (k, v): (v, k)):
		result_string += key

	return result_string


def extract_iso_subgraph(graph_obj, isolate):
	iso_graph = nx.MultiDiGraph()
	iso_graph.graph['isolate'] = 'isolate'

	for node, data in graph_obj.nodes_iter(data=True):
		if isolate in data['ids'].split(','):
			#print node
			iso_graph.add_node(node, data)
			# For NetworkX 2
			#iso_graph.add_node(node)
			#nx.set_node_attributes(iso_graph, {node: data})

	iso_graph = link_nodes(iso_graph, isolate)

	return iso_graph


def extract_gene(seq_locus_id, seq_isolate_origin, graph_obj, annotation_path_dict):


	iso_anno_obj = input_parser(annotation_path_dict[3][seq_isolate_origin])


	tar_gene_anno = 'Not found'

	for entry in iso_anno_obj:
		if entry[2] == 'gene':
			if entry[8]['locus_tag'] == seq_locus_id:
				tar_gene_anno = entry

			if 'old_locus_tag' in entry[8].keys():
				if entry[8]['old_locus_tag'] == seq_locus_id:
					tar_gene_anno = entry


	if tar_gene_anno != 'Not found':

		logging.info(tar_gene_anno[3], tar_gene_anno[4])
		logging.info(int(tar_gene_anno[4]) - int(tar_gene_anno[3]))
		logging.info(tar_gene_anno[6])

		out_seq = extract_original_seq_region(graph_obj, tar_gene_anno[3], tar_gene_anno[4], seq_isolate_origin)

		if tar_gene_anno[6] == '-':
			out_seq = reverse_compliment(out_seq)

		return out_seq

	else:
		return tar_gene_anno

	logging.info('in function')


# For heaviest path function

def get_neighbour_most_iso(list_of_nodes, graph_obj, weight_matrix):

	# Returns the neighbouring node of the current node that contains the most isolates or the highest weight.
	longest_list_node = 'nope'
	longest_list_length = 0

	isolate_list = graph_obj.graph['isolates'].split(',')


	if len(weight_matrix) == 0:
		for neigh_node in list_of_nodes:

			if len(graph_obj.node[neigh_node]['ids'].split(',')) > longest_list_length:
				longest_list_node = neigh_node
				longest_list_length = len(graph_obj.node[neigh_node]['ids'].split(','))

	else:
		# Using the weight matrix

		# Calculating the arverage distance (Maybe move if too time consuming)
		ave_dist_dict = calc_average_distance_dict(weight_matrix, isolate_list)
		largest_node_weight = 0

		for neigh_node in list_of_nodes:
			#print neigh_node
			if 'visited' not in graph_obj.node[neigh_node].keys():
				#print 'checking node'
				node_weight = 0

				for iso_name in graph_obj.node[neigh_node]['ids'].split(','):
					node_weight = node_weight + ave_dist_dict[iso_name]

				if node_weight > largest_node_weight:
					longest_list_node = neigh_node
					largest_node_weight = node_weight
			else:
				1 == 1
				#print 'visited'

	if longest_list_node != 'nope':
		graph_obj.node[longest_list_node]['visited'] = 'yes'

	return longest_list_node


def calc_average_distance_dict(weight_matrix, a_list_of_isolates):
	# Return the average distance from

	res_dict = {}

	sum_of_distances = 0

	for an_iso in a_list_of_isolates:
		#print '\n'
		#print an_iso
		#print weight_matrix[an_iso]
		#print sum(weight_matrix[an_iso]) - weight_matrix[an_iso][an_iso]

		res_dict[an_iso] = (sum(weight_matrix[an_iso]) - weight_matrix[an_iso][an_iso]) / (float(len(a_list_of_isolates)) - 1.0)


	return res_dict


def extract_heaviest_path(graph_obj, start_node, stop_node, weight_matrix=''):
	# setup
	out_graph = nx.MultiDiGraph()

	# Adding, and let's go!

	node_list = [start_node]

	curr_node = start_node

	while curr_node != stop_node:

		neighbors_out = graph_obj.successors(curr_node)

		new_node = get_neighbour_most_iso(neighbors_out, graph_obj, weight_matrix)

		curr_node = new_node

		logging.info(curr_node)

		if curr_node != 'nope':
			node_list.append(new_node)


		if curr_node == 'nope':
			curr_node = node_list[-2]
		logging.info(curr_node)


	node_list = node_list[:-1]

	for heavy_node in node_list:
		out_graph.add_node(heavy_node, graph_obj.node[heavy_node])
		#out_graph.add_node(heavy_node)
		#nx.set_node_attributes(out_graph, {heavy_node: graph_obj.node[heavy_node]})



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

# Extracting a branch sequence file

# Shandu's code


def extract_max(graph_obj, node, extract_size, region):
	'''Extract the largest possible sequence from the beginning or end of a node.

	Arguments for the region parameter:
	   beginning: extract sequence from beginning of node
	   end: extract sequence from end of node

	If the size to be extracted from the node (extract_size) is longer than the sequence of the node, this returns the entire sequence of the node.

	Otherwise, this returns only the length of the sequence specified by extract_size.
	'''
	if len(graph_obj.node[node]['sequence']) >= extract_size:
		if region == 'end':
			return graph_obj.node[node]['sequence'][-extract_size:]
		else:
			return graph_obj.node[node]['sequence'][:extract_size]
	else:
		return graph_obj.node[node]['sequence']


def extract_branch_seq(graph_obj, out_file_name, extract_size):
	'''Create a file listing the different sequence versions at each branch of the graph.

	extract_size specifies the maximum sequence length to be extracted from each node of the graph.
	'''

	# For python 2.7
	import Queue as queue
	#For python 3
	#import queue

	# Changed error in formatting fasta file (missing newline and extra space in the header)

	start_node = graph_obj.graph['start_node']		# This should be an attribute of the graph_obj (currently uses the start node of new_fun_3_genome)

	size = extract_size*2
	nx.set_node_attributes(graph_obj, 'tag', 'unreached')
	out_file = open(out_file_name, 'w')
	ready = queue.Queue()
	graph_obj.node[start_node]['tag'] = 'reached'
	ready.put(start_node)
	while not ready.empty():
		current_node = ready.get()
		successors = graph_obj.successors(current_node)
		for suc_node in successors:
			if graph_obj.node[suc_node]['tag'] == 'unreached':
				graph_obj.node[suc_node]['tag'] = 'reached'
				ready.put(suc_node)
			details = current_node + '-' + suc_node
			sequence = extract_max(graph_obj, current_node, extract_size, 'end')
			remaining_size = size - len(sequence)
			sequence += extract_max(graph_obj, suc_node, remaining_size, 'beginning')
			# Determine if sequence is long enough
			if len(sequence) < size:
				# Sequence is too short
				details_initial = details
				sequence_initial = sequence

				more_paths = list()
				# List of node paths from which sequences will be extracted to be added to sequence_initial

				paths = list()
				# List of node paths to which additional nodes will be added until the path contains enough nodes to yield a long enough sequence. At this point the path will be added to more_paths.
				paths.append([suc_node])

				while len(paths) != 0:
					path = paths.pop()
					sequence = ''
					size_remain = remaining_size
					for node in path:
						sequence += extract_max(graph_obj, node, size_remain, 'beginning')
						size_remain = remaining_size - len(sequence)
					# Determine if sequence is long enough
					if size_remain > 0:
						# Sequence is too short
						if graph_obj.successors(path[-1]) == []:
							# Reached the end of graph
							more_paths.append(path)
						else:
							# Extend the path to include the next successor (a new path is created for each successor of the last node in the current path)
							successors = graph_obj.successors(path[-1])
							for suc_node in successors:
								new_path = path + [suc_node]
								paths.append(new_path)
					else:
						# Sequence is long enough
						more_paths.append(path)

				# For each path in more_paths, extract the required sequence and append it to sequence_initial
				for path in more_paths:
					details = details_initial
					sequence = sequence_initial
					remaining = remaining_size
					for node in path[1:]:	# path[0] is the successor (suc_node) of current_node so its sequence has already been extracted
						remaining = size - len(sequence)
						if remaining > 0:
							details += '-' + node
							sequence += extract_max(graph_obj, node, remaining, 'beginning')
					out_file.write('>Branch-' + details + '\n')
					out_file.write(sequence + '\n\n')
			else:
				# Sequence is long enough
				out_file.write('>Branch-' + details + '\n')
				out_file.write(sequence + '\n\n')
	out_file.close()


def get_branch_mapping_dict(path_to_edge_file):
	''' Take the file generated by samtools and convert to a dict for the edge pairs '''
	res_dict = {}

	aln_path_obj = open(path_to_edge_file)

	for line in aln_path_obj:
		if len(line) > 1:
			splitline = line.split('\t')
			node_info = splitline[0]
			node_coverage = splitline[2]

			splitNodeInfo = node_info[7:]

			res_dict[splitNodeInfo] = node_coverage

	return res_dict

'''
def find_best_aln_subpaths_old(edge_aln_dict):

	logging.info("Finding best subpaths")

	best_aln_paths = []

	for edgeLabel in edge_aln_dict.keys():
		other_path_available = False
		edgeLabel_split = edgeLabel.split('-')
		print edgeLabel
		for other_edgeLabel in edge_aln_dict.keys():
			other_edgeLabel_split = other_edgeLabel.split('-')

			if edgeLabel_split[0] == other_edgeLabel_split[0] and edgeLabel_split[-1] == other_edgeLabel_split[-1]:
				if edgeLabel != other_edgeLabel:
					other_path_available = True
					if int(edge_aln_dict[edgeLabel]) > int(edge_aln_dict[other_edgeLabel]):
						best_aln_paths.append(edgeLabel)
		if other_path_available == False:
			best_aln_paths.append(edgeLabel)

	return best_aln_paths
'''


def find_best_aln_subpaths(edge_aln_dict, coverage_threshold):

	logging.info("Finding best subpaths")

	best_aln_paths = []

	best_aln_paths_dict = {}

	# Pre-filtering
	edge_aln_dict_filtered = {}

	for edgeLabel in edge_aln_dict.keys():
		if int(edge_aln_dict[edgeLabel]) > int(coverage_threshold):
			edge_aln_dict_filtered[edgeLabel] = edge_aln_dict[edgeLabel]

	logging.info('Pre-filter = ', str(len(edge_aln_dict)))
	logging.info('Post-filter = ', str(len(edge_aln_dict_filtered)))
	timeCount = 0
	total_time = len(edge_aln_dict_filtered)


	for edgeLabel in edge_aln_dict_filtered.keys():

		edgeLabel_split = edgeLabel.split('-')

		edge_key = edgeLabel_split[0] + '-' + edgeLabel_split[-1]

		if edge_key not in best_aln_paths_dict.keys():

			best_aln_paths_dict[edge_key] = edgeLabel

		else:
			# If this edge already exists, compare weights and only add the largest one.

			if int(edge_aln_dict_filtered[edgeLabel]) > int(edge_aln_dict_filtered[best_aln_paths_dict[edge_key]]):

				best_aln_paths_dict[edge_key] = edgeLabel

		timeCount += 1

	# Now add the best paths to the list
	for edgeLabel in best_aln_paths_dict.keys():
		best_aln_paths.append(best_aln_paths_dict[edgeLabel])


	return best_aln_paths


def get_path_weight(path_list, aGraph):

	total_weight = 0

	test_path = aGraph.subgraph(path_list)

	for anEdge in list(test_path.edges_iter(data=True)):

		total_weight = total_weight + anEdge[2]['weight']

	return total_weight


def retrieve_genomic_sequence(bp_start, bp_stop, fasta_object):

	if int(bp_start) > 0:
		a_sequence = fasta_object[0]['DNA_deq'][bp_start - 1:bp_stop]

	else:
		a_sequence = fasta_object[0]['DNA_deq'][abs(bp_start) - 1:abs(bp_stop)]


	return a_sequence


def create_new_graph_from_aln_paths(graph_obj, aln_path_obj, path_dict, trim=True):

	# This needs to go
	startNode = 'Aln_61_1'
	endNode = 'Aln_387_1'

	logging.info('Creating new graph from edge mapping results')
	newIsoGraph = nx.Graph()

	logging.info("Add new edges")

	for node_path in aln_path_obj:

		node_path_list = node_path.split('-')

		newIsoGraph.add_path(node_path_list, aln_isolate='shortPath',weight=int(path_dict[node_path]))


	# Filter out alignments starting from a high threshold, and lowering it untill there is a path from the start to the stop node
	# Test trimming all one degree nodes, and re-adding the start - stop nodes
	# Move the "Add sequences" step to the end to lower memory usage
	# If that doesn't work, do the heuristic stepwise graph traversal

	if trim == True:
		# find nodes with 3 edges, then kill the one with 1

		remove_node_list = []

		for a_Node, data in newIsoGraph.nodes_iter(data=True):
			if newIsoGraph.degree(a_Node) == 3:
				# Get neighbours
				all_neighbours_of_node = nx.all_neighbors(newIsoGraph, a_Node)

				for node_neighbour in all_neighbours_of_node:

					if newIsoGraph.degree(node_neighbour) == 1:
						remove_node_list.append(node_neighbour)

		for trimnode in remove_node_list:
			newIsoGraph.remove_node(trimnode)


	# Extract path

	nx.write_graphml(newIsoGraph, 'newPaths.xml')

	logging.info("calculate heaviest path")

	heaviest_path_weight = 0

	for path in nx.all_simple_paths(newIsoGraph, startNode, endNode):
		logging.info('start')
		pathWeight = get_path_weight(path, newIsoGraph)
		logging.info(pathWeight)
		if pathWeight > heaviest_path_weight:
			heaviest_path_weight = pathWeight
			heaviest_path = path


	heavy_graph = newIsoGraph.subgraph(heaviest_path)

	nx.write_graphml(heavy_graph, 'heaviestPath.xml')

	#Add the new path to the old graph


	logging.info("Add sequences")

	for seqNode, data in heavy_graph.nodes_iter(data=True):


		heavy_graph.node[seqNode]['sequence'] = graph_obj.node[seqNode]['sequence']



	# Extract sequence

	nodesInOrder = nx.shortest_path(heavy_graph, source=startNode, target=endNode)

	heavtSeq = ''

	seqDict = nx.get_node_attributes(heavy_graph,'sequence')
	for node in nodesInOrder:
		logging.info(node)
		heavtSeq = heavtSeq + seqDict[node]

	return heavtSeq


# Ancestral genome creation

def get_end_node_dict(graph_obj):
	end_dict_value = {}
	end_dict_node = {}

	for an_isolate in graph_obj.graph['isolates'].split(','):
		end_dict_value[an_isolate] = 0

	for node, data in graph_obj.nodes_iter(data=True):
		for an_isolate in graph_obj.graph['isolates'].split(','):
			if an_isolate in data['ids'].split(','):
				if abs(int(data[an_isolate + '_rightend'])) > end_dict_value[an_isolate]:
					end_dict_value[an_isolate] = abs(int(data[an_isolate + '_rightend']))
					end_dict_node[an_isolate] = node

	return end_dict_node


def generate_ancesteral_genome(graph_obj, weight_matrix=''):

	# Remove all nodes with only one isolate in them (Simplify graph)

	out_path = extract_heaviest_path(graph_obj, graph_obj.graph['start_node'], 'Aln_700_1', weight_matrix=weight_matrix)
	#print len(out_path)
	out_path.graph['start_node'] = graph_obj.graph['start_node']
	nx.write_graphml(out_path, 'ancestor_path.xml')

	seq_string = ''
	for node, data in out_path.nodes_iter(data=True):
		seq_string = seq_string + data['sequence']
	#print len(seq_string)

	return out_path

	'''
	quit()
	print get_end_node_dict(graph_obj)
	max_nodes = 5000

	#print([p for p in nx.all_shortest_paths(graph_obj,source=graph_obj.graph['start_node'],target='Aln_700_1')])

	p=nx.shortest_path(graph_obj)
	print p[graph_obj.graph['start_node']]['Aln_700_1']

	quit()
	path_count = 0
	for a_path in nx.all_simple_paths(graph_obj, graph_obj.graph['start_node'], 'Aln_700_1', cutoff=max_nodes):
		path_count += 1

	print path_count
	'''


def add_ancestral_path(old_graph_obj, anc_graph_obj):

	iso_node_count = {}

	for a_node, data in anc_graph_obj.nodes_iter(data=True):

		for a_iso in old_graph_obj.node[a_node]['ids'].split(','):
			if a_iso in iso_node_count.keys():
				iso_node_count[a_iso] = iso_node_count[a_iso] + 1
			else:
				iso_node_count[a_iso] = 1

		old_graph_obj.node[a_node]['ids'] = data['ids'] + ',' + 'ancestral'
		#print old_graph_obj.node[a_node]['ids']


	#print iso_node_count
	return old_graph_obj


def get_panTrans_stats(in_annoTransCSV):
	''' Get general stats from the output of the pan-transcriptome generation '''
	csv_obj = open(in_annoTransCSV, 'r')

	header_line = True

	total_genes_count = 0
	core_genes_count = 0

	for line in csv_obj:

		is_core_gene = True

		if header_line == True:
			header = line.split(',')[1:]
			is_core_gene = False
			header_line = False
		else:

			line_list = line.split(',')[1:]

			for list_item in line_list:
				if list_item == '0':
					is_core_gene = False

			'''
			for list_item in line_list:
				if list_item == '1':
					is_core_gene = False
			'''

			for list_item in line_list:
				if list_item == 'partial':
					is_core_gene = False



		if is_core_gene == True:
			core_genes_count = core_genes_count + 1

		total_genes_count = total_genes_count + 1

	logging.info('Total genes', total_genes_count)
	logging.info('Core genes', core_genes_count)

# ----------------------------------------------------- # Development code


def split_all_long_nodes(a_in_graph, max_length):
	logging.info('splitting all nodes')

	node_list = a_in_graph.nodes()

	for a_node in node_list:
		#print 'current node: ' + a_node
		a_in_graph = split_node(a_in_graph, a_node, max_length)

	return a_in_graph


def split_node(in_graph, node, max_length):

	logging.info('splitting node: ' + node)

	#print 'testing'

	# Extract info from node

	#print in_graph.node[node]

	graph_iso_list = in_graph.node[node]['ids'].split(',')

	len_dict = {}

	longest_seq = 0

	for a_isolate in graph_iso_list:
		a_length = abs(abs(int(in_graph.node[node][a_isolate + '_rightend'])) - abs(int(in_graph.node[node][a_isolate + '_leftend'])))

		len_dict[a_isolate] = a_length

		if a_length > longest_seq:
			longest_seq = a_length

	#print len_dict
	#print longest_seq

	# determine split

	if longest_seq < max_length:
		return in_graph

	new_nodes = int(longest_seq) // int(max_length) + 1
	#print new_nodes

	new_node_dict = {}
	count = 1
	while count <= new_nodes:

		new_node_dict[node + '_' + str(count)] = {'ids':in_graph.node[node]['ids']}

		count += 1

	#print new_node_dict

	for a_isolate in graph_iso_list:

		if int(in_graph.node[node][a_isolate + '_rightend']) < 0:
			# dealing with reverse comp
			logging.info('next')


		else:
			# Not reverse comp
			#print a_isolate

			node_interval = len_dict[a_isolate] // new_nodes
			node_length_remainder = len_dict[a_isolate] % new_nodes
			#print node_interval
			#print node_length_remainder

			currnode = 1

			# Initialise with the start position
			prev_node_stop_plus_one = int(in_graph.node[node][a_isolate + '_leftend'])

			#print prev_node_stop_plus_one

			while currnode <= new_nodes:

				curr_node_start = prev_node_stop_plus_one

				curr_node_stop = curr_node_start + node_interval - 1

				# Taking care of the remainders
				if currnode <= node_length_remainder + 1:
					curr_node_stop = curr_node_stop + 1


				prev_node_stop_plus_one = curr_node_stop + 1

				new_node_dict[node + '_' + str(currnode)][a_isolate + '_leftend'] = str(curr_node_start)
				new_node_dict[node + '_' + str(currnode)][a_isolate + '_rightend'] = str(curr_node_stop)

				currnode += 1


			#print curr_node_start, curr_node_stop
			#print int(in_graph.node[node][a_isolate + '_rightend'])
			#print new_node_dict




	# Create new nodes

	for a_new_node in new_node_dict.keys():
		#print a_new_node
		in_graph.add_node(a_new_node, new_node_dict[a_new_node])
		#in_graph.add_node(a_new_node)
		#nx.set_node_attributes(in_graph, {a_new_node: new_node_dict[a_new_node]})


	# Delete old node

	in_graph.remove_node(node)


	# Link nodes
	#for a_isolate in graph_iso_list:
	#	in_graph = link_nodes_2(in_graph, a_isolate)

	return in_graph

# ---------------------------------------------------- # Testing functions


def seq_addition_test(in_graph, node_ID, seq_fasta_paths_dict):

	# Get seq from node
	pass_seq_test = True

	the_node_seq = in_graph.node[node_ID]['sequence']
	logging.info('node seq')
	logging.info(the_node_seq)

	node_seq_iso_list = in_graph.node[node_ID]['ids'].split(',')

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
			logging.info(seq_isolate + ' - Passed')
			logging.info(sub_iso_seq)
		else:
			pass_seq_test = False
			logging.info(seq_isolate + ' - Failed ' + str(start_pos) + ' - ' + str(end_pos))
			logging.info(sub_iso_seq)

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

	for n, d in in_graph.nodes_iter(data=True):
		if 'sequence' in d.keys():
			comp_len = comp_len + len(d['sequence'])

			for node_iso in d['ids'].split(','):
				len_dict[node_iso] = len_dict[node_iso] + len(d['sequence'])


	#print comp_len

	report_file.write("\nCompressed sequence length: " + str(comp_len))

	#print len_dict

	for iso_name in len_dict.keys():
		report_file.write('\n' + iso_name + " length: " + str(len_dict[iso_name]))

	#print "Density: " + str(nx.density(in_graph))

	report_file.write("\nDensity: " + str(nx.density(in_graph)))

	return nx_summary





