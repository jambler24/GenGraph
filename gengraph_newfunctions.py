# This is where new untested functions sit until ready for deploy
import itertools
import gengraph as geng
import numpy as np
from scipy.stats import rankdata as rd
import pandas as pd
import os
from multiprocessing import Pool, cpu_count

# --------------------- PanGenome related

def extract_pan_genome(graph_obj, gtf_dict, out_file_name):
	"""Currently not functioning"""

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
	"""Extract the largest possible sequence from the beginning or end of a node.

	Arguments for the region parameter:
		beginning: extract sequence from beginning of node
		end: extract sequence from end of node

	If the size to be extracted from the node (extract_size) is longer than the sequence of the node, this returns the entire sequence of the node.

	Otherwise, this returns only the length of the sequence specified by extract_size.
	"""
	if len(graph_obj.node[node]['sequence']) >= extract_size:
		if region == 'end':
			return graph_obj.node[node]['sequence'][-extract_size:]
		else:
			return graph_obj.node[node]['sequence'][:extract_size]
	else:
		return graph_obj.node[node]['sequence']


def extract_branch_seq(graph_obj, out_file_name, extract_size):
	"""Create a file listing the different sequence versions at each branch of the graph.

	extract_size specifies the maximum sequence length to be extracted from each node of the graph.
	"""

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
	""" Take the file generated by samtools and convert to a dict for the edge pairs """
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


# ------------------------------------------------- Ancestral genome creation

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


# ----------------------------------------------------- # Development code - old and unused

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


def calc_GC(a_sequence):

	gc_perc = (a_sequence.count("G") + a_sequence.count("C")) / len(a_sequence)

	return gc_perc

def hash_seq(a_sequence):

	import hashlib


	return hash_val


def create_encode_dict():
	import string

	#string.printable

	index_count = 0

	length = 4

	encoding_dict = {}

	possible_chars = ['A', 'C', 'G', 'T']

	end_char = '-'

	for comb in itertools.combinations_with_replacement(possible_chars, length):
		encoding_dict[''.join(comb)] = string.printable[index_count]
		index_count += 1




	return encoding_dict


def char_encode(a_sequence, encode_dict):



	return 'dvdjvjvf'


def create_freagment_similarity_table(sequence_file, window_size, step_size):

	from gengraph import input_parser
	import pandas as pd

	seq_info = input_parser(sequence_file)

	print(seq_info)

	enc_dict = create_encode_dict()

	print(enc_dict)

	quit()

	column_names = ['fragment_id', 'GC', 'hash', 'char_encode']

	fragment_table = pd.DataFrame(columns = column_names)

	window_size = 50
	stepsize = 20

	for a_seq_file in seq_info:

		# Open file and get windows

		seq_strings = input_parser(a_seq_file['seq_path'])

		for a_contig in seq_strings:

			it = iter(a_contig['DNA_seq'])

			count = 0

			result = tuple(islice(it, window_size))

			if len(result) == window_size:
				print(''.join(result))
				print('----')

				# The first fragment

			for elem in it:

				result = result[1:] + (elem,)

				if count == stepsize:
					count = 0
					result = result[1:] + (elem,)
					kmer_string = ''.join(result)

					# Here we are dealing with a fragment

					fragment_id = 's2f5'
					GC_content = calc_GC(kmer_string)
					hashval = 'oolkjh'
					char_encoded_value = char_encode(kmer_string, )

					new_df_row = pd.DataFrame([[fragment_id, GC_content, hashval, char_encoded_value]], columns=column_names)

					fragment_table = fragment_table.append(new_df_row, ignore_index=True)

					print(fragment_table)

				count += 1




	return 'dataframe'


def global_align_by_composition(sequence_file):


	return 'k'

#----------------------------------------------------------------------------------------------------------------
#Phillip code
#Alignment Functions

def simpleIsSmaller(x, y):
	'''
	Small helper function for suffix_array_skew. Lexicographic comparison for strings of length <= 2.
	:x: First string being compared.
	:y: Second string.
	:return: Boolean, True if x < y, False otherwise.
	'''	
	if x[0] < y[0]:
		return True
	elif (x[0] == y[0] and x[1] < y[1]):		
			return True		
	return False


def suffix_array_skew(s, check = False):
	'''
	Skew algorithm for linear time suffix array construction, courtesy of Juha Karkkainen and Peter Sanders.
	Notes referenced: http://www.mi.fu-berlin.de/wiki/pub/ABI/SS13Lecture3Materials/script.pdf
	:s: String or numpy char array that suffix array is being constructed from.	
	:check: True if array contains integers (ranks), False if strings. There is probably a smarter way to do this.
	:return: Suffix array of s.
	'''
	n = len(s)
	#Checks if s is a string, if so convert to numpy char array. Adds stop characters.
	if (isinstance(s, str)):	
		s = s + "$$$"
		s = np.char.array([x for x in s])
	else:
		s = np.concatenate((s, np.array([-1,-1,-1])), axis=0)	
	
	#parameters used later in algorithm.
	n0 = int((n+2)/3)
	n1 = int((n+1)/3)
	n2 = int(n/3)
	n02 = n0+n2
	#starting position of substrings of length = 3, where T[i,i+1,i+2] with i != 0 (mod 3).
	#Basically all 3-mers, but only if start position of 3-mer is not divisible by 3.
	triple_positions = np.array([x for x in range(1, n+(n0-n1)) if x%3 != 0])

	#Functions for retrieving triples from array.
	if not (check):
		y = lambda a: (s[a:a+3])
	else:
		y = lambda a: [int(x) for x in (s[a:a+3])]

	#Using pandas to retrieve the rankings (lexicographic) of the triples. Conserves order of original array.
	#Ties are allowed and ranking starts at 0. For example: Rankings = [0,1,1,0,0]
	triples = pd.DataFrame.from_records(np.array([y(r) for r in triple_positions]))		
	cols = triples.columns.tolist()	
	triples = (triples[cols].apply(tuple, axis=1).rank(method="dense").astype(int)).values-1

	#This is just re-ordering the ranking array. Required for recursion.
	#Thus: T = triples starting at i%3 == 1, and triples starting at i%3 == 2.	
	t1 = triples[::2]
	t2 = triples[1::2]
	T = np.concatenate((t1,t2), axis = 0)

	#if all ranking name not unique (ties): recursively call algorithm again, else we can calc suffix array of T.	
	if not (np.unique(T).size == T.size):			
		print(n0)
		T = suffix_array_skew(T, True)			
		suff = [1+3*x if x < n0 else 2+3*(x-n0) for x in T]			
	else:
		#construct suffix array of T
		#value of triples = new index of values of triple_positions at same index		
		n_trip = len(triples)
		suff = [0]*n_trip		
		for i in range(n_trip):
			suff[triples[i]] = triple_positions[i]				

	#inducing last third of suffixes (those with starting position i%3==0), and sorting them.
	d0 = {x-1: s[x-1] for x in suff if x%3 == 1}	
	d_sorted = {k: v for k, v in sorted(d0.items(), key=lambda item: item[1])}	
	suff2 = list(d_sorted.keys())	

	#We are storing the indices of the triples in triple_positions.
	t = lambda a: int((a-1)/3) if a%3==1 else int((a-1)/3) + n0		
	for i in range(n02):
		triple_positions[t(suff[i])] = i

	#now merge suff and suff2 into one suffix array
	#renaming terms to keep in line with notes (so its easier to understand following code if you reference notes)
	a0, a12, og_s, n, n0, ind = suff2, suff, s, n, int((n+2)/3), triple_positions

	#u, v = maximum values for pointers iterating over a0 and a12
	u = len(a0)
	v = len(a12)
	#defining inverse function, inverse of a12, such that inv(a12(i)) = i
	r12 = lambda a: ind[int((a-1)/3)] if (a%3==1 and a < n) else (ind[int((a-1)/3) +n0] if a < n else 0)

	#x = pointer for a0
	x = 0
	#y = pointer for a12. Skip first value of 12 if n%3 == 1
	if n%3 == 1:
		y = 1
	else:
		y = 0

	out = []
	#iterating over a0 and 12.
	while (x < u) and (y < v):
		i = a0[x]
		j = a12[y]			
		#see notes for details. this code determines if i or j should be appended to suffix array first.
		if j%3 == 1:
			if (og_s[i] < og_s[j]) or (og_s[i] == og_s[j] and r12(i+1) < r12(j+1)):				
				out.append(i)
				x += 1
			else:											
				out.append(j)
				y += 1
		else:
			t1 =  og_s[i:i+2]
			t2 =  og_s[j:j+2]

			if (simpleIsSmaller(t1,t2)) or ((np.array_equal(t1, t2)) and r12(i+2) < r12(j+2)):
				out.append(i)
				x += 1
			else:
				out.append(j)
				y += 1
	#appending the remaining suffixes.
	for e in range(x, u):
		out.append(a0[e])
	for f in range(y, v):
		out.append(a12[f])

	final_suff = out
	return final_suff	

def bin_search(T, sa, q_list):
	'''
	Simple binary search algorithm to find exact matches using suffix arrays.
	Prints out locations of exact matches for each query.
	:param T: Text being searched
	:param sa: Suffix array of text
	:param q_list: List containing query sequences as strings. Ex: ["ATGTC","GTGTCA"]
	'''
	sa_len = len(sa)
	#starting location for binary search = middle.
	current_p = int(sa_len /2)	

	#iterating over list of queries.
	for q in q_list:
		print("-------------------")
		print(q)
		l = len(q)
		#Starting range = entire array
		max_p = len(sa)
		min_p = 0

		#does search until min and max are equal
		while (abs(max_p-min_p)>1):
			#c_suff = current suffix being compared (extracted from T)
			c_suff = T[sa[current_p]:]
			#r = string being compared to q, same length
			r = c_suff[:l]		

			#if string matches prefix of suffix: print location of this suffix, and
			#extend upwards and downwards for additional exact matches (since suffix array is sorted all matches will be neighbours)
			if q == r:								
				print(sa[current_p])
				#move up while matches hit and print match locations
				temp_p = current_p+1				
				while (temp_p < sa_len):														
					r = T[sa[temp_p]:][:l]
					if q == r:
						print(sa[temp_p])
						temp_p += 1
					else:
						break
				#move down while matches hit and print match locations
				temp_p = current_p-1
				while(temp_p > -1):										
					r = T[sa[temp_p]:][:l]
					if q == r:
						print(sa[temp_p])
						temp_p -= 1
					else:
						break
				break

			#If query is smaller: make current position the maximum of new search range. 
			elif q < r:				
				max_p = current_p
				current_p = int((max_p+min_p)/2)
			#If query bigger: make current position minimum of new search range.
			elif q > r:				
				min_p = current_p
				current_p = int((max_p+min_p)/2)


def exact_match(graph, genome_name, query_list):
	'''
	Prints locations (in linear sequence of specific genome) of exact matches with list of queries, using suffix arrays of reference.
	This could be used for finding initial seeds for a seed-and-extend approach to alignment.
	This should later be changed to return some dictionary or list containing the locations.
	:param graph: GenGraph graph object
	:param genome_name: Name of genome of interest as string. Ex: "Beijing"
	:param query_list: List containing query sequences as strings. Ex: ["ATGTC","GTGTCA"]
	'''
	print("loading graph, extracting seq")	
	T = geng.extract_original_seq(graph, genome_name)
	T_array = np.array([x for x in T])
	print("seq extracted")
	#this prevents command prompt from freezing (in windows this is a problem when code takes >5 minutes to run)
	os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = '1'
	print("Calculating suffix array of seq (this will take a while, +- 5-10 mins):")
	sa = suffix_array_skew(T_array)
	del T_array
	print("Suffix array done, starting exact match search: ")
	bin_search(T, sa, query_list)

#Functions for inexact match following BWA backward search algorithm
#Algorithm can be found here https://academic.oup.com/bioinformatics/article/25/14/1754/225615
#Authors: Heng Li, Richard Durbin

def bwt_from_sa(T, sa):
	'''
	Calculates Burrows Wheeler transform from Text and Suffix array.
	:param T: Input text
	:param sa: Suffix Array of T
	:return: BWT(text), string.
	'''
	bw = []
	for s in sa:
		if s == 0:
			bw.append('$')
		else:
			bw.append(T[s-1])
	return ''.join(bw)

def BWT_sa(text):
	'''
	Calculates Burrows Wheeler transform from Text.
	:param text: Input text
	:return: BWT(text), string.
	'''
	text = text + "$"	
	input_array = np.array([x for x in text])
	sa = suffix_array_skew(input_array)
	b  = bwt_from_sa(text, sa)
	del sa

	return b

def calc_FM_index(B):
	'''
	Calculates fm index and offset array. Offset array = c. c[a] = number of characters smaller than a in text.

	:param B: Burrows-Wheeler transform of a string.
	'''
	set_b = set(B)
	len_b = len(B)
	alphabet_size = len(set_b)
	alph_temp = {a:0 for a in set_b}
	d = {k: v for k, v in sorted(alph_temp.items(), key=lambda item: item[0])}			
	fm = np.zeros((len_b, alphabet_size), dtype= np.uint32)	

	#counting occurences of character a in B[0,i]
	for i in range(len_b):
		d[B[i]] += 1		
		fm[i,] = list(d.values())	
	
	c = np.zeros(alphabet_size, dtype=np.uint32)	
	for i in range(1, alphabet_size):		
		c[i] = sum(fm[-1,:i])	

	out = (fm, c)	
	
	return out

def InexRecur(W,i,z,k,l, D, c, fm, insertion_penalty = 1, deletion_penalty = 1):
	'''
	Recursive function used by inexact_match algorithm.

	:param W: Substring being search for
	:param i: Current index of W being compared. Since this is a backward search, i starts at end of W.
	:param z: Mismatches, starts at max differences, and gets reduced each time mismatch is found.
	:param k: Initial floor of search range
	:param l: Initial ceiling of search range
	:param D: Array of size |W|, used to prune search tree and reduce amount of comparisons.
	:param fm: FM-index of text being searched.
	:param insertion_penalty: Difference score for an insertion, can be tweaked to reflect biological mutation rates. Default = 1.
	:param deletion_penalty: Difference score for a deletion, can be tweaked to reflect biological mutation rates. Default = 1.

	:return: Set of suffix array indices which contain matches.
	'''
	tempset = set()
	
	#D[i] = 0 if i < 0, else D[i]. Pruning search.
	if (i < 0):
		maxDifferences = 0
	else:
		maxDifferences = D[i]

	#Base cases for recursion, stops if z (mismatches) or i (index) reach allowed minimum.
	#Return empty if maxDifferences have been reached, i.e. no match.
	if (z < maxDifferences):		
		return set()
	#i < 0 implies we have backward searched through whole query. We are done, return range k->l+1.
	if i < 0:	
		for m in range(k,l+1):
			tempset.add(m)								
		return tempset

	result = set()
	#if insertions allowed, add results from searching i-1, with insertion penalty added.
	result = result.union(InexRecur(W, i-1, z-insertion_penalty, k, l, D, c, fm))

	order = {'$':0, 'A':1,'C':2,'G':3,'T':4}
	t = lambda y: order[y]	

	#If W[i] matches b (characters: $,A,C,G,T), go to next letter without mismatch penalty. Else add mismatch penalty.
	for b in list(order.keys())[1:]:
		#New search range, changing k (floor) and l (ceiling)
		#fm[i] = 0 for i < 0, hence the two if statements.
		if (k-1 < 0):
			newMin = c[t(b)] + 1
		else:
			newMin = c[t(b)] + fm[k-1, t(b)] + 1
		if (l < 0):
			newMax = c[t(b)] 
		else:
			newMax = c[t(b)] + fm[l, t(b)]		

		if newMin <= newMax:				
			result = result.union(InexRecur(W, i, z-deletion_penalty, newMin, newMax, D, c, fm))
			if b == W[i]:
				#Match found - go to next letter			
				result = result.union(InexRecur(W, i-1, z, newMin, newMax, D, c, fm))
			else:
				#No Match found - go to next letter but substract one from max differences allowed.				
				result = result.union(InexRecur(W, i-1, z-1, newMin, newMax, D, c, fm))

	return result


def inexact_search(X, W, z):
	'''
	Basic implementation of inexact match algorithm found in BWA, accounting for mismatches/substitutions, as well as insertions + deletions.
	Courtesy of Heng Li and Richard Dubin. Paper found here: https://academic.oup.com/bioinformatics/article/25/14/1754/225615
	This implementation is complete but slow. The full BWA algorithm uses various heuristics and tricks to reduce computation time.
	One rather easy improvement would be to implement seeding. This can be done by writing a function which considers the first m characters of a query only.
	Where m = seed length. Then allow only up to n differences for the seed. Then only align the full reads which have seeds that find matches are aligned.
	
	:param X: Text (genome sequence) being searched
	:param W: Query. Substring (read sequence) being searched for (aligned/mapped)
	:param z: Maximum amount of mismatches being allowed for (substitutions, indels)
	:return: Set containing indices of Suffix Array that match query.
	'''
	#Preprocessing steps. 
	#Using multiprocessing to calculate BWT of Text and Reversed Text simultaneously.
	p = Pool(cpu_count())
	data = p.map(BWT_sa, [X, X[::-1]])
	p.close()	
	print("BWT done")	
	B = data[0]
	B_rev = data[1]	

	#Calculaying FM-Index and offset-array for BWT(X) and BWT(X_reversed)
	p = Pool(cpu_count())
	data = p.map(calc_FM_index, [B, B_rev])
	p.close()
	print("FM index done")

	fm = data[0][0]
	c = data[0][1]-1
	fm_r = data[1][0]		

	#Calculating D (array used to prune search-tree and reduce runtime)
	#Variables used for calculating D
	k = 0	#Floor
	l = len(B)-1	#Ceiling
	zz = 0
	alph = set(W)	#Alphabet
	D = [0]*len(W)

	#Alphabet dictionary. Hardcoded for genomes, but can calculate it from set(W) to generalise.
	#The alternative is to just use pandas so u have named columns in the first place.
	order = {'$':0, 'A':1,'C':2,'G':3,'T':4}	
	t = lambda y: order[y]
	
	for i in range(0, len(W)):
		if (k-1 < 0):
			k = c[t(W[i])] + 1
		else:
			k = c[t(W[i])] + fm_r[k-1, t(W[i])] + 1	
		if (l < 0):
			l = c[t(W[i])]
		else:
			l = c[t(W[i])] + fm_r[l, t(W[i])]

		if k > l:
			k = 0
			l = len(B)-1
			zz = zz+1
		D[i] = zz

	#Now that we have D, can use recursion function to find indices in Suffix Array that contains matches.
	sa_range = InexRecur(W,len(W)-1,z,0,len(B)-1, D, c, fm)	

	return sa_range