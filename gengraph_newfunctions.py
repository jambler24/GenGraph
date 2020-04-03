# This is where new untested functions sit until ready for deploy


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
