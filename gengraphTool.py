
from gengraph import *

import pkg_resources

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='''Welcome to GenGraph v0.1''', epilog="""Tools for the creation and use of graph genomes""")

	parser.add_argument('toolkit', type=str, default='test_mode', help='Select the tool you would like to use')

	parser.add_argument('--out_file_name', type=str, help='Prefix of the created file')

	parser.add_argument('--out_file_path', type=str, default='./', help='Output file destination')

	parser.add_argument('--alignment_file', nargs=1, help='The path to the alignment file')

	parser.add_argument('--backbone_file', nargs=1, default='default', help='The path to the backbone file')

	parser.add_argument('--out_format', nargs=1, default='default', help='Format for the output')

	parser.add_argument('--block_aligner', nargs=1, default=['progressiveMauve'], help='Block aligner to use')

	parser.add_argument('--progressiveMauve_path', nargs=1, default=['progressiveMauve'], help='Path to progressiveMauve if not in PATH')

	parser.add_argument('--node_msa_tool', nargs=1, default='mafft', help='MSA tool to use')

	parser.add_argument('--seq_file', type=str, help='Tab delimited text file with paths to the aligned sequences')

	parser.add_argument('--no_seq', dest='should_add_seq', action='store_false',
						help='Create a graph genome with no sequence stored within')

	parser.add_argument('--make_circular', type=str, default='No',
						help='To circularise the graph for a sequence, give the name of that sequence')

	parser.add_argument('--recreate_check', dest='rec_check', action='store_true',
						help='Set to True to attempt to recreate the input sequences from the graph and compare to the originals')

	parser.add_argument('--extract_sequence', type=str, default='some_isolate',
						help='Returns the sequence of the selected isolate')

	parser.add_argument('--isolate', type=str, default='some_isolate', help='pass the isolate variable. For graph generation, this should be the genome that best represents the ancesteral state.')

	parser.add_argument('--extract_sequence_range', nargs=2, default=['all', 'all'],
						help='Extract sequence between two positions')

	parser.add_argument('--graph_file', type=str, help='Give the path to the graph file')

	parser.add_argument('--max_node_length', type=int, default=-1, help='Max sequence length that can be aligned')

	parser.add_argument('--input_file', type=str, help='Generic input')

	parser.add_argument('--locus_ID', type=str, help='The name of the gene or feature')

	parser.add_argument('--scratch_path', type=str, default='./', help='Path for temporary data')

	parser.set_defaults(should_add_seq=True)

	parser.set_defaults(rec_check=False)

	args = parser.parse_args()

	print('Running GenGraph Toolkit')

	# Setting up logging

	if not args.out_file_name:
		logging.basicConfig(filename='new_run.log', level=logging.DEBUG)
	else:
		logging.basicConfig(filename=args.out_file_name + '.log', level=logging.DEBUG)
	# Check NetworkX version
	nx_version = pkg_resources.get_distribution("networkx").version

	logging.info('NetworkX version used:' + nx_version)

	if args.toolkit == 'test_mode':
		print("Test functions here")

		test_aln_graph = nx.read_graphml(args.graph_file)

		parsed_input_dict = parse_seq_file(args.seq_file)

		result = seq_recreate_check(test_aln_graph, parsed_input_dict)

		print(result)

	if args.toolkit == 'make_genome_graph':
		# Requires:
		# --out_file_name
		# --seq_file

		# optional:
		# --recreate_check
		# --no_seq

		print('Creating genome graph')

		# Clean this up
		global_aligner = args.block_aligner
		local_aligner = args.node_msa_tool

		path_to_progressiveMauve = args.progressiveMauve_path[0]

		start_time = time.time()

		if not args.seq_file:
			print('A sequence file needs to be specified with the --seq_file flag for graph creation.')
			quit()

		if not args.out_file_name:
			args.out_file_name = 'default'

		parsed_input_dict = parse_seq_file(args.seq_file)

		check_result = input_file_check(parsed_input_dict)

		if len(check_result) > 0:
			for input_error in check_result:
				logging.error(input_error)
				print(input_error)
			quit()

		# --------------------------------------------------------------------------------- Initial global alignment

		if args.block_aligner[0] == 'progressiveMauve' and args.backbone_file == 'default':
			print('Conducting progressiveMauve')

			logging.info(parsed_input_dict)

			progressiveMauve_alignment(path_to_progressiveMauve, parsed_input_dict[2], args.out_file_name, scratch=args.scratch_path)

			logging.info('progressiveMauve Complete')
			print('progressiveMauve Complete')

		# --------------------------------------------------------------------------------- Conversion to block graph

		if args.backbone_file == 'default':
			bbone_file = 'globalAlignment_' + args.out_file_name + '.backbone'
		else:
			print('Using existing BBone file')
			bbone_file = args.backbone_file[0]

		logging.info('Running bbone_to_initGraph')

		genome_aln_graph = bbone_to_initGraph(bbone_file, parsed_input_dict)

		refine_initGraph(genome_aln_graph)

		add_missing_nodes(genome_aln_graph)

		nx.write_graphml(genome_aln_graph, 'intermediate_Graph.xml')

		# --------------------------------------------------------------------------------- node splitting

		if args.max_node_length != -1:
			genome_aln_graph = split_all_long_nodes(genome_aln_graph, args.max_node_length)

			print('Nodes split')

			nx.write_graphml(genome_aln_graph, 'intermediate_split_Graph.xml')

		# --------------------------------------------------------------------------------- Local node realignment

		print('Conducting local node realignment')

		genome_aln_graph = realign_all_nodes(genome_aln_graph, parsed_input_dict)

		add_graph_data(genome_aln_graph)

		genome_aln_graph = link_all_nodes(genome_aln_graph)

		print('Genome graph created')

		# --------------------------------------------------------------------------------- adding annotation data to the graph
		# TEMP FIX FOR NOW

		parsed_input_dict

		# --------------------------------------------------------------------------------- adding sequence to the graph

		if args.should_add_seq:

			if args.isolate[0] == 'some_isolate':
				ref_isolate = graph_obj.graph['isolates'].split(',')[0]

			print('Sequence dict used:')
			print(parsed_input_dict[1])

			genome_aln_graph = add_sequences_to_graph(genome_aln_graph, parsed_input_dict)
			print('Sequence added')

		if args.make_circular != 'No':
			make_circular(genome_aln_graph, args.make_circular)
			print('Graph circularised')

		# Do the recreate check to see if the original sequence is correctly recalled from the graph
		if args.rec_check:
			print('Doing recreate check')
			seq_recreate_check(genome_aln_graph, parsed_input_dict)

		# Saving output
		out_put_dir = args.out_file_path
		if args.out_file_path[-1] != '/':
			out_put_dir += '/'

		if args.out_format == 'default':

			out_filename_created = out_put_dir + args.out_file_name + '.xml'
			nx.write_graphml(genome_aln_graph, out_filename_created)

		if args.out_format[0] == 'serialize':
			print('Writing to serialized file')
			pickle.dump(genome_aln_graph, open(out_put_dir + args.out_file_name + '.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

		end_time = (time.time() - start_time)

		print("run time: " + str(end_time))
		generate_graph_report(genome_aln_graph, out_put_dir + args.out_file_name)

	if args.toolkit == 'make_graph_from_fasta':

		# Requires:
		# --input_file
		# --out_file_name



		fasta_object = input_parser(args.input_file)

		print('Adding', len(fasta_object), 'sequences')

		seqStartDict = {}

		for seq_entry in fasta_object:
			seqStartDict[seq_entry['gene_details']] = 1

		new_graph = fasta_alignment_to_subnet(args.input_file, true_start=seqStartDict, add_seq=True)

		nx.write_graphml(new_graph, 'intermediate_virus_Graph.xml')

		# new_graph = add_sequences_to_graph(new_graph, fasta_object)

		nx.write_graphml(new_graph, args.out_file_name + '.xml')

	if args.toolkit == 'region_alignment_score':

		graph_obj = nx.read_graphml(args.graph_file)


		result = check_isolates_in_region(graph_obj, args.extract_sequence_range[0], args.extract_sequence_range[1],
										  args.isolate, threshold=1.0, return_dict=True)

		print('result is')
		print(result)

	if args.toolkit == 'extract_fasta_file':

		out_fasta = open(args.out_file_name, 'w')

		graph_obj = nx.read_graphml(args.graph_file)

		extracted_seq = extract_original_seq(graph_obj, args.isolate)

		fasta_headder = '>' + args.isolate + '\n'

		out_fasta.write(fasta_headder)

		n = 70
		for seq_line in [extracted_seq[i:i + n] for i in range(0, len(extracted_seq), n)]:
			out_fasta.write(seq_line + '\n')

		out_fasta.close()

	if args.toolkit == 'extract_region':
		imported_genome = nx.read_graphml(args.graph_file)

		print("Extracting from " + str(args.extract_sequence_range[0]) + " to " + str(args.extract_sequence_range[1]))

		print(extract_original_seq_region(imported_genome, args.extract_sequence_range[0], args.extract_sequence_range[1],
									args.isolate))

	if args.toolkit == 'extract_ancesteral_genome':

		# Requires:
		# --out_file_name
		# --graph_file

		print("CHANGES MADE TO THE MATRIX ARE PROBLEMATIC AND WILL CAUSE A INCORRECT TREE TO BE MADE")

		imported_genome = nx.read_graphml(args.graph_file)

		sim_matrix = calc_simmilarity_matrix(imported_genome)

		# sim_matrix = (sim_matrix - 1) * -1

		plotDend = True
		add_to_GG = True

		print(sim_matrix.index.values.tolist())

		print(sim_matrix.as_matrix())

		if plotDend == True:
			from scipy.cluster import hierarchy
			import matplotlib.pyplot as plt

			Z = hierarchy.linkage(sim_matrix.as_matrix(), 'single')

			plt.figure()
			dn = hierarchy.dendrogram(Z, labels=sim_matrix.index.values.tolist())
			plt.show()

		print(sim_matrix)

		quit()

		anc_genome_obj = generate_ancesteral_genome(imported_genome, weight_matrix=sim_matrix)

		# ancesteral_genome = extract_heaviest_path(imported_genome, args.isolate, weight_matrix=sim_matrix)

		if add_to_GG == True:
			fresh_imported_genome = nx.read_graphml(args.graph_file)
			anc_genome_added_graph = add_ancestral_path(fresh_imported_genome, anc_genome_obj)
			nx.write_graphml(anc_genome_added_graph, args.out_file_name + 'Ancesteral.xml')

		ancesteral_genome_seq = extract_seq_heavy(anc_genome_obj)

		export_to_fasta(ancesteral_genome_seq, args.out_file_name, args.out_file_name)

	if args.toolkit == 'extract_pan_transcriptome':

		# This needs to deal with the required GTF dict and seq_file_dict

		# Needs
		# --graph_file
		# --seq_file
		# --out_file_name

		print('Extracting...')

		parsed_input_dict = parse_seq_file(args.seq_file)

		graph_obj = nx.read_graphml(args.graph_file)

		test_gtf_dict = parsed_input_dict[3]

		if args.isolate == 'some_isolate':
			ref_isolate = ''
		else:
			ref_isolate = args.isolate

		# Extracting the csv from the graph
		print('Extracting annotated pan genome csv')
		extract_anno_pan_genome_csv(graph_obj, test_gtf_dict, args.out_file_name, sim_threshold=0.95)
		print('Extracting pan genome csv')
		extract_pan_genome_csv(graph_obj, test_gtf_dict, args.out_file_name, hom_threshold=0.95, refseq=ref_isolate)
		print('Extracting pan genome transcriptome')
		create_fasta_from_pangenome_csv(args.out_file_name + '.csv', test_gtf_dict, parsed_input_dict,
										args.out_file_name)

	# Converting the csv to a fasta file of transcripts
	# create_fasta_from_pangenome_csv(args.input_file, test_gtf_dict, parsed_input_dict, args.out_file_name)

	if args.toolkit == 'extract_gene':
		'''Return the sequence of a gene'''

		# --seq_file
		# --graph_file
		# --locus_ID
		# --isolate

		gene_isolate = args.isolate
		gene_name = args.locus_ID
		create_fasta_file = False

		print('Here we go')

		parsed_seq_obj = parse_seq_file(args.seq_file)

		imp_genome_obj = nx.read_graphml(args.graph_file)

		print(extract_gene(gene_name, gene_isolate, imp_genome_obj, parsed_seq_obj))

	if args.toolkit == 'map_to_graph':
		# Requires:
		# --out_file_name
		# --graph_file

		print("Creating branch mapping file")

		imported_graph_obj = nx.read_graphml(args.graph_file)

		imported_graph_obj.graph['start_node'] = 'Aln_61_1'

		# extract_branch_seq(imported_graph_obj, args.out_file_name, 20)

		# Link the above to bellow later

		res = get_branch_mapping_dict('/Volumes/HDD/Genomes/M_tuberculosis/gg_genomes/alnRestult.txt')
		# print res
		aln_path_list = find_best_aln_subpaths(res, 20)

		print(str(len(aln_path_list)) + ' paths extracted')

		imported_genome = nx.read_graphml(args.graph_file)

		newGraphSeq = create_new_graph_from_aln_paths(imported_genome, aln_path_list, res)

		export_to_fasta(newGraphSeq, args.out_file_name, args.out_file_name)

	if args.toolkit == 'generate_report':
		imported_genome = nx.read_graphml(args.graph_file)

		generate_graph_report(imported_genome, args.out_file_name)
