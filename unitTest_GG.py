# Dependencies

import networkx as nx

# Built in


import unittest

from subprocess import call

import argparse

import sys

import csv

import copy


# Import GenGraph functions


from gengraph import *


''' Basic tests '''


class BasicFunctionsTestCase(unittest.TestCase):
	"""Tests for basic functions in `phil_0_4.py`."""

	def test_reverse(self):
		"""Is the string correctly reversed?"""
		self.assertEqual(reverse('ATGTTA'),'ATTGTA')

	def test_compliment_Base(self):
		"""Is the base correctly complimented?"""
		self.assertEqual(compliment_Base('A'),'T')
		self.assertEqual(compliment_Base('T'),'A')
		self.assertEqual(compliment_Base('G'),'C')
		self.assertEqual(compliment_Base('C'),'G')
		self.assertEqual(compliment_Base('N'),'N')

	def test_compliment_DNA(self):
		"""Is the string correctly complimented?"""
		self.assertEqual(compliment_DNA('ATGTTA'),'TACAAT')

	def test_reverse_compliment(self):
		"""Is the string correctly reverse complimented?"""
		self.assertEqual(reverse_compliment('ATGTTA'),'TAACAT')

	def test_is_even(self):
		"""Is the int even?"""
		self.assertTrue(is_even(8))
		self.assertFalse(is_even(7))
		self.assertTrue(is_even(-8))
		self.assertFalse(is_even(-7))


''' Graph utility function tests'''

'''
class GraphUtilityFunctionsTestCase(unittest.TestCase):
	"""Tests for graph utility functions in `phil_0_4.py`."""

	#test_graph = './test_files/unit_tests/aln_add3_seq.xml'
	#imported_genome = nx.read_graphml(test_graph)

	def test_nodes_connected(self):
		"""Are linked nodes correctly detected? (directly linked neighbours)"""
		test_graph = './test_files/unit_tests/aln_add3_seq.xml'
		imported_genome = nx.read_graphml(test_graph)
		self.assertTrue(nodes_connected('X_20', 'X_18', imported_genome))
		self.assertFalse(nodes_connected('X_18', 'X_20', imported_genome))
		self.assertFalse(nodes_connected('X_21', 'X_18', imported_genome))
		self.assertFalse(nodes_connected('X_20', 'X_20', imported_genome))
'''

''' Graph generation function tests'''


class GraphGeneratingFunctionsTestCase(unittest.TestCase):
	"""Tests for graph generation functions in `phil_0_4.py`."""

	def test_alignment_to_graph_conversion(self):
		"""Are fasta alignment files correctly converted to graphs?"""

		test_fasta_aln = './unitTestFiles/alignedSeq.fa'

		test_start_dict = {'seq1':1,'seq2':-50,'seq3':1}

		fasta_object = input_parser(test_fasta_aln)

		new_graph = fasta_alignment_to_subnet(test_fasta_aln, true_start=test_start_dict, node_prefix='X', orientation={}, re_link_nodes=True, add_seq=True)

		nx.write_graphml(new_graph, './unitTestFiles/alignedSeq.xml')

		seq1_seq = extract_original_seq(new_graph, 'seq1')
		seq2_seq =  extract_original_seq(new_graph, 'seq2')
		seq3_seq =  extract_original_seq(new_graph, 'seq3')

		does_pass = True

		print seq1_seq
		print 'GAGATTAGGAGTAGATAGATAGATATTTAGAGCCCGGAAAATTTATATTATTTAAT'
		print seq2_seq
		print reverse_compliment('GATTAGGAGTAGATAGATAGATATTTAGAGAGAGAAAATTTATATTATTT')
		print seq3_seq
		print 'GAGATTAGGAGTAGATAGATAGTATTTAGAGAGAGAAAATTTATATTATTTAAT'

		if seq1_seq != 'GAGATTAGGAGTAGATAGATAGATATTTAGAGCCCGGAAAATTTATATTATTTAAT':
			does_pass = False
			print 'seq1 fails'

		if seq2_seq != reverse_compliment('GATTAGGAGTAGATAGATAGATATTTAGAGAGAGAAAATTTATATTATTT'):
			does_pass = False
			print 'seq2 fails'

		if seq3_seq != 'GAGATTAGGAGTAGATAGATAGTATTTAGAGAGAGAAAATTTATATTATTTAAT':
			does_pass = False
			print 'seq3 fails'		

		self.assertTrue(does_pass)

	def test_create_genome_alignment_graph(self):
		"""Is the correct block graph created?"""

		test_bbone_file = './unitTestFiles/globalAlignment_phyloTest.backbone'

		parsed_input_dict = parse_seq_file('./unitTestFiles/multiGenome.txt')



		created_test_graph = create_genome_alignment_graph(test_bbone_file, parsed_input_dict[0], parsed_input_dict[1])

		#nx.write_graphml(created_test_graph, './test_files/unit_tests/bbone_to_graph_test_out.xml')

		test_pass = True

		#print test_pass

		self.assertTrue(test_pass)

	'''

	def test_split_node(self):
		"""Are nodes being correctly split?"""

		curr_in_graph = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/test_files/two_genome/test_graph_2.xml'

		node_to_split = 'H37Rv_seq_19'

		test_length = 10000

		imported_genome = nx.read_graphml(curr_in_graph)

		print 'genome import complete'

		#split_g = split_node(imported_genome, node_to_split, test_length)

		#split_g = split_all_long_nodes(imported_genome, test_length)


		Rv_seq = input_parser('/Volumes/HDD/Genomes/M_tuberculosis/H37Rv/mycobacterium_tuberculosis_h37rv_2_supercontigs.fasta')
		Ra_seq = input_parser('/Users/panix/M_tuberculosis/H37Ra/mycobacterium_tuberculosis_h37ra_1_contigs.fasta')

		print '\n--------------------- Recall test ---------------------'
		print len(Rv_seq[0]['DNA_seq'])
		print len(Ra_seq[0]['DNA_seq'])

		print Rv_seq[0]['DNA_seq'][483399:483409]
		#print Ra_seq[0]['DNA_seq'][2132526:2132536].upper()

		seq_added_file = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/seq_added_TEST_LALALA.xml'

		seq_added_genome = nx.read_graphml(seq_added_file)

		print 'Origional length:'
		print len(Rv_seq[0]['DNA_seq'])
		print 'extracted length:'
		#extracted_seq = extract_seq(seq_added_genome, 'H37Rv')
		#print len(extracted_seq)

		
		base_count = 0
		while base_count < len(extracted_seq):
			if Rv_seq[0]['DNA_seq'][base_count:base_count+1] != extracted_seq[base_count:base_count+1]:
				print 'Error: ' + str(base_count)
				base_count += 1


			else:
				base_count += 1

		

		#if extracted_seq == Rv_seq[0]['DNA_seq']:
		#	print 'WIN!!!'


		
		for iso in imported_genome.graph['isolates'].split(','):
			print iso

			imported_genome = link_nodes_2(imported_genome, iso)

		nx.write_graphml(imported_genome, 'linked_test_all.xml')
		

		prob_node_file = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/test_files/two_genome/intermediate_split_Graph.xml'
		seq_data_file = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/test_files/two_genome/mtb_2genome_file.txt'

		prob_genome_obj = nx.read_graphml(prob_node_file)
		seq_fasta_paths_dict = parse_seq_file(seq_data_file)[1]

		prob_node = 'Aln_region_2_42'

		realign_prob_g = local_node_realign(prob_genome_obj, prob_node, seq_fasta_paths_dict)

		nx.write_graphml(realign_prob_g, 'prob_graph_2_42_test.xml')

		is_split = True

		self.assertTrue(is_split)

	def fasta_alignment_to_subnet(self):
		"""Are the fasta files correctly converted to graphs?"""

		print 'HDFGSDFGADFGADFGASFGASDGASDF'

		curr_in_graph = '/Users/panix/Dropbox/Programs/tools/genome_alignment_graph_tool/test_files/two_genome/temp_aln.fasta'

		new_graph = fasta_alignment_to_subnet(curr_in_graph)

		nx.write_graphml(new_graph, 'graph_fasta_file.xml')

		is_working = True

		self.assertTrue(is_working)



	def test_link_nodes(self):
		"""Are nodes correctly linked?"""

		test_graph_unlinked = './test_files/unit_tests/unlinked_node_test.xml'
		test_graph_linked = './test_files/unit_tests/bbone_to_graph_test_comparison.xml'

		imported_unlinked_genome = nx.read_graphml(test_graph_unlinked)
		imported_linked_genome = nx.read_graphml(test_graph_linked)

		iso_list = ['Alpha','Beta','Cappa']

		for isolate in iso_list:
			imported_unlinked_genome = link_nodes(imported_unlinked_genome, isolate)

		test_pass = nx.is_isomorphic(imported_unlinked_genome, imported_linked_genome)

		self.assertTrue(test_pass)
'''


if __name__ == '__main__':
	unittest.main()

