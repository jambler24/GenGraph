# Dependencies

import networkx as nx
from gengraphTool import *

# Built in

from subprocess import call
from operator import itemgetter

import argparse
import sys
import os
import csv
import copy
import time
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align
import pandas as pd
import difflib

try:
	import cPickle as pickle
except ModuleNotFoundError:
	import pickle

import logging


path_to_muscle = 'muscle3.8.31_i86darwin64'

path_to_clustal = 'clustalo'

path_to_mafft = 'mafft'

path_to_kalign = 'kalign'

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

		# Old slow version
		#seq_string = extract_original_seq_region(self, region_start, region_stop, seq_name)

		# Changed to fast version
		seq_string = extract_original_seq_region_fast(self, region_start, region_stop, seq_name)

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

		total_annos = len(anno_dict)
		success_count = 0

		for line in anno_dict:

			conv_coords = convert_coordinates(self, line[3], line[4], source_seq_id, target_seq_id)

			if conv_coords is not None:

				success_count += 1

				line[0] = chr_name
				line[3] = str(conv_coords[target_seq_id + '_leftend'])
				line[4] = str(conv_coords[target_seq_id + '_rightend'])
				info_dict = line[8]
				line[8] = 'ID=' + info_dict['ID'] + ';' + 'Name=' + info_dict['Name']

				line_string = '\t'.join(line)
				out_file.write(line_string)

		res_string = str(success_count) + '/' + str(total_annos) + ' annotations transferred'

		return res_string

	def ids(self):
		"""
		Returns the IDs of the sequences found in the graph.
		:return: list of ID strings
		"""

		# Need to change to ids to be consistant with graph
		ID_list = self.graph['isolates'].split(',')

		ID_list = [x.encode('UTF8') for x in ID_list]

		return ID_list

	def get_region_subgraph(self, region_start, region_stop, seq_name, neighbours=0):
		"""
		Return a sub-graph from a defined region with positions relative to some of the sequences.
		:param region_start: Start position (bp)
		:param region_stop: Stop position (bp)
		:param seq_name: Sequence ID that the bp positions are relative to.
		:param neighbours: Number of neighbours to be included. Currently testing, only 1 level available.
		:return: A GenGraph sub-network representing the region.
		"""
		sub_graph = extract_region_subgraph(self, region_start, region_stop, seq_name, neighbours=neighbours)

		return sub_graph

	def plot_subgraph(self, region_start, region_stop, seq_name, neighbours=0):
		"""
		Return a plot of a sub-graph from a defined region with positions relative to some of the sequences.
		:param region_start: Start position (bp)
		:param region_stop: Stop position (bp)
		:param seq_name: Sequence ID that the bp positions are relative to.
		:param neighbours: Number of neighbours to be included. Currently testing, only 1 level available.
		:return: Plots a netowrkx + matplotlib figure of the region subgraph.
		"""

		sub_graph = extract_region_subgraph(self, region_start, region_stop, seq_name, neighbours=neighbours)
		pos = nx.spring_layout(sub_graph)
		nx.draw_networkx_nodes(sub_graph, pos, cmap=plt.get_cmap('jet'), node_size=300)
		nx.draw_networkx_labels(sub_graph, pos)
		nx.draw_networkx_edges(sub_graph, pos, arrows=True)
		plt.show()

	def region_similarity(self, region_start, region_stop, ref_name, alt_name, method='shared'):
		"""
		Info text
		:param region_start:
		:param region_stop:
		:param ref_name:
		:param alt_name:
		:param method:
		:return:
		"""
		if method == 'shared':

			sub_graph = extract_region_subgraph(self, region_start, region_stop, ref_name, neighbours=0)

			total_ref_length = region_stop - region_start
			total_lone_ref_length = 0

			for a_node in sub_graph.nodes:

				if alt_name not in sub_graph.node[a_node]['ids'].split(','):
					print(a_node)
					print(sub_graph.node[a_node]['ids'])

					total_lone_ref_length += len(sub_graph.node[a_node]['sequence'])

			print(total_ref_length)
			print(total_lone_ref_length)
			sim_perc = (total_ref_length - total_lone_ref_length) / (total_ref_length) * 100

			return sim_perc

		else:

			print(method + ' not implimented.')

			return 0

	# ---------------------------------------------------- New functions under testing

	# ---------------------------------------------------- Read Alignment functions
	kmerDict = {}
	global firstPath
	firstPath = True

	def fast_kmer_create(self, kmerLength):
		"""
		:param kmerLength: A number that indicated the length of the kmers that you want as an output
		:return: A dictionary with all the kmers in the graph, each entry in the dictionary has a kmer name, the nodes that the kmer covers, the sequence of the kmer and
				the starting position of the kmer on the graph
		"""
		# Start Positions stored in the dictionary are character list positions, so the first character in the node sequence is position 0 not position 1
		graphNodes = list(self.nodes)
		kmerNumber = 1

		# Runs through the graph nodes creating kmers and checking if each of the nodes has an inversion in it
		for i in graphNodes:
			currentNodeName = str(i)
			idsCounter = str(self.node[currentNodeName]['ids']).count(',')
			idsCounter += 1
			negCounter = str(self.node[currentNodeName]).count('-')
			if idsCounter * 2 == negCounter:
				ids = str(self.node[currentNodeName]['ids']).split(',')
				for f in ids:
					self.node[currentNodeName][str(f) + '_leftend'] = int(
						self.node[currentNodeName][str(f) + '_leftend']) * -1
					self.node[currentNodeName][str(f) + '_rightend'] = int(
						self.node[currentNodeName][str(f) + '_rightend']) * -1
				self.node[currentNodeName]['sequence'] = str(self.node[currentNodeName]['sequence'])[::-1]

			nodeSequence = self.nodes[currentNodeName]['sequence']
			revNodeSequence = nodeSequence[::-1]
			nodeKeys = list(self.nodes[currentNodeName])
			kmerStartPos = 0
			values = list(dict(self.node[currentNodeName]).values())
			inv = False
			for num in values:
				if isinstance(num, int):
					if num < 0:
						inv = True
			if inv:
				for j in range(len(nodeSequence)):
					if kmerStartPos + kmerLength <= len(nodeSequence):
						kmerSeq = nodeSequence[kmerStartPos:kmerStartPos + kmerLength]
						kmerRevSeq = revNodeSequence[kmerStartPos:kmerStartPos + kmerLength]
						dictStartPos = kmerStartPos
						kmerStartPos += 1
						tempDict = {'nodes': [currentNodeName], 'sequence': kmerSeq, 'inversesequence': kmerRevSeq,
									'startposition': dictStartPos}
						newkmer = {'kmer_' + str(kmerNumber): tempDict}
						kmerNumber += 1
						GgDiGraph.kmerDict.update(newkmer)
					else:
						kmerPartseq = nodeSequence[kmerStartPos:]
						dictStartPos = kmerStartPos
						kmerStartPos += 1
						currentNode = [currentNodeName]
						self.kmer_over_multiple_nodes(kmerPartseq, currentNode, kmerLength, kmerNumber, dictStartPos,
													  True)
						kmerNumber += 1
			else:
				for j in range(len(nodeSequence)):
					if kmerStartPos + kmerLength <= len(nodeSequence):
						kmerSeq = nodeSequence[kmerStartPos:kmerStartPos + kmerLength]
						dictStartPos = kmerStartPos
						kmerStartPos += 1
						tempDict = {'nodes': [currentNodeName], 'sequence': kmerSeq, 'startposition': dictStartPos}
						newkmer = {'kmer_' + str(kmerNumber): tempDict}
						kmerNumber += 1
						GgDiGraph.kmerDict.update(newkmer)
					else:
						kmerPartseq = nodeSequence[kmerStartPos:]
						dictStartPos = kmerStartPos
						kmerStartPos += 1
						currentNode = [currentNodeName]
						self.kmer_over_multiple_nodes(kmerPartseq, currentNode, kmerLength, kmerNumber, dictStartPos,
													  False)
						kmerNumber += 1
		return GgDiGraph.kmerDict

	def kmer_over_multiple_nodes(self, kmerPart, currentNode, kmerLength, kmerNumber, kmerStartPos, hasInv):
		"""
		You do not call this function directly. This function is used when creating kmers from the graph and is called by
		the fast_kmer_create function. Recursively goes through nodes until the kmer has been successfully created.
		:param kmerPart: the current part of the kmer sequence
		:param currentNode: the current node that the kmer is traversing over
		:param kmerLength: the length of the desired outcome for the kmer
		:param kmerNumber: the current number of this kmer
		:param kmerStartPos: The starting position of this current kmer
		:param hasInv: Boolean as to whether the node has an inversion
		:return: Adds to the kmerDict of all kmers in the graphs. This is called by the fast kmer create when a kmer spans over
				more than one node.
		"""
		remainingNumber = kmerLength - len(kmerPart)
		connectedNodes = list(self[currentNode[-1]])
		connectedNodeNumber = len(connectedNodes)

		for k in range(connectedNodeNumber):
			nextNode = self.nodes[connectedNodes[k]]
			nextNodeSeq = nextNode['sequence']
			nodeKeys = list(nextNode)
			nextNodeReverseSeq = nextNode['sequence'][::-1]
			# If the kmer can be fully created by the current node
			if len(nextNodeSeq) >= remainingNumber:
				vals = list(dict(nextNode).values())
				invs = False
				for num in vals:
					if isinstance(num, int):
						if num < 0:
							invs = True
				if invs == True:
					kmerSeq = kmerPart + nextNodeSeq[:remainingNumber]
					kmerRevSeq = kmerPart + nextNodeReverseSeq[:remainingNumber]
					allNodes = []
					pos = 0
					for x in range(len(currentNode)):
						allNodes.append(currentNode[pos])
						pos += 1
					allNodes.append(connectedNodes[k])
					tempDict = {'nodes': allNodes, 'sequence': kmerSeq, 'inversesequence': kmerRevSeq,
								'startposition': kmerStartPos}
					kmerKey = 'kmer_' + str(kmerNumber) + '_' + str(k + 1)
					if kmerKey in GgDiGraph.kmerDict.keys():
						lastDigit = kmerKey[-1]
						lastDigit = int(lastDigit) + 1
						newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
						kmerKey = newkmerKey
					newkmer = {kmerKey: tempDict}
					GgDiGraph.kmerDict.update(newkmer)
				else:
					# If the node is not longenough to make a full kmer you will need to recursively move onto its neighbours until
					# all of the kmers have been successfully created
					if hasInv == True:
						kmerSeq = kmerPart + nextNodeSeq[:remainingNumber]
						revKmerSeq = self.nodes[currentNode[-1]]['sequence']
						revKmerSeq = revKmerSeq[0:len(kmerPart)]
						revKmerSeq = revKmerSeq[::-1]
						revKmerSeq = revKmerSeq + nextNodeSeq[:remainingNumber]
						allNodes = []
						pos = 0
						for x in range(len(currentNode)):
							allNodes.append(currentNode[pos])
							pos += 1
						allNodes.append(connectedNodes[k])
						tempDict = {'nodes': allNodes, 'sequence': kmerSeq, 'inversesequence': revKmerSeq,
									'startposition': kmerStartPos}
						kmerKey = 'kmer_' + str(kmerNumber) + '_' + str(k + 1)
						if kmerKey in GgDiGraph.kmerDict.keys():
							lastDigit = kmerKey[-1]
							lastDigit = int(lastDigit) + 1
							newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
							kmerKey = newkmerKey
						newkmer = {kmerKey: tempDict}
						GgDiGraph.kmerDict.update(newkmer)
					else:
						kmerSeq = kmerPart + nextNodeSeq[:remainingNumber]
						allNodes = []
						pos = 0
						for x in range(len(currentNode)):
							allNodes.append(currentNode[pos])
							pos += 1
						allNodes.append(connectedNodes[k])
						tempDict = {'nodes': allNodes, 'sequence': kmerSeq, 'startposition': kmerStartPos}
						kmerKey = 'kmer_' + str(kmerNumber) + '_' + str(k + 1)
						if kmerKey in GgDiGraph.kmerDict.keys():
							lastDigit = kmerKey[-1]
							lastDigit = int(lastDigit) + 1
							newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
							kmerKey = newkmerKey
						newkmer = {kmerKey: tempDict}
						GgDiGraph.kmerDict.update(newkmer)
			else:
				kmerTemp = kmerPart + nextNodeSeq
				nodesInvolved = []
				posx = 0
				for x in range(len(currentNode)):
					nodesInvolved.append(currentNode[posx])
					posx += 1
				nodesInvolved.append(connectedNodes[k])
				# print(kmerTemp + " " +str(nodesInvolved))
				self.kmer_over_multiple_nodes(kmerTemp, nodesInvolved, kmerLength, kmerNumber, kmerStartPos, False)

	# recursively goes over multiple nodes
	# will need the currentNode that is parsed at this step so instead you have the two previous nodes the kmer goes over so you have all three nodes saved into the dictionary

	def kmers_of_node(self, nodeName):
		"""

		:param nodeName: the Node that you want to use create kmers from
		:return: returns a dictionary of all of the kmers involved with the given node.
		"""
		kmersOfNode = {}
		if len(GgDiGraph.kmerDict) == 0:
			# default kmer length is 4
			GgDiGraph.fast_kmer_create(self, 4)

		keyList = list(GgDiGraph.kmerDict.keys())
		# creates the kmers of the current node if the dictionary of all of the kmers is not already made
		for i in keyList:
			currentKmer = GgDiGraph.kmerDict.get(i)
			kmerNodeList = currentKmer['nodes']
			for j in kmerNodeList:
				if j == nodeName:
					kmer = {i: currentKmer}
					kmersOfNode.update(kmer)
		return kmersOfNode

	def create_query_kmers(self, querySequence, kmerLength):
		"""

		:param querySequence: the sequence of the query you want to align to the graph
		:param kmerLength: The length of the kmers that you want to make from the query sequence.
				These should be the same as the graph kmer length
		:return:Returns a dictionary with all of the kmers from the query sequence. The dictionary contains the query
				kmer name and query kmer sequence.
		"""
		startPos = 0
		queryKmers = {}
		kmerNumber = 1
		for i in range(len(querySequence)):
			if startPos + kmerLength <= len(querySequence):
				kmerSeq = querySequence[startPos:startPos + kmerLength]
				kmer = {'qkmer_' + str(kmerNumber): kmerSeq}
				kmerNumber += 1
				queryKmers.update(kmer)
				startPos += 1
		return queryKmers

	# need to fix readgaps. It is sometimes still not giving the correct values for gaps in between
	def readGaps(self, nodeNeighbours, query, firstReadInfo, secondReadInfo, distanceBetween):
		'''
		This function is not actually used as a more efficient way was found that bypasses this function. In any
		case this function should not be called and was called by the read aligner.
		:param nodeNeighbours: neighbours of the current node
		:param query: query sequence
		:param firstReadInfo: first aligned read dictionary information
		:param secondReadInfo: second aligned read dictionary information
		:param distanceBetween: the current distance between the two aligned reads
		:return: Returns the total distance between the two aligned reads
		'''
		global dist
		for x in nodeNeighbours:
			if x == secondReadInfo['nodescoveredbyread'][0]:
				if self.firstPath == True:
					self.dist = distanceBetween + (
							(len(self.nodes[firstReadInfo['nodescoveredbyread'][-1]]['sequence'])) - (
							firstReadInfo['alignementendpositioninlastnode'] + 1)) + secondReadInfo[
									'alignmentstartpositioninfirstnode']
					self.firstPath = False
					# print(self.dist)
					return self.dist
				else:
					return self.dist
		for i in nodeNeighbours:
			if len(self.nodes[i]['sequence']) + int(distanceBetween) < len(query):
				distanceBetween = int(distanceBetween) + len(self.nodes[i]['sequence'])
				newNeighbours = list(self[i])
				self.readGaps(newNeighbours, query, firstReadInfo, secondReadInfo, distanceBetween)
			else:
				return distanceBetween
		return distanceBetween

	def debruin_read_alignment(self, queryseq, kmerLength, outPutFilename):
		"""
		Running this Function will print out blocks of information on the aligned reads. Each block of 7 lines will correspond to a single aligned read.
		Author: Campbell Green
		Line 1: The Aligned read sequence
		Line 2: The Nodes that the aligned read covers
		Line 3: All nodes that the read aligns to inversely(Nodes with an aligned inversion)
		Line 4: The starting alignment position on the first node involved in the alignment
		Line 5: The ending alignment position on the last node involved in the alignment
		Line 6: Percentage of the initial query that has aligned to the graph
		Line 7:Percentage of the aligned read that has successfully aligned to the graph
		:param query: This is the read, as a string saved in a .txt file, that you want to try and align to the graph.
		:param kmerLength: The length of the kmers that will be created
		:return: returns a dictionary with the sequences that align, x's represent sections of the sequence that doesnt correctly align
				It also returns the nodes that the sequence aligns to and the position on the first node the sequence starts aligning to and the position on the last node that the sequence ends aligning to.
				These positions are in string positions, so the first base pair is position 0.
		"""
		queryKmerDict = self.create_query_kmers(queryseq, kmerLength)
		referenceKmerDict = self.fast_kmer_create(kmerLength)
		# creates the query and graph kmers with the desired kmer length
		finalkmerGroups = []
		groupedListOftwos = []
		matchedKmers = []
		finalAlignedReads = {}
		inversionNodes = []
		# list of lists of matched kmers
		# first pos is query kmer name, second pos is reference nodes that the kmer covers, third is the start position of reference kmer, fourth is reference kmer name, 5th is qkmerseq

		queryKmerKeys = list(queryKmerDict.keys())
		referenceKmerKeys = list(referenceKmerDict.keys())

		# mathces the each query kmer to refernce graph kmers with the same sequence
		for i in queryKmerKeys:
			currentQueryKmer = queryKmerDict.get(i)
			for j in referenceKmerKeys:
				currentReferenceKmer = referenceKmerDict.get(j)
				if currentReferenceKmer['sequence'] == currentQueryKmer:
					tempList = [i, currentReferenceKmer['nodes'], currentReferenceKmer['startposition'], j,
								currentQueryKmer, 'normal']
					matchedKmers.append(tempList)
				elif 'inversesequence' in currentReferenceKmer:
					if currentReferenceKmer['inversesequence'] == currentQueryKmer:
						tempList = [i, currentReferenceKmer['nodes'], currentReferenceKmer['startposition'], j,
									currentQueryKmer, 'inverse']
						matchedKmers.append(tempList)
		# results in a list of list with the matched kmers
		# has each match between query kmer and reference graph kmer

		# runs through the matched kmers and groups the kmers into groups of twos by the kmers that line up next to each other
		for k in matchedKmers:
			for x in matchedKmers:
				tempGroupedList = []
				currentQKmer = str(k[0])
				underscorePos = currentQKmer.find('_')
				qkmerNumber = currentQKmer[int(underscorePos) + 1:]
				nextQkmerNumber = int(qkmerNumber) + 1
				nextQkmer = 'qkmer_' + str(nextQkmerNumber)

				if nextQkmer == x[0]:
					if k[1][0] == x[1][0] and k[2] + 1 == x[2]:
						addOne = str(k[0]) + '-' + str(k[3])
						addTwo = str(x[0]) + '-' + str(x[3])
						if len(k[1]) == 1 and k[5] == 'inverse':
							node = referenceKmerDict[addOne[addOne.find('-') + 1:]]['nodes'][0]
							if node not in inversionNodes:
								inversionNodes.append(node)
							addOne = 'i' + addOne
						if len(x[1]) == 1 and x[5] == 'inverse':
							addTwo = 'i' + addTwo
							node = referenceKmerDict[addTwo[addTwo.find('-') + 1:]]['nodes'][0]
							if node not in inversionNodes:
								inversionNodes.append(node)
						tempGroupedList.append(addOne)
						tempGroupedList.append(addTwo)
						groupedListOftwos.append(tempGroupedList)
					elif len(self.nodes[k[1][0]]['sequence']) - 1 == k[2]:
						if len(k[1]) > 1:
							if k[1][1] == x[1][0] and x[2] == 0:
								addOne = str(k[0]) + '-' + str(k[3])
								addTwo = str(x[0]) + '-' + str(x[3])
								if (len(k[1]) == 1 and k[5] == 'inverse'):
									node = referenceKmerDict[addOne[addOne.find('-') + 1:]]['nodes'][0]
									if node not in inversionNodes:
										inversionNodes.append(node)
									addOne = 'i' + addOne
								if (len(x[1]) == 1 and x[5] == 'inverse'):
									node = referenceKmerDict[addTwo[addTwo.find('-') + 1:]]['nodes'][0]
									if node not in inversionNodes:
										inversionNodes.append(node)
									addTwo = 'i' + addTwo
								tempGroupedList.append(addOne)
								tempGroupedList.append(addTwo)
								groupedListOftwos.append(tempGroupedList)

		# end up with list of lists with groups of two kmers that are next to each other

		# takes the groups of two lists and groups all o the kmers that are in a line and next to each other on the graph
		for t in groupedListOftwos:
			continueCheck = True
			tempFinal = []

			for c in finalkmerGroups:
				if t[0] in c and t[1] in c:
					continueCheck = False

			if continueCheck == True:
				for l in t:
					tempFinal.append(l)

				for p in groupedListOftwos:
					if tempFinal[-1] == p[0]:
						tempFinal = tempFinal + list(set(p) - set(tempFinal))
				finalkmerGroups.append(tempFinal)
		# result in final groups with all kmers that are next to each other on the graphs grouped up
		AlignNumber = 1

		# creates a dictionary with all the final aligned kmers that align 100% perfectly to the reference graph
		for n in finalkmerGroups:
			matchedSquence = ''
			nodesCovered = []
			startpos = 0
			endpos = 0
			finalKmerName = n[-1][n[-1].find('-') + 1:]
			finalNode = referenceKmerDict[finalKmerName]['nodes'][-1]
			for m in n:
				inverse = False
				if m[0] == 'i':
					inverse = True
				kmerName = m[m.find('-') + 1:]
				# get kmerseq from matchedkmers not from the reference kmer list
				kmerseq = ''
				if inverse == True:
					qkmerName = m[1:m.find('-')]
				else:
					qkmerName = m[:m.find('-')]
				for r in matchedKmers:
					# shouldnt be finalkmerName, got to get each kmerName
					if r[0] == qkmerName and kmerName == r[3]:
						kmerseq = r[4]
				kmernode = referenceKmerDict[kmerName]['nodes']
				nodesCovered.extend(kmernode)
				positionsLeft = 0
				nodesCovered = list(dict.fromkeys(nodesCovered))
				if m == n[0]:
					startpos = referenceKmerDict[kmerName]['startposition']
				if m == n[-1]:
					if len(referenceKmerDict[kmerName]['nodes']) == 1 and referenceKmerDict[kmerName]['nodes'][
						-1] == finalNode:
						endpos = referenceKmerDict[kmerName]['startposition'] + kmerLength - 1
					else:
						for b in (referenceKmerDict[kmerName]['nodes']):
							if b == referenceKmerDict[kmerName]['nodes'][0]:
								tempnumber = len(self.nodes[b]['sequence']) - referenceKmerDict[kmerName][
									'startposition']
								positionsLeft += tempnumber
							elif b != referenceKmerDict[kmerName]['nodes'][0] and b != \
									referenceKmerDict[kmerName]['nodes'][-1]:
								positionsLeft += len(self.nodes[b]['sequence'])
						endpos = len(referenceKmerDict[kmerName]['sequence']) - positionsLeft - 1
				if len(matchedSquence) == 0:
					matchedSquence = matchedSquence + kmerseq
				else:
					matchedSquence = matchedSquence + kmerseq[-1]
			match = False
			for s in finalAlignedReads:
				if matchedSquence in str(finalAlignedReads[s]['sequence']) and set(nodesCovered).issubset(
						finalAlignedReads[s]['nodescoveredbyread']):
					match = True
			if match == False:
				tempValues = {'sequence': matchedSquence, 'nodescoveredbyread': nodesCovered,
							  'alignmentstartpositioninfirstnode': startpos, 'alignementendpositioninlastnode': endpos}
				tempAlign = {'Aligned_Read_' + str(AlignNumber): tempValues}
				AlignNumber += 1
				finalAlignedReads.update(tempAlign)
		# At the moment this returns reads that matched that are any size larger than one kmer length.
		# finalAlignedReads Are all reads that map 100 percent accurately
		# Start grouping the final aligned read dictionary to show bigger reads and not only 100% alignment
		finalGroupedReads = {}
		readGroups = []
		unlinkedReads = []

		# then finds the distances between the 100 percent mapped reads and if those distances are small enough it groups those 100 percent mapped reads together
		# this represents the final aligned reads which may not be 100 percent accurately mapped to each other.
		for z in finalAlignedReads:
			endNode = finalAlignedReads[z]['nodescoveredbyread'][-1]
			links = 0
			for c in finalAlignedReads:
				distanceBetween = 0
				startNode = finalAlignedReads[c]['nodescoveredbyread'][0]
				if z != c and list(finalAlignedReads.keys()).index(z) < list(finalAlignedReads.keys()).index(c):
					if int(z[13:]) + 1 == int(c[13:]):
						if endNode == startNode:
							distanceBetween = finalAlignedReads[c]['alignmentstartpositioninfirstnode'] - 1 - \
											  finalAlignedReads[z]['alignementendpositioninlastnode']
						else:
							if nx.has_path(self, endNode, startNode):
								path = nx.all_shortest_paths(self, endNode, startNode)
								pathList = list(path)
								for v in range(len(pathList)):
									tdistanceBetween = 0
									pointer = v
									startDist = 0
									for x in pathList[v]:
										if x == pathList[v][0]:
											startDist = len(self.nodes[x]['sequence']) - finalAlignedReads[z][
												'alignementendpositioninlastnode'] - 1
											tdistanceBetween = tdistanceBetween + startDist
										elif x == pathList[v][-1]:
											tdistanceBetween = tdistanceBetween + finalAlignedReads[c][
												'alignmentstartpositioninfirstnode']
										else:
											tdistanceBetween = tdistanceBetween + len(self.nodes[x]['sequence'])
									if distanceBetween == 0:
										distanceBetween = tdistanceBetween
									elif tdistanceBetween < distanceBetween and distanceBetween != 0:
										distanceBetween = tdistanceBetween
						if distanceBetween < len(queryseq) and distanceBetween != 0:
							readGroups.append([z, distanceBetween, c])
							links += 1
			if links == 0:
				contains = False
				for o in readGroups:
					if o[0] == z or o[2] == z:
						contains = True
				if contains == False:
					unlinkedReads.append(z)

		# Alignement read groups will shop the start and end point of the reads but will hva egaps in between that can be larger then gaps between the sections in the reads
		finalReadGroups = []
		intfinalReadGroups = []
		tempReadGroups = []
		for u in readGroups:
			temp = [u[0], u[-1]]
			tempReadGroups.append(temp)
		# print(tempReadGroups)
		out = []
		while len(tempReadGroups) > 0:
			# first, *rest = tempReadGroups
			first = tempReadGroups[0]
			rest = tempReadGroups[1:]
			first = set(first)

			lf = -1
			while len(first) > lf:
				lf = len(first)

				rest2 = []
				for r in rest:
					if len(first.intersection(set(r))) > 0:
						first |= set(r)
					else:
						rest2.append(r)
				rest = rest2

			out.append(first)
			tempReadGroups = rest
		for h in out:
			readNo = []
			intermediate = []
			for y in h:
				readNo.append(int(y[str(y).rfind('_') + 1:]))
			readNo.sort()
			for o in readNo:
				for p in h:
					if o == int(p[str(p).rfind('_') + 1:]):
						intermediate.append(p)
			intfinalReadGroups.append(intermediate)
		# print(intfinalReadGroups)

		for x in intfinalReadGroups:
			readsWithDist = []
			for c in range(len(x) - 1):
				for z in readGroups:
					if z[0] == x[c] and z[-1] == x[c + 1]:
						if c == 0:
							readsWithDist.append(x[c])
							readsWithDist.append(z[1])
							readsWithDist.append(x[c + 1])
						else:
							readsWithDist.append(z[1])
							readsWithDist.append(x[c + 1])
			finalReadGroups.append(readsWithDist)

		# all the final reads are now grouped together with the distances between the reads represented in the list of reads
		# Now the dictionary of final reads are created with the sequence of the read,the start and end positions of the reads
		# and the nodes that the reads covers
		# it also shows the percentage of the read that aligns to the graph and the percentage of the read that has aligned to the graph
		allAlignedReads = {}
		AlignNo = 1
		for z in finalReadGroups:
			sequence = ''
			checkPercentSeq = ''
			xCounter = 0
			nodesCoveredByRead = []
			startAlignPos = 0
			endAlignPos = 0
			percentOfQuery = 0
			percentAligned = 0
			for i in range(len(z)):
				if i == 0:
					readInfo = finalAlignedReads[z[i]]
					startAlignPos = readInfo['alignmentstartpositioninfirstnode']
					sequence = sequence + str(readInfo['sequence'])
					checkPercentSeq = checkPercentSeq + str(readInfo['sequence'])
					nodesCoveredByRead.extend(x for x in readInfo['nodescoveredbyread'] if x not in nodesCoveredByRead)
				elif i == len(z) - 1:
					readInfo = finalAlignedReads[z[i]]
					endAlignPos = readInfo['alignementendpositioninlastnode']
					sequence = sequence + str(readInfo['sequence'])
					checkPercentSeq = checkPercentSeq + str(readInfo['sequence'])
					nodesCoveredByRead.extend(x for x in readInfo['nodescoveredbyread'] if x not in nodesCoveredByRead)
				elif i % 2 == 0:
					readInfo = finalAlignedReads[z[i]]
					sequence = sequence + str(readInfo['sequence'])
					checkPercentSeq = checkPercentSeq + str(readInfo['sequence'])
					nodesCoveredByRead.extend(x for x in readInfo['nodescoveredbyread'] if x not in nodesCoveredByRead)
				else:
					if int(z[i]) < 0:
						endremove = int(z[i]) - 1
						sequence = sequence[:len(sequence) + endremove]
					else:
						for v in range(int(z[i])):
							xCounter += 1
							sequence = sequence + 'x'
			nodeNos = []
			newNodesCoveredByRead = []
			for r in nodesCoveredByRead:
				nodeNos.append(str(r[str(r).find('_') + 1:str(r).rfind('_')]) + str(r[str(r).rfind('_') + 1:]))
			for a in nodeNos:
				for b in nodesCoveredByRead:
					if a == str(b[str(b).find('_') + 1:str(b).rfind('_')]) + str(b[str(b).rfind('_') + 1:]):
						newNodesCoveredByRead.append(b)
			startAlignPos += 1
			endAlignPos += 1
			if newNodesCoveredByRead[-1] in inversionNodes:
				endAlignPos = endAlignPos * -1
			if newNodesCoveredByRead[0] in inversionNodes:
				startAlignPos = startAlignPos * -1
			percentOfQuery = difflib.SequenceMatcher(None, queryseq, checkPercentSeq).ratio() * 100
			percentAligned = (len(sequence) - xCounter) / len(sequence) * 100
			tempVal = {'sequence': sequence, 'nodescoveredbyread': newNodesCoveredByRead,
					   'alignmentstartpositioninfirstnode': startAlignPos,
					   'alignementendpositioninlastnode': endAlignPos,
					   'percentageofqueryaligned': percentOfQuery, 'percentageofalignedreadtograph': percentAligned}
			key = {'Alignment_' + str(AlignNo): tempVal}
			AlignNo += 1
			allAlignedReads.update(key)

		for w in unlinkedReads:
			read = finalAlignedReads.get(w)
			percentOfQuery = difflib.SequenceMatcher(None, queryseq, read['sequence']).ratio() * 100
			percentAligned = (len(read['sequence'])) / len(read['sequence']) * 100
			s = read['alignmentstartpositioninfirstnode']
			e = read['alignementendpositioninlastnode']
			s += 1
			e += 1
			if read['nodescoveredbyread'][-1] in inversionNodes:
				e = e * -1
			if read['nodescoveredbyread'][0] in inversionNodes:
				s = s * -1
			tempVals = {'sequence': read['sequence'], 'nodescoveredbyread': read['nodescoveredbyread'],
						'alignmentstartpositioninfirstnode': s, 'alignementendpositioninlastnode': e,
						'percentageofqueryaligned': percentOfQuery, 'percentageofalignedreadtograph': percentAligned}
			key = {'Alignment_' + str(AlignNo): tempVals}
			AlignNo += 1
			allAlignedReads.update(key)

		# print out the blocks of information for each aligned read
		print('7 line block format for each aligned read')
		print('Line 1: The Aligned read sequence ')
		print('Line 2: The Nodes that the aligned read covers')
		print('Line 3: All nodes that the read aligns to inversely(Nodes with an aligned inversion)')
		print('Line 4: The starting alignment position on the first node involved in the alignment')
		print('Line 5: The ending alignment position on the last node involved in the alignment')
		print('Line 6: Percentage of the initial query that has aligned to the graph')
		print('Line 7:Percentage of the aligned read that has successfully aligned to the graph')
		print()

		for t in list(allAlignedReads.keys()):
			print(allAlignedReads[t]['sequence'])
			print(allAlignedReads[t]['nodescoveredbyread'])
			if len(inversionNodes) == 0:
				print('None')
			else:
				currentInversionNodes = []
				for x in allAlignedReads[t]['nodescoveredbyread']:
					if x in inversionNodes:
						currentInversionNodes.append(x)
				print(currentInversionNodes)
			print(allAlignedReads[t]['alignmentstartpositioninfirstnode'])
			print(allAlignedReads[t]['alignementendpositioninlastnode'])
			print(allAlignedReads[t]['percentageofqueryaligned'])
			print(allAlignedReads[t]['percentageofalignedreadtograph'])
			print()
		self.create_alignment_JSON(allAlignedReads,inversionNodes,outPutFilename)
		return allAlignedReads

	def create_alignment_JSON(self, AlignedReads, inversionNodes, fileName):
		file = open(fileName + '.JSON','w+')
		file.write('[\n')
		for i in list(dict(AlignedReads)):
			file.write('\t{ \n')
			file.write('\t\t"alignmentId": "'+i+'",\n')
			file.write('\t\t"sequence": "' + str(AlignedReads[i]['sequence']) + '",\n')
			file.write('\t\t"alignedNodes": "' + str(AlignedReads[i]['nodescoveredbyread']) + '",\n')
			currentInversionNodes = []
			for x in AlignedReads[i]['nodescoveredbyread']:
				if x in inversionNodes:
					currentInversionNodes.append(x)
			file.write('\t\t"invertedAlignedNodes": "' + str(currentInversionNodes) + '",\n')
			file.write('\t\t"alignmentStartPosition": ' + str(AlignedReads[i]['alignmentstartpositioninfirstnode']) + ',\n')
			file.write('\t\t"alignmentEndPosition": ' + str(AlignedReads[i]['alignementendpositioninlastnode']) + ',\n')
			file.write('\t\t"queryAlignedPercent": "' + str(AlignedReads[i]['percentageofqueryaligned']) + '",\n')
			file.write('\t\t"readToGraphPercent": "' + str(AlignedReads[i]['percentageofalignedreadtograph']) + '",\n')
			file.write('\t}\n')
		file.write(']')

	def extract_JSON(self,JSONfile):
		print('7 line block format for each aligned read')
		print('Line 1: The Aligned read sequence ')
		print('Line 2: The Nodes that the aligned read covers')
		print('Line 3: All nodes that the read aligns to inversely(Nodes with an aligned inversion)')
		print('Line 4: The starting alignment position on the first node involved in the alignment')
		print('Line 5: The ending alignment position on the last node involved in the alignment')
		print('Line 6: Percentage of the initial query that has aligned to the graph')
		print('Line 7:Percentage of the aligned read that has successfully aligned to the graph')
		print()

		AlignedRead = {}
		file = open(JSONfile,'r')
		file.readline()
		checkLine = file.readline()
		while checkLine != ']':
			id = file.readline()
			seq = file.readline()
			nodes = file.readline()
			invert = file.readline()
			start = file.readline()
			end = file.readline()
			qPercent = file.readline()
			rPercent = file.readline()
			file.readline()
			checkLine = file.readline()

			tempVals = {'sequence': seq[int(seq.find(':'))+3:-3], 'nodescoveredbyread': nodes[nodes.find(':')+3: -3],
						'alignmentstartpositioninfirstnode': start[start.find(':')+2: -2], 'alignementendpositioninlastnode': end[end.find(':')+2: -2],
						'percentageofqueryaligned': qPercent[qPercent.find(':')+3:-3], 'percentageofalignedreadtograph': rPercent[rPercent.find(':')+3: -3]}
			key = {id[id.find(':')+3:-3]: tempVals}
			AlignedRead.update(key)

			print(seq[seq.find(':')+3:-3])
			print(nodes[nodes.find(':')+3:-3])
			print(invert[invert.find(':')+3:-3])
			print(start[start.find(':')+2:-2])
			print(end[end.find(':')+2:-2])
			print(qPercent[qPercent.find(':')+3:-3])
			print(rPercent[rPercent.find(':')+3:-3])
			print()

		return AlignedRead


# ---------------------------------------------------- New functions under testing


def import_gg_graph(path):
	"""
	Import a GG graph
	:param path: file path
	:return: a GG graph genome object
	"""
	graph_obj = nx.read_graphml(path)

	out_graph = GgDiGraph(graph_obj)

	return out_graph


def get_neighbours_context(graph, source_node, label, dir='out'):
	'''Retrieves neighbours of a node using only edges containing a certain label'''
	'''Created by shandu mulaudzi'''

	# Changes sequence to ids
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


def extract_region_subgraph(graph, region_start, region_stop, seq_name, neighbours=0):

	"""
	Extracts the sub-graph of a region between two coordinates for a sequence ID. An example may be
	returning the sub-graph that represents a gene for a particular isolate. Specifying a value for
	neighbours will result in the graph plus any adjacent nodes to be returned.
	:param graph: A input gengraph object
	:param region_start: The start position (bp) of the region.
	:param region_stop: The stop position (bp) of the region.
	:param seq_name: The sequence ID that the base positions are referencing.
	:param neighbours: Whether to include neighbouring nodes. This currently only considers adjacent nodes, and
	will be expanded to include neighbours of neighbours etc.
	:return: A subgraph as a gengraph object.
	"""


	node_list = []
	pre_node_list = []

	for node, data in graph.nodes(data=True):

		if seq_name in data['ids'].split(','):

			node_start = abs(int(data[seq_name + '_leftend']))
			node_stop = abs(int(data[seq_name + '_rightend']))

			# Check node orientation for sequence
			if int(data[seq_name + '_leftend']) < 0:
				orientation = '-'
			else:
				orientation = '+'

			# Check if whole sequence found in one node
			if region_start > node_start and region_stop < node_stop:

				if neighbours != 0:
					node_list.append(node)
					for a_neighbour in graph.neighbors(node):
						node_list.append(a_neighbour)
						print(a_neighbour)
					for a_neighbour in graph.predecessors(node):
						print(a_neighbour)
						node_list.append(a_neighbour)
				else:
					node_list = [node]

				sub_graph = graph.subgraph(node_list)

				return sub_graph

			# For the middle nodes (nodes containing a whole sub-sequence)
			elif region_start < node_start and region_stop > node_stop:
				pre_node_list.append(node)

			# For the start node
			elif region_start < node_stop < region_stop:
				pre_node_list.append(node)

			# For the end node
			elif region_start < node_start < region_stop:
				pre_node_list.append(node)

	# This could be better, the neighbors function also returns the query node.
	if neighbours != 0:
		for a_hit_node in pre_node_list:
			node_list.append(a_hit_node)
			for a_neighbour in graph.neighbors(a_hit_node):
				node_list.append(a_neighbour)
			for a_neighbour in graph.predecessors(a_hit_node):
				node_list.append(a_neighbour)
	else:
		node_list = pre_node_list

	sub_graph = graph.subgraph(node_list)

	return sub_graph


# ---------------------------------------------------- General functions 

def input_file_check(input_dict):
	"""
	Checkto see if the input test file is formatted correctly, and the fasta files will work.
	:param input_dict: The input_dict returned from input_parser
	:return: A list of any errors
	"""
	# TODO: Catch invisible chars

	errors_list = []

	# Check if the header is correct

	# Check if fasta files have multi chromosomes
	for isolate, a_fasta_path in input_dict[1].items():

		# Check if files are where they should be
		if os.path.isfile(a_fasta_path) is False:
			errors_list.append('Fasta file not found - ' + a_fasta_path)
		else:
			file = open(a_fasta_path, 'r').read()
			chrom_count = file.count('>')

			if chrom_count != 1:
				errors_list.append('Input file fail - ' + str(chrom_count) + ' chromosomes seen in fasta file: ' + a_fasta_path)

	return errors_list


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
		
		gene_ID_dict = {"gene_details": sequence_details[1:], "DNA_seq" : sequence}

		output_list.append(gene_ID_dict)

		return output_list

	if file_path[-4:] == ".bim":
		input_file = open(file_path, "r")
		return input_file
		
	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result = numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".ssv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=';')
		data_matrix = list(data_table)
		result = numpy.array(data_matrix)
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
		for node, data in a_graph.nodes(data=True):
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
				new_node_dict = {isolate + '_leftend':sorted_list[count][2] + 1, isolate + '_rightend':sorted_list[count + 1][1] - 1, 'ids':isolate}
				new_node_dict['sequence'] = isolate_Seq[sorted_list[count][2]:sorted_list[count + 1][1] - 1]

				node_ID = isolate + "_" + str(count)
				att_node_dict = {node_ID: new_node_dict}

				a_graph.add_node(node_ID)
				nx.set_node_attributes(a_graph, att_node_dict)

			count += 1


def node_check(a_graph):
	from operator import itemgetter
	logging.info('checking nodes')
	doespass = True

	iso_list = a_graph.graph['isolates'].split(',')

	for isolate in iso_list:
		isolate_node_list = []
		for node, data in a_graph.nodes(data=True):
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
		for node, data in a_graph.nodes(data=True):
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

	# Because MAUVE does not use the sequence names in the bbone file, you need to replace the generic names like 'seq0'
	# with the identifiers like 'H37Rv'
	for header_item in header_line:
		for seqID in input_dict[0].keys():
			if seqID == header_item.split('_')[0]:
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

		found_in_string = ','.join(found_in_list)

		node_dict['ids'] = found_in_string

		node_ID = 'Aln_' + str(node_count)

		node_dict['name'] = node_ID

		# Error here from networkx
		att_node_dict = {node_ID: node_dict}
		genome_network.add_node(node_ID)
		nx.set_node_attributes(genome_network, att_node_dict)

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
					genome_network.add_node(node_ID, **node_dict)
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
	for node, data in inGraph.nodes(data=True):
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

	for node, data in graph_obj.nodes(data=True):

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

	all_node_list = sorted(pos_lol, key=lambda start_val: start_val[0])


	# -------------------------------------------------------- Here we add edges to the graph, weaving in the new nodes

	count = 0

	edges_obj = graph_obj.edges()

	while count < (len(all_node_list) - 1):

		node_1 = all_node_list[count][2]
		node_2 = all_node_list[count+1][2]

		if (node_1, node_2) in edges_obj:

			print(node_1)
			print(node_2)
			nx.write_graphml(graph_obj, 'problemG')
			print(graph_obj.get_edge_data(node_1, node_2))
			print(graph_obj.get_edge_data(node_1, node_2)[0]['ids'])
			print(graph_obj.get_edge_data(node_1, node_2)[0]['ids'].split(','))

			if sequence_name not in graph_obj.get_edge_data(node_1, node_2)[0]['ids'].split(','):

				new_seq_list = graph_obj.get_edge_data(node_1, node_2)[0]['ids'] + ',' + sequence_name

			else:
				print('Error found.')
				quit()
				new_seq_list = graph_obj.edges[node_1][node_2]['ids']

			nx.set_edge_attributes(graph_obj, {'ids': new_seq_list})

		else:
			graph_obj.add_edge(node_1, node_2, ids=sequence_name)

		count += 1

	logging.info('Edges added')

	return graph_obj


def link_all_nodes(graph_obj):
	print('Linking nodes')
	logging.info('Running link_all_nodes')
	for isolate in graph_obj.graph['isolates'].split(','):
		logging.info(isolate)
		graph_obj = link_nodes(graph_obj, isolate)
	return graph_obj


def add_sequences_to_graph(graph_obj, paths_dict):

	logging.info('Adding sequences')

	for node, data in graph_obj.nodes(data=True):

		seq_source = data['ids'].split(',')[0]
		is_reversed = False

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
				node_seq = reverse_compliment(node_seq)

			graph_obj.node[node]['sequence'] = node_seq

	return graph_obj


def add_sequences_to_graph_fastaObj(graph_obj, imported_fasta_object):

	logging.info('Adding sequences')

	seqObj = reshape_fastaObj(imported_fasta_object)

	for node, data in graph_obj.nodes(data=True):

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

	graph_isolate_list = graph_obj.graph['isolates'].split(',')

	expected_ref_length = int(stop_pos) - int(start_pos)
	expected_ref_length = expected_ref_length + 1

	if expected_ref_length < 0:
		print('POSSIBLE ERROR, start larger than stop')
		print(int(stop_pos), int(start_pos))
		print(expected_ref_length)

	# Extract subgraph
	node_leftend_label = reference_name + '_leftend'
	node_rightend_label = reference_name + '_rightend'

	# Identify the nodes that contain the start and stop positions

	print(reference_name)

	for node, data in graph_obj.nodes(data=True):
		if reference_name in data['ids'].split(','):

			if int(data[node_leftend_label]) > 0:
				#print 'positive'
				if abs(int(data[node_leftend_label])) <= int(start_pos) <= abs(int(data[node_rightend_label])):
					start_node = node
					print('Ping')
					print(start_pos)
					start_overlap = bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) <= int(stop_pos) <= abs(int(data[node_rightend_label])) or abs(int(data[node_leftend_label])) >= int(stop_pos) >= abs(int(data[node_rightend_label])):
					print('ping')
					stop_node = node
					print(stop_pos)
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])

			if int(data[node_leftend_label]) < 0:
				if abs(int(data[node_leftend_label])) <= int(start_pos) <= abs(int(data[node_rightend_label])):
					start_node = node
					print('Ping')
					print('neg node')
					print(start_pos)
					start_overlap =  bp_distance(data[node_rightend_label], start_pos)

				if abs(int(data[node_leftend_label])) <= int(stop_pos) <= abs(int(data[node_rightend_label])):
					logging.info('neg node')
					stop_node = node
					logging.info(stop_pos)
					stop_overlap = bp_distance(stop_pos, data[node_leftend_label])

	# Dealing with genes that occur at the start and stop nodes of the graph... needs a propper solution
	# Caused by the lack of a 'length' attribute in the start and stop node
	# Temp fix = just return a nope.

	if start_node.split('_')[-1] == 'start' or stop_node.split('_')[-1] == 'stop':
		if return_dict:
			# THIS WILL FAIL maybe...
			return {reference_name: 1}
		else:
			return [reference_name]

	# Calculating overlap

	start_node_pb_length = bp_length_node(graph_obj.node[start_node])
	stop_node_pb_length = bp_length_node(graph_obj.node[stop_node])

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

	for path_node in nodes_in_path:
		if reference_name in graph_obj.node[path_node]['ids']:
			ref_nodes_in_path.append(path_node)
		else:
			alt_nodes_in_path.append(path_node)

	total_alt_length = 0
	total_ref_length = 0

	for alt_node in alt_nodes_in_path:
		total_alt_length = total_alt_length + bp_length_node(graph_obj.node[alt_node])



	for ref_node in ref_nodes_in_path:
		if ref_node != start_node and ref_node != stop_node:
			total_ref_length = total_ref_length + bp_length_node(graph_obj.node[ref_node])

	total_ref_length = total_ref_length + start_overlap + stop_overlap

	iso_diff_dict = {}
	iso_sim_dict = {}

	for node_iso in graph_isolate_list:
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

	if return_dict:
		return iso_sim_score_dict
	else:
		return ret_list


def convert_coordinates(graph_obj, q_start, q_stop, ref_iso, query_iso):
	"""
	Convert coordinates from one sequence to another
	:param graph_obj:
	:param q_start:
	:param q_stop:
	:param ref_iso:
	:param query_iso:
	:return:
	"""

	# Find the right node

	node_leftend_label = ref_iso + '_leftend'
	node_rightend_label = ref_iso + '_rightend'

	query_leftend_label = query_iso + '_leftend'
	query_rightend_label = query_iso + '_rightend'

	ref_node_left = ''
	ref_node_right = ''

	query_node_left = ''
	query_node_right = ''

	for node, data in graph_obj.nodes(data=True):
		if ref_iso in data['ids'].split(',') and query_iso in data['ids'].split(','):
			if int(data[node_leftend_label]) < int(q_start) and int(q_stop) < int(data[node_rightend_label]):
				ref_node_left = data[node_leftend_label]
				ref_node_right = data[node_rightend_label]
				query_node_left = data[query_leftend_label]
				query_node_right = data[query_rightend_label]

				ref_low_num = abs(int(ref_node_left))
				if abs(int(ref_node_left)) > abs(int(ref_node_right)):
					ref_low_num = abs(int(ref_node_right))

				query_low_num = abs(int(query_node_left))
				if abs(int(query_node_left)) > abs(int(query_node_right)):
					query_low_num = abs(int(query_node_right))

				conversion_factor = ref_low_num - query_low_num

				new_r_start = int(q_start) - conversion_factor
				new_r_stop = int(q_stop) - conversion_factor

				return {query_iso + '_leftend':new_r_start, query_iso + '_rightend':new_r_stop}


def convert_coordinate(graph_obj, q_position, ref_iso, query_iso):

	# Find the right node

	node_leftend_label = ref_iso + '_leftend'
	node_rightend_label = ref_iso + '_rightend'

	query_leftend_label = query_iso + '_leftend'
	query_rightend_label = query_iso + '_rightend'


	for node, data in graph_obj.nodes(data=True):

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

	last_bline = {'blocklist': []}

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
	"""
	This used a multiple sequence aligner to realign the nodes that represent blocks of sequences.
	:param in_graph: input graph object that represents the blocks of co-linear sequences
	:param node_ID: The node to realign. String.
	:param seq_fasta_paths_dict: Dict containing the paths to the fasta sequences.
	:return:
	"""
	# TODO: This fails if there are multiple chromosomes in the input fasta file.
	# TODO: load fasta into memory instead of reading each time.

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

		''' 
		Currently only rev comp sequences are seen in the BBone file, represented by a - but not reversed start / stop 
		'''

		if int(node_data_dict[node_isolate + '_leftend']) > 0:
			orientation_dict[node_isolate] = '+'
		else:
			orientation_dict[node_isolate] = '-'

		if orientation_dict[node_isolate] == '+':
			iso_node_seq = iso_full_seq[int(node_data_dict[node_isolate + '_leftend']) - 1:int(node_data_dict[node_isolate + '_rightend']) ]
			node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])

		else:
			slice_start = abs(int(node_data_dict[node_isolate + '_leftend'])) - 1
			slice_end = abs(int(node_data_dict[node_isolate + '_rightend']))
			iso_node_seq = iso_full_seq[slice_start: slice_end]
			node_seq_start_pos[node_isolate] = int(node_data_dict[node_isolate + '_leftend'])
			iso_node_seq = reverse_compliment(iso_node_seq)

		temp_fasta_file.write('>' + node_isolate + '\n')
		temp_fasta_file.write(iso_node_seq + '\n')
		node_seq_len_est = len(iso_node_seq)

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
	for node, data in graph_obj.nodes(data=True):

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


def progressiveMauve_alignment(path_to_progressiveMauve, fasta_path_list, out_aln_name):
	"""
	A wrapper for progressiveMauve for use in GenGraph for the identification of co-linear blocks
	:param path_to_progressiveMauve: Absolute path to progressiveMauve executable
	:param fasta_path_list: List of paths to fasta files
	:param out_aln_name: Name for alignment file, added to mauve output
	:return:

	"""
	# Maybe add --skip-gapped-alignment flag?

	logging.info(path_to_progressiveMauve)
	progressiveMauve_call = [path_to_progressiveMauve, '--output=globalAlignment_' + out_aln_name, '--scratch-path-1=./mauveTemp', '--scratch-path-2=./mauveTemp'] + fasta_path_list

	return call(progressiveMauve_call, stdout=open(os.devnull, 'wb'))

# ---------------------------------------------------- Utility functions


def nodes_connected(u, v, graph_obj):
	return u in graph_obj.neighbors(v)


def convert_transcriptome_to_aln_input(trans_fasta_path, out_name):

	chrom_name = 'MTB_Pan1'

	break_length = 120

	fasta_lol = input_parser(trans_fasta_path)

	# Make masked break string

	break_str = 'N' * break_length

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

	for node, data in graph_obj.nodes(data=True):

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

	# REMOVE AND PUT AS TEST

	last_node_line = [0,0,0]

	for segment in pos_lol:

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
	Retrieve a sub-sequence from the graph for an isolate given the positions. This can be used for example to extract
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

	for node, data in graph_obj.nodes(data=True):

		if seq_name in data['ids'].split(','):

			node_start = abs(int(data[seq_name + '_leftend']))
			node_stop = abs(int(data[seq_name + '_rightend']))

			# Check node orientation for sequence
			if int(data[seq_name + '_leftend']) < 0:
				orientation = '-'
			else:
				orientation = '+'
			# Check if whole sequence found in one node
			if region_start >= node_start and region_stop <= node_stop:

				if orientation == '-':
					# Reverse comp seq stuff

					node_seq = reverse_compliment(data['sequence'])
					slce_start = region_start - node_start
					slce_end = region_stop - node_start
					return node_seq[slce_start:slce_end + 1]

				else:
					slce_start = region_start - node_start
					slce_end = region_stop - node_start

					return data['sequence'][slce_start:slce_end + 1]

			elif node_start <= region_start <= node_stop:

				#print('start node')

				if orientation == '-':
					# Reverse comp seq stuff
					rev_seq = reverse_compliment(data['sequence'])
					sub_seq = rev_seq[region_start - node_start:]
					multi_node_dict[node_start] = sub_seq

				else:
					sub_seq = data['sequence'][region_start - node_start:]
					multi_node_dict[node_start] = sub_seq

			elif node_start <= region_stop <= node_stop:

				#print('end node')

				if orientation == '-':
					# Reverse comp seq stuff
					rev_seq = reverse_compliment(data['sequence'])
					sub_seq = rev_seq[0: region_stop - node_start + 1]
					multi_node_dict[node_start] = sub_seq

				else:
					sub_seq = data['sequence'][0: region_stop - node_start + 1]
					multi_node_dict[node_start] = sub_seq

			elif region_start <= node_start and node_stop <= region_stop:

				#print('mid node')

				if orientation == '-':
					# Reverse comp seq stuff
					rev_seq = reverse_compliment(data['sequence'])
					multi_node_dict[node_start] = rev_seq

				else:

					multi_node_dict[node_start] = data['sequence']

	result_string = ''
	#print(multi_node_dict)
	for key in sorted(multi_node_dict.keys()):
		result_string += multi_node_dict[key]

	return result_string


def extract_iso_subgraph(graph_obj, isolate):
	iso_graph = nx.MultiDiGraph()
	iso_graph.graph['isolate'] = 'isolate'

	for node, data in graph_obj.nodes(data=True):
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

	for n, d in in_graph.nodes(data=True):
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


# ---------------------------------------------------- # read alignment
# Functions to align reads to the GenGraph genome graph.

# ------------------ Traditional approaches.

# Break down into k-mers to either create a hash table, or to create de-bruijn graphs.



