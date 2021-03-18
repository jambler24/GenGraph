import unittest
from gengraphTool import *
import networkx as nx
import matplotlib.pyplot as plt
from gengraph import *
import networkx.algorithms.isomorphism as iso


# - No mutations, all sequences the same, only one node
# - One mutation in one of the strains, a single substitution.
# - Same for a single deletion
# - And an insertion
# - A large rearrangement where 300bp are reverse complimented in one sequence.
# - A large rearrangement where 300bp are reverse complimented in one sequence containing a SNP in the reversed sequence.
# - One with a general mix of SNPS

# TO run a specific test: python3 -m unittest -v unitTest_Alignment.MultipleSequenceAlignmentTest.align_one_read
# To profile do: python3 -m cProfile -o gengraphTool_6genome.prof ./gengraphTool.py make_genome_graph --seq_file ./test_files/multiGenome6.txt --out_file_name test6genome

test_graph_dir = './TestGraphs/'

class BasicTests(unittest.TestCase):
    noMutationRead = "GCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGA"
    noMutationGraph = import_gg_graph(test_graph_dir + 'no_mutation.xml')

    def test_read_Length(self):
        self.assertEqual(readLength(self.noMutationRead), 100)

    def test_sequence_Length(self):
        self.assertEqual(sequenceLength(self.noMutationGraph), 1000)

    def test_no_of_Isolates(self):
        self.assertEqual(isolateNumber(self.noMutationGraph), 3)



class NoMutationsTests(unittest.TestCase):
    noMutationRead = "GCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGA"
    noMutationGraph = import_gg_graph(test_graph_dir + 'no_mutation.xml')

    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.noMutationGraph), 1)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.noMutationGraph, self.noMutationRead), 1)

    def test_alignment_Start(self):

        self.assertEqual(alignStartPosition(self.noMutationGraph, self.noMutationRead), 100)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.noMutationGraph, self.noMutationRead), 199)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(self.noMutationGraph, self.noMutationRead), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.noMutationGraph, self.noMutationRead), "100%")


class OneSubstitutionTest(unittest.TestCase):

    oneSubstitutionGraph = import_gg_graph(test_graph_dir + 'one_sub.xml')
    readWithSubstitution = "CACAGTGTGGCACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAG"
    readWithoutSubstitution = "CACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAG"

    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.oneSubstitutionGraph), 4)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.oneSubstitutionGraph, self.readWithoutSubstitution), 3)
        self.assertEqual(coverage(self.oneSubstitutionGraph, self.readWithSubstitution), 3)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.oneSubstitutionGraph, self.readWithoutSubstitution), 30)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.oneSubstitutionGraph, self.readWithoutSubstitution), 129)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(self.oneSubstitutionGraph, self.readWithoutSubstitution), "10%")
        self.assertEqual(coverage(self.oneSubstitutionGraph, self.readWithSubstitution), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.oneSubstitutionGraph, self.readWithoutSubstitution), "100%")
        self.assertEqual(queryAlignmentPercentage(self.oneSubstitutionGraph, self.readWithSubstitution), "100%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv'), self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv_2')), 0)
        self.assertEqual(isolateSequenceDifferences(self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv'), self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv_One_Sub')), 1)
        self.assertEqual(isolateSequenceDifferences(self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv_2'), self.oneSubstitutionGraph.get_sequence(0, 999, 'H37Rv_One_Sub')), 1)

class OneDeletionTest(unittest.TestCase):

    oneDeletionGraph = import_gg_graph(test_graph_dir + 'one_deletion.xml')
    readWithDeletion = "CACAGTGTGGAACGCGGTCGCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGC"
    readWithoutDeletion = "CACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAG"

    def test_no_of_nodes(self):
        self.assertEqual(nodeNumber(self.oneDeletionGraph), 3)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.oneDeletionGraph, self.readWithoutDeletion), 3)
        self.assertEqual(coverage(self.oneDeletionGraph, self.readWithDeletion), 3)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.oneDeletionGraph, self.readWithoutDeletion), 30)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.oneDeletionGraph, self.readWithoutDeletion), 129)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(self.oneDeletionGraph, self.readWithoutDeletiond), "10%")
        self.assertEqual(coverage(self.oneDeletionGraph, self.readWithDeletion), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.oneDeletionGraph, self.readWithoutDeletion), "100%")
        self.assertEqual(queryAlignmentPercentage(self.oneDeletionGraph, self.readWithDeletion), "10%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv'), self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv_2')), 0)
        self.assertEqual(isolateSequenceDifferences(self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv'), self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv_One_Deletion')), 1)
        self.assertEqual(isolateSequenceDifferences(self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv_2'), self.oneDeletionGraph.get_sequence(0, 999, 'H37Rv_One_Deletion')), 1)

class OneInsertionTest(unittest.TestCase):

    oneInsertionGraph = import_gg_graph(test_graph_dir + 'one_insertion.xml')
    readWithInsertion = "CACAGTGTGGAACGCGGTCGTCTCCGAACTTAGACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCA"
    readWithoutInsertion = "CACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAG"
    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.oneInsertionGraph), 3)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.oneInsertionGraph, self.readWithoutInsertion), 2)
        self.assertEqual(coverage(self.oneInsertionGraph, self.readWithInsertion), 3)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.oneInsertionGraph, self.readWithoutInsertion), 30)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.oneInsertionGraph, self.readWithoutInsertion), 129)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(self.oneInsertionGraph, self.readWithoutInsertion), "10%")
        self.assertEqual(coverage(self.oneInsertionGraph, self.readWithInsertion), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.oneInsertionGraph, self.readWithoutInsertion), "100%")
        self.assertEqual(queryAlignmentPercentage(self.oneInsertionGraph, self.readWithInsertion), "100%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv'), self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv_2')), 0)
        self.assertEqual(isolateSequenceDifferences(self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv'), self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv_One_Deletion')), 1)
        self.assertEqual(isolateSequenceDifferences(self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv_2'), self.oneInsertionGraph.get_sequence(0, 999, 'H37Rv_One_Deletion')), 1)

'''
class ReverseComplimentTest(unittest.TestCase):

    reverseComplimentGraph = import_gg_graph(test_graph_dir + 'reverse.xml')
    read = "TCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGGTTAGAGCAGCAACAGCCAACACCACAGACCGCTACACCATCGTCCTAAA"
    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.reverseComplimentGraph), 245)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.reverseComplimentGraph, self.read), 29)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.reverseComplimentGraph, self.read), 50)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.reverseComplimentGraph, self.read), 149)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(graphObject, read), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.reverseComplimentGraph, self.read), "100%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(seq1, seq2), 0)
        self.assertEqual(isolateSequenceDifferences(seq1, seq3), 300)
        self.assertEqual(isolateSequenceDifferences(seq2, seq3), 300)

    def test_contain_reverse(self):
        self.assertTrue(containReverse(self.reverseComplimentGraph))

class ReverseComplimentWithSNPTest(unittest.TestCase):


    reverseComplimentGraph = import_gg_graph(test_graph_dir + 'reverse_with_snp.xml')
    read = "TCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATTGTTAGAGCAGCAACAGCCAACACCACAGACCGCTACACCATCGTCCTAAA"
    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.reverseComplimentGraph), 245)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.reverseComplimentGraph, self.read), 29)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.reverseComplimentGraph, self.read), 50)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.reverseComplimentGraph, self.read), 149)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(graphObject, read), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.reverseComplimentGraph, self.read), "100%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(seq1, seq2), 0)
        self.assertEqual(isolateSequenceDifferences(seq1, seq3), 300)
        self.assertEqual(isolateSequenceDifferences(seq2, seq3), 300)

    def test_contain_reverse(self):
        self.assertTrue(containReverse(self.reverseComplimentGraph))
'''

class MixOfSNPTest(unittest.TestCase):

    mixOfSNPGraph = import_gg_graph(test_graph_dir +  'mix_of_snps.xml')
    read = "TTGCCCGATGACCCCGGTTCAGGCTTCACCAGAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATG"
    def test_no_of_nodes(self):

        self.assertEqual(nodeNumber(self.mixOfSNPGraph), 19)

    def test_nodes_covered_by_read(self):
        self.assertEqual(coverage(self.mixOfSNPGraph, self.read), 5)

    def test_alignment_Start(self):
        self.assertEqual(alignStartPosition(self.mixOfSNPGraph, self.read), 0)

    def test_alignment_End(self):
        self.assertEqual(alignEndPosition(self.mixOfSNPGraph, self.read), 99)

    def test_coverage_Percentage(self):
        self.assertEqual(coverage(self.mixOfSNPGraph, self.read), "10%")

    def test_query_alignment_Percentage(self):
        self.assertEqual(queryAlignmentPercentage(self.mixOfSNPGraph, self.read), "100%")

    def test_sequence_Differences(self):
        self.assertEqual(isolateSequenceDifferences(seq1, seq2), 0)
        self.assertEqual(isolateSequenceDifferences(seq1, seq3), 6)
        self.assertEqual(isolateSequenceDifferences(seq2, seq3), 6)

    def test_contain_reverse(self):
        self.assertFalse(containReverse(self.mixOfSNPGraph))

    def test_no_of_SNPS(self):
        self.assertEqual(noOfSNPS(self.mixOfSNPGraph), 6)

    def test_list_of_SNPS(self):
        self.assertEqual(listOfSNPS(self.mixOfSNPGraph), ["A,C,4", "C,G,32", "C,A,118", "A,C,166", "C,G,240", "A,T,993"])

class GraphSimilarityTest(unittest.TestCase):
    # Can remove section graph aligns to and insert the read, then check whether graphs are similar or the same
    oneSubstitutionGraph = import_gg_graph(test_graph_dir + 'one_sub.xml')
    readWithSubstitution = "CACAGTGTGGCACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAG"
    readwithOneDifference = "CACAGTGTGGCACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGTGCTCCGCTGACCCCTCAG"

    def test_isomorphism(self):
        # changedGraph has matched read inserted to replace sequence it matches to
        # if it is ismorphic then it proves that the read aligns perfectly to a section of the GenGraph
        changedGraph = insertAlignedRead(self.oneSubstitutionGraph, self.readWithSubstitution)
        self.assertTrue(isIsomorphic(self.oneSubstitutionGraph , changedGraph))
        self.assertEqual(calculateGED(self.oneSubstitutionGraph, changedGraph), 0.0)

    def test_GED(self):
        # The graph edit distance is the number of changes to nodes/edges to make two graphs isomorphic(identical)
        # This can show how similar two graphs are and how well the read aligns to the graph

        changedGraph = insertAlignedRead(self.oneSubstitutionGraph, self.readwithOneDifference)

        self.assertFalse(isIsomorphic(self.oneSubstitutionGraph , changedGraph))
        self.assertEqual(calculateGED(self.oneSubstitutionGraph, changedGraph), 1.0)

class MultipleSequenceAlignmentTest(unittest.TestCase):

    def global_align_test(self):

        from gengraph_newfunctions import global_align_by_composition, create_freagment_similarity_table

        seq_file = './test_files/s_pneumoniae.txt'

        output = create_freagment_similarity_table(seq_file, 100, 50)

        self.assertEqual(output, 'dataframe')

    def align_one_read(self):



        self.assertEqual('this', 'this')


class HashTesting(unittest.TestCase):

    def test_one(self):

        self.assertEqual('dataframe', 'dataframe')



if __name__ == '__main__':
    unittest.main()
