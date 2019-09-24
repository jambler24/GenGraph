from gengraph import *
from gengraphTool import *
from readAligner import *
import networkx as nx
import matplotlib.pyplot as plt

graph_obj = import_gg_graph('./test_kmer.xml')
kzn_obj = import_gg_graph('./kzn_2_pangenome.xml')
print(graph_obj.ids())


print(graph_obj.number_of_nodes())

print(graph_obj.nodes)


alignerGraph = Aligner(graph_obj)



#GED = nx.graph_edit_distance(graph_obj,onesub)

#print(graph_obj.get_sequence(0, 1000, 'H37Rv'))


#kznNodes = list(kznsub.nodes)
#node1 = kznNodes[0]
#print(kznNodes)

#testGraph = nx.Graph()

#testGraph.add_node(0, sequence = ["Aln_1","TTG"])
#testGraph.add_node(1,sequence = ["Aln_2","TTTTAAAAAG"])
#testGraph.add_edge(0,1)
#sequence = ["Aln_2","TTAAGACT"]
#nx.set_node_attributes(testGraph,sequence,'sequence')

#print(testGraph.nodes[0]['sequence'])

#print(graph_obj.ids())
#subgraph = graph_obj.get_region_subgraph(0, 10, 'H37Rv')


nx.draw(graph_obj, with_labels=True)
plt.show()

#teststring = 'TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAGCAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCCCCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGATCCAACTCGGGGTCCGCATCGCTCCGCCGGCGACCGACGAAGCCGACGACACTACCGTGCCGCCTTCCGAAAATCCTGCTACCACATCGCCAGACACCACAACCGACAACGACGAGATTGATGACAGCGCTGCGGCACGGGGCGATAACCAGCACAGTTGGCCAAGTTACTTCACCGAGCGCCCGCACAATACCGATTCCGCTACCGCTGGCGTAACCAGCCTTAACCGTCGCTACACCTTTGATACGTTCGTTATCGGCGCCTCCAACCGGTTCGCGCACGCCGCCGCCTTGGCGATCGCAGAAGCACCCGCCCGCGCTTACAACCCCCTGTTCATCTGGGGCGAGTCCGGTCTCGGCAAGACACACCTGCTACACGCGGCAGGCAACTATGCCCAACGGTTGTTCCCGGGAATGCGGGTCAAATATGTCTCCACCGAGGAATTCACCAACGACTTCATTAACTCGCTCCGCGATGACCGCAAGGTCGCATTCAAACGCAGCTACCGCGACGTAGACGTGCTGTTGGTCGACGACATCCAATTCATTGAAGGCAAAGAGGGTATTCAAGAGGAGTTCTTCCACACCTTCAACACCTTGCACAATGCCAACAAGCAAATCGTCATCTCATCTGACCGCCCACCCAAGCAGCTCGCCACCCTCGAGGACCGGCTGAGAACCCGCTTTGAGTGGGGGCTGATCACTGACGTACAACCACCCG'
#first100 = teststring[0:100]
#onehundredtofourhundred = teststring[100:400]
#reverse = onehundredtofourhundred[::-1]
#fourhundredtoend = teststring[400:]
#final = first100 + reverse + fourhundredtoend
#print(final)