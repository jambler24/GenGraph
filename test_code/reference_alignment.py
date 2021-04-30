# This is a test file


# Importing needed things

from GenGraph import *
from module_testing import *




# Import test graph

test_graph_path = "latest2genome.xml"

in_graph = import_gg_graph(test_graph_path)

# Create k-mer dict from the test graph
# The create_kmer_dict function is from module_testing

kmer_dict_out = create_kmer_dict(graph_obj, 20)


print('ok')
