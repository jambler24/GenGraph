from gengraph import *

from itertools import combinations, combinations_with_replacement
import multiprocessing as mp

from itertools import product, islice
import datetime
import string

'''

Format for index
[

    [['CATG...CTAG'], [['Aln_79_27', 2094], ['Aln_79_28', 12], ['Aln_79_29', 5]]], 
    [['CATG...CAAT'], [['Aln_79_27', 2094], ['Aln_79_28', 12], ['Aln_79_31', 5]]], 

]

New format for index?
[

    ['CATG...CTAG'], [['Aln_79_27', 2094], ['Aln_79_28', 12], ['Aln_79_29', 5]], 
    ['CATG...CAAT'], [['Aln_79_27', 2094], ['Aln_79_28', 12], ['Aln_79_31', 5]], 

]

'''

# Import of vg graphs


def import_gfa(file_path):

    gfa_file_obj = open(file_path, 'r')

    gg_network = nx.MultiDiGraph()

    for line in gfa_file_obj:

        if line[0] == 'H':
            # Header line
            header_list = line.split('\t')
            format_version = header_list[1].strip()

            print(format_version)

        if line[0] == 'S':
            # Node line
            node_list = line.split('\t')
            node_id = node_list[1]
            node_seq = node_list[2].strip()

            gg_network.add_node(node_id, sequence=node_seq)

    gfa_file_obj = open(file_path, 'r')
    for line in gfa_file_obj:

        if line[0] == 'L':
            # link / edge node
            edge_list = line.split('\t')
            edge_source = edge_list[1]
            edge_target = edge_list[3]
            gg_network.add_edge(edge_source, edge_target)

        if line[0] == 'P':
            # Add path info
            path_info = line.split('\t')
            seq_info = path_info[2]
            path_nodes = path_info[1]

    return gg_network


def get_next_base_older(kmer_length, kmer_matrix, graph_obj):
    """
    Starting at a nucleotide in a node, get all kmers of the required length
    :param kmer_length:
    :param kmer_matrix:
    :param graph_obj:
    :return:
    """

    for a_kmer_list in kmer_matrix:
        #print('---------')
        #print(kmer_matrix)

        if len(a_kmer_list[0][0]) < kmer_length:
            # This k-mer is not at the required length
            #print(a_kmer_list[1])
            last_node = a_kmer_list[1][-1][0]

            last_node_seq = graph_obj.nodes[last_node]['sequence']
            last_node_length = len(last_node_seq)

            node_start_index = a_kmer_list[1][-1][1]

            # Add nucleotides from the current node until either the kmer is long enough or the node sequence ends
            while len(a_kmer_list[0][0]) < kmer_length and node_start_index < last_node_length:

                #print(last_node_seq[node_start_index])
                #print(a_kmer_list[0])
                a_kmer_list[0][0] += last_node_seq[node_start_index]

                node_start_index += 1

            # If it is still not long enough, move onto the next node.
            if len(a_kmer_list[0][0]) < kmer_length:
                # Need to move to next node

                # Find which nodes are next
                neighbors = graph_obj.neighbors(last_node)

                # Need to replace the current one, then add more?
                is_first = True
                for neighbor in neighbors:
                    #print(neighbor)

                    print(a_kmer_list)
                    old_list = pickle.loads(pickle.dumps(a_kmer_list, -1))
                    print('=====')
                    print(old_list)

                    if is_first:
                        # Add new node
                        a_kmer_list[1].append([neighbor, 1])

                        # Add first base of new node
                        a_kmer_list[0][0] += graph_obj.nodes[neighbor]['sequence'][0]

                        is_first = False

                        #print(a_kmer_list)
                        #print(old_list)

                    else:
                        #print('sa')

                        new_list = pickle.loads(pickle.dumps(old_list[1].append([neighbor, 1])))
                        print('---')
                        print(graph_obj.nodes[neighbor])
                        print(neighbor)
                        print(new_list)
                        new_list[0][0] += graph_obj.nodes[neighbor]['sequence'][0]

                        kmer_matrix.append(new_list)

                #print('-=--')
                #print(kmer_matrix)

            #print(kmer_matrix)

    is_complete = True

    # Check if pass here
    #print('pass check')
    for kmer, a_path in kmer_matrix:
        #print(kmer)
        #print(len(kmer[0]))
        #print(a_path)
        if len(kmer[0]) < kmer_length:
            is_complete = False

    #print(is_complete)
    #print(kmer_matrix)

    if is_complete is False:

        get_next_base(kmer_length, kmer_matrix, graph_obj)

    return kmer_matrix


def get_next_base_old(kmer_length, kmer_matrix, graph_obj):
    """
    Starting at a nucleotide in a node, get all kmers of the required length
    :param kmer_length:
    :param kmer_matrix:
    :param graph_obj:
    :return:
    """

    for a_kmer_list in kmer_matrix:
        #print('---------')
        #print(kmer_matrix)

        kmer_seq = a_kmer_list[0][0]
        kmer_last_node = a_kmer_list[1][-1][0]

        if len(kmer_seq) < kmer_length:
            # This k-mer is not at the required length
            #print(a_kmer_list[1])

            last_node_seq = graph_obj.nodes[kmer_last_node]['sequence']
            last_node_length = len(last_node_seq)

            node_start_index = a_kmer_list[1][-1][1]

            # Add nucleotides from the current node until either the kmer is long enough or the node sequence ends
            while len(a_kmer_list[0][0]) < kmer_length and node_start_index < last_node_length:

                #print(last_node_seq[node_start_index])
                #print(a_kmer_list[0])
                a_kmer_list[0][0] += last_node_seq[node_start_index]

                node_start_index += 1

            # If it is still not long enough, move onto the next node.
            if len(a_kmer_list[0][0]) < kmer_length:
                # Need to move to next node

                # Find which nodes are next
                neighbors = graph_obj.neighbors(kmer_last_node)


                # Need to replace the current one, then add more?
                is_first = True
                for neighbor in neighbors:
                    print(neighbor)

                    print(a_kmer_list)
                    old_list = pickle.loads(pickle.dumps(a_kmer_list, -1))
                    print('=====')
                    print(old_list)

                    if is_first:
                        # Add new node
                        a_kmer_list[1].append([neighbor, 1])

                        # Add first base of new node
                        a_kmer_list[0][0] += graph_obj.nodes[neighbor]['sequence'][0]

                        is_first = False

                        #print(a_kmer_list)
                        #print(old_list)

                    else:
                        #print('sa')

                        new_list = pickle.loads(pickle.dumps(old_list[1].append([neighbor, 1])))
                        print('---')
                        print(graph_obj.nodes[neighbor])
                        print(neighbor)
                        print(new_list)
                        new_list[0][0] += graph_obj.nodes[neighbor]['sequence'][0]

                        kmer_matrix.append(new_list)

                #print('-=--')
                #print(kmer_matrix)

            #print(kmer_matrix)

    is_complete = True

    # Check if pass here
    #print('pass check')
    for kmer, a_path in kmer_matrix:
        #print(kmer)
        #print(len(kmer[0]))
        #print(a_path)
        if len(kmer[0]) < kmer_length:
            is_complete = False

    #print(is_complete)
    #print(kmer_matrix)

    if is_complete is False:

        get_next_base(kmer_length, kmer_matrix, graph_obj)

    return kmer_matrix


def get_next_base(kmer_length, kmer_matrix, graph_obj):
    """
    Starting at a nucleotide in a node, get all kmers of the required length
    :param kmer_length:
    :param kmer_matrix:
    :param graph_obj:
    :return:
    """

    new_matrix = []

    for a_kmer_list in kmer_matrix:

        kmer_seq = a_kmer_list[0][0]
        kmer_last_node = a_kmer_list[1][-1][0]

        if len(kmer_seq) < kmer_length:
            # This k-mer is not at the required length

            last_node_seq = graph_obj.nodes[kmer_last_node]['sequence']
            last_node_length = len(last_node_seq)

            node_start_index = a_kmer_list[1][-1][1]

            # Add nucleotides from the current node until either the kmer is long enough or the node sequence ends
            while len(a_kmer_list[0][0]) < kmer_length and node_start_index < last_node_length:

                a_kmer_list[0][0] += last_node_seq[node_start_index]

                node_start_index += 1

            # If it is still not long enough, move onto the next node.
            if len(a_kmer_list[0][0]) < kmer_length:

                # Need to move to next node

                # Find which nodes are next
                neighbors = graph_obj.neighbors(kmer_last_node)

                # Need to replace the current one, then add more?
                for neighbor in neighbors:

                    new_kmer_list = pickle.loads(pickle.dumps(a_kmer_list, -1))

                    new_kmer_list[1].append([neighbor, 1])

                    new_kmer_list[0][0] += graph_obj.nodes[neighbor]['sequence'][0]

                    new_matrix.append(new_kmer_list)
            else:

                new_matrix.append(a_kmer_list)

        else:

            new_matrix.append(a_kmer_list)

    is_complete = True

    # Check if pass here
    for kmer, a_path in new_matrix:

        if len(kmer[0]) < kmer_length:
            is_complete = False

    if is_complete is False:

        return get_next_base(kmer_length, new_matrix, graph_obj)

    else:
        return new_matrix


def get_node_kmers(a_node, graph_obj, kmer_length, return_structure):
    """
    Returns a list structure of all possible kmers from a node, including the positions, and nodes that the sequences
    extend into.
    :param a_node: A string node ID from the graph object
    :param graph_obj: A GenGraph graph object
    :param kmer_length:
    :param return_structure:
    :return:
    """

    # Reversed nodes?
    node_length = len(graph_obj.nodes[a_node]['sequence'])
    current_base_pos = 1

    kmer_matrix = []

    while current_base_pos <= node_length:

        kmer_matrix.append([
                [graph_obj.nodes[a_node]['sequence'][current_base_pos - 1]],
                [[a_node, current_base_pos]]
            ])

        kmer_matrix = get_next_base(kmer_length, kmer_matrix, graph_obj)

        current_base_pos += 1

    kmer_matrix = get_next_base(kmer_length, kmer_matrix, graph_obj)

    if return_structure == 'list':

        return kmer_matrix

    elif return_structure == 'kmer_dict':

        kmer_dict = {}

        for a_kmer in kmer_matrix:
            if a_kmer[0][0] in kmer_dict.keys():
                kmer_dict[a_kmer[0][0]] += [a_kmer[1]]
            else:
                kmer_dict[a_kmer[0][0]] = [a_kmer[1]]

        return kmer_dict







'''
[
    [
        ['ACTTCGACGACTTCGACGAT'], 
        [
            ['Aln_79_26', 1], ['Aln_79_27', 1]
        ]
    ], 
    [['CTTCGACGACTTCGACGATA'], [['Aln_79_26', 2], ['Aln_79_27', 1]]], 
    [['TTCGACGACTTCGACGATAA'], [['Aln_79_26', 3], ['Aln_79_27', 1]]], 
    [['TCGACGACTTCGACGATAAG'], [['Aln_79_26', 4], ['Aln_79_27', 1]]], 
    [['CGACGACTTCGACGATAAGG'], [['Aln_79_26', 5], ['Aln_79_27', 1]]], 
    [['GACGACTTCGACGATAAGGG'], [['Aln_79_26', 6], ['Aln_79_27', 1]]], 
    [['ACGACTTCGACGATAAGGGC'], [['Aln_79_26', 7], ['Aln_79_27', 1]]], 
    [['CGACTTCGACGATAAGGGCC'], [['Aln_79_26', 8], ['Aln_79_27', 1]]], 
    [['GACTTCGACGATAAGGGCCG'], [['Aln_79_26', 9], ['Aln_79_27', 1]]]

]


'''


def create_query_kmers_OLD_FAST(q_sequence, kmer_size):

    q_kmers = [q_sequence[x:y] for x, y in combinations(range(len(q_sequence) + 1), r=2)
               if len(q_sequence[x:y]) == kmer_size]

    return q_kmers


def create_query_kmers(q_sequence, kmer_size):

    # TODO: Deal with lost tailing sequence. Lasy kmer size chunk not used.

    q_kmers = {}
    window_start = 0
    window_stop = kmer_size
    q_seq_length = len(q_sequence['sequence'])

    while window_stop <= q_seq_length:
        q_kmers[q_sequence['sequence'][window_start:window_stop]] = window_start
        window_start += 1
        window_stop += 1

    # Try clean up
    #del window_stop
    #del window_start
    #del q_seq_length

    return q_kmers


def create_kmer_dict(in_graph_obj, kmer_size):
    """
    Takes in the genome graph and extracts all possible k-mers including their positions
    :param in_graph_obj: Graph created by GenGraph
    :param kmer_size: The size of the k-mers to be created
    :return: a dict of kmers, containing lists that describe where the k-mer is found in the reference graph.

    Structure of returned dict:

    {
        k-mer: [
            [
                ['nodes', 'align start position'],
            ]
        ]
    }

    EG:
    {'TAAACAACGGCCCCGACCCC':
        [
            [['Aln_69_59', 1676], ['Aln_69_60', 1], ['H37Rv_104', 1], ['Aln_70_1', 1]],
            [['Aln_69_59', 1676], ['Aln_230', 1]]
        ]
    }


    """

    total_nodes = len(in_graph_obj.nodes())

    count = 0

    all_kmer_positions = {}

    pool = mp.Pool(mp.cpu_count())

    per_node_kmer_list = pool.starmap(get_node_kmers, [(a_node, in_graph_obj, kmer_size, 'kmer_dict') for a_node in in_graph_obj.nodes()])

    for a_node_kmers in per_node_kmer_list:
        for key, value in a_node_kmers.items():
            if key in all_kmer_positions.keys():
                all_kmer_positions[key] += value
            else:
                all_kmer_positions[key] = value

    return all_kmer_positions


def create_kmer_graph(in_graph_obj, kmer_size):

    # TODO: Multi processing here

    total_nodes = len(in_graph_obj.nodes())

    count = 0

    all_kmer_positions = {}

    pool = mp.Pool(mp.cpu_count())

    per_node_kmer_list = pool.starmap(get_node_kmers, [(a_node, in_graph_obj, kmer_size, 'kmer_dict') for a_node in in_graph_obj.nodes()])

    for a_node_kmers in per_node_kmer_list:
        for key, value in a_node_kmers.items():
            if key in all_kmer_positions.keys():
                all_kmer_positions[key] += value
            else:
                all_kmer_positions[key] = value

    return all_kmer_positions


def align_seq_hash(q_sequence, ref_hash_dict, kmer_size, use_qual=True):
    """
    Breaks up the query sequence into k-mers of given size matching that of the reference hash dictionary. Finds the
    kmer in the ref_hash_dict, returns the k-mers and any matching positions.
    :param q_sequence: Query sequence as produced by the process_fastq_lines() function representing the 4 lines of a
    fastq file in a list.
    :param ref_hash_dict: The reference genome hash created by create_kmer_dict()
    :param kmer_size: The size of the k-mers to be used
    :param use_qual:
    :return: Dictionary of mapped positions

    {kmer: [
        'k-mer_first_nucleotode_seq',
        'quality_string',
            [
                [aligned_node, aligned_node_position],
                [aligned_node2, aligned_node_position2]
            ]
        ]
    }

    """

    exact_align = 0
    multi_align = 0
    no_align = 0
    reversed_align = 0

    q_kmers = create_query_kmers(q_sequence, kmer_size)

    q_kmer_dict = {k: [v] for v, k in enumerate(q_kmers)}

    if use_qual:
        for a_seq_kmer, val_list in q_kmer_dict.items():
            q_kmer_dict[a_seq_kmer].append(q_sequence['quality'][val_list[0]])

    for kmer, position in q_kmer_dict.items():

        try:
            kmer_graph_pos = ref_hash_dict[kmer]

        except KeyError:
            kmer_graph_pos = []

        q_kmer_dict[kmer] += kmer_graph_pos

    # Try clean up a bit
    #del q_kmers

    return q_kmer_dict


def calcGCcontent():

    return GCprecentage


def testHash(nuc_string):

    import hashlib

    hashVal = abs(hash(nuc_string[0])) % (10 ** 5)

    return hashVal


def create_hash_info(in_matrix, method='testHash'):

    out_hash_dict = {}

    for k_mer in in_matrix:

        if method == 'testHash':

            hash_val = testHash(k_mer[0])

            out_hash_dict[hash_val] = k_mer

    return out_hash_dict


def create_encoding_dict(kmer_length, mode='default'):
    encoding_dicts = {'encode': {}, 'decode': {}}
    characters = 'ACGT'
    end_char = '-'
    list_of_combinations = []

    for a_kmer in combinations_with_replacement(characters, kmer_length):
        kmer_string = ''.join(a_kmer)
        count = 1
        while count < kmer_length:
            #kmer_string_end = ''.join(a_kmer[count:])
            #kmer_string_end = count * '-' + kmer_string_end
            #list_of_combinations.append(kmer_string_end)

            kmer_string_end = ''.join(a_kmer[:-1 * count])
            kmer_string_end = kmer_string_end + count * end_char
            if kmer_string_end not in list_of_combinations:
                list_of_combinations.append(kmer_string_end)

            count += 1
        list_of_combinations.append(kmer_string)

    kmer_count = 0

    while kmer_count < len(list_of_combinations):

        encoding_dicts['encode'][list_of_combinations[kmer_count]] = string.printable[kmer_count]
        encoding_dicts['decode'][string.printable[kmer_count]] = list_of_combinations[kmer_count]

        kmer_count += 1

    return encoding_dicts


def encode_nucleotides(nuc_string, encoding_dict, kmer_size):

    encoded_result = ""

    for a_kmer in islice(nuc_string, kmer_size, None, kmer_size):
        try:
            print(encoding_dict[a_kmer])
            encoded_result += encoding_dict[a_kmer]
        except KeyError:
            while len(a_kmer) < kmer_size:
                a_kmer += '-'
            encoded_result += encoding_dict[a_kmer]

    return encoded_result


'''
a_matrix = get_node_kmers('Aln_79_27', graph_obj, 20)
print(a_matrix)
print(len(a_matrix))
quit()
'''
# ----------------------------------------------------------------- ><

''' 
# ------------------------- Testing the use of encoding to compress the created graph -------------------------

encode_dict = create_encoding_dict(4)

test_string = 'GCAGATCGAGCCTACGGCTACGGACGCGGCGGCGGCATATACGCATACGACTACTCTATACTCGG'

encoded_test_string = encode_nucleotides(test_string, encode_dict['encode'], 4)

print(sys.getsizeof(test_string))
print(sys.getsizeof(encoded_test_string))

print(test_string)
print(encoded_test_string)
print(encode_dict['decode']['p'])
print(encode_dict)
quit()
'''

path_to_GG_file = 'test_files/latest2genome.xml'

#path_to_reads = 'test_files/minSRR1144793.fastq'
path_to_reads = '/Users/panix/data/gg_data/MTB__S507_LFO46Pool91_3128__L8_GTCCGC_L008_R2_001.fastq'
#path_to_reads = '/Users/panix/data/gg_data/160k_subset.fastq'

graph_obj = import_gg_graph(path_to_GG_file)

#begin_time = datetime.datetime.now()

#kmer_dict_out = create_kmer_dict(graph_obj, 20)

#print(datetime.datetime.now() - begin_time)

#print(kmer_dict_out)

#with open('kmer_dict_multi.pickle', 'wb') as handle:
#    pickle.dump(kmer_dict_out, handle, protocol=pickle.HIGHEST_PROTOCOL)

kmer_dict_out = pickle.load(open("kmer_dict_multi.pickle", "rb"))

begin_time = datetime.datetime.now()


def process_fastq_lines(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def align_fastq_to_kmer_graph(fastq_file, reference_kmer_dict):
    """
    This function seeks to align reads from a fastq file into a k-mer graph.
    :param fastq_file: Path to the fastq file
    :param reference_kmer_dict: A k-mer dict created before hand, used to kmow the position of k-mers in the reference
    genome. For example using the create_kmer_dict function.
    :return:
    """

    set_kmer = 20

    graph_write_limit = 16000

    graph_write_count = 0

    graph_read_count = 0

    in_fasta = open(fastq_file, 'r')

    # networkx, none, neo4j
    graph_db = 'networkx'

    slim_graph = True

    if graph_db == 'networkx':

        out_graph = nx.MultiGraph()

    n = 4
    lines = []



    for a_line in in_fasta:

        lines.append(a_line.rstrip())
        if len(lines) == n:

            graph_read_count += 1

            fastq_record = process_fastq_lines(lines)
            lines = []

            # Align the sequence to hash

            align_res = align_seq_hash(fastq_record, reference_kmer_dict, set_kmer)
            # Returns the aligmnent positions of all kmers of the fastq_record

            # process alignment result
            previous_node = False
            for key, val in align_res.items():

                try:
                    # Check if an alignment to the ref was found.
                    kmer_node_name = val[2][0][0] + '-' + str(val[2][0][1])

                    if kmer_node_name not in out_graph.nodes():
                        ref_node_pos = val[2][0][0] + '-' + str(val[2][0][1])
                        if slim_graph:
                            out_graph.add_node(kmer_node_name)
                        else:
                            node_dict = {'nuc': key[0], 'ref': ref_node_pos, 'kmer': key, 'qual': val[1]}
                            out_graph.add_node(kmer_node_name, **node_dict)
                    else:
                        if not slim_graph:
                            old_qual = out_graph.nodes[kmer_node_name]['qual']
                            out_graph.nodes[kmer_node_name]['qual'] = old_qual + ',' + val[1]

                except IndexError:
                    # Create a new node for the query seq if no alignment to the ref was found.
                    kmer_node_name = key

                    if kmer_node_name not in out_graph.nodes():
                        if slim_graph:
                            out_graph.add_node(kmer_node_name)
                        else:
                            ref_node_pos = 'alt' + str(val[0])
                            node_dict = {'nuc': key[0], 'ref': ref_node_pos, 'kmer': key, 'qual': val[1]}
                            out_graph.add_node(kmer_node_name, **node_dict)

                    else:
                        if not slim_graph:
                            old_qual = out_graph.nodes[kmer_node_name]['qual']
                            out_graph.nodes[kmer_node_name]['qual'] = old_qual + ',' + val[1]

                if previous_node is not False:
                    # This links k-mers from the same read
                    if out_graph.has_edge(previous_node, kmer_node_name):
                        current_weight = out_graph[previous_node][kmer_node_name]['a']['weight']
                        out_graph[previous_node][kmer_node_name]['a']['weight'] = current_weight + 1
                    else:
                        out_graph.add_edge(previous_node, kmer_node_name, key='a', weight=1)

                previous_node = kmer_node_name

            if graph_read_count == graph_write_limit:

                print('read align limit')

                graph_write_count += 1

                #with open('temp/' + str(graph_write_count) + '_aligned_kmer.pickle', 'wb') as handle:
                #    pickle.dump(out_graph, handle, protocol=pickle.HIGHEST_PROTOCOL)

                del out_graph

                print(graph_write_count)

                out_graph = nx.MultiGraph()

                graph_read_count = 0


    return out_graph


res_of_the_thing = align_fastq_to_kmer_graph(path_to_reads, kmer_dict_out)

print(datetime.datetime.now() - begin_time)

#nx.write_graphml(res_of_the_thing, 'aligned_kmer.xml')

# --------------------- Working up to here --------------------- --------------------- ---------------------

'''
# Strange output 
for kmer, list in kmer_dict_out.items():
    if len(list) > 1:
        print(kmer_dict_out[kmer], kmer)
'''

quit()

node_dbg = nx.Graph()

total = len(graph_obj.nodes())

count = 0

for a_node in graph_obj.nodes():

    print(a_node)
    print(str(count), '/', str(total))
    count += 1

    if count == 100:
        quit()

    a_matrix = get_node_kmers(a_node, graph_obj, 20)

    #index_file = open('graphIndex.txt', 'w')

    #hash_matrix = create_hash_info(a_matrix)

    last_node = ''

    for kmer in a_matrix:

        #index_file.write(str(kmer) + '\n')

        node_name = ''

        for kmer_node in kmer[1]:

            #print(kmer_node)
            kmer_node_string = kmer_node[0] + '_' + str(kmer_node[1]) + ":"
            node_name += kmer_node_string

        node_name = node_name[:-1]
        node_dbg.add_node(node_name, sequence=kmer[0][0])

        if len(last_node) > 0:
            node_dbg.add_edge(last_node, node_name)

        last_node = node_name

        #print(kmer[0][0])

for a_node in graph_obj.nodes():

    node_start_name = str(a_node) + '_' + '1'

    node_len = len(graph_obj.nodes[a_node]['sequence'])

    node_end_name




nx.write_graphml(node_dbg, 'dbg_test.xml')

quit()

#path_to_w148 = '/Users/panix/Desktop/genomes/w148/W_148_NCBI.fa'
#path_to_h37rv = '/Users/panix/Desktop/genomes/H37Rv/sequence.fasta'

#w148_genome = input_parser(path_to_w148)
#h37rv_genome = input_parser(path_to_h37rv)

for a_node in graph_obj.nodes():
    print(a_node)

graph_obj.plot_subgraph(1557101, 1565448, 'H37Rv', neighbours=1)

quit()

'''
G = graph_obj.get_region_subgraph(3232, 13033, 'W_148', neighbours=1)

pos=nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_size = 500)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, arrows=True)
plt.show()

'''


#a_seq = graph_obj.get_sequence(1289, 1289, 'W_148')

# pos node, all in
# a_seq = graph_obj.get_sequence(120, 130, 'H37Rv')

# neg node, all in


window_size = 100
step_size = 100
end_bp = 200000

#start_pos = 1191351
#stop_pos = 1191450
iso = 'H37Rv'

a_seq = graph_obj.get_sequence(1, 11, iso)
print(a_seq)
print('-----')

count = 1
while count < end_bp:
    print('-------')
    print(count)
    window_end = count + window_size
    print(window_end)

    graph_ex_seq = graph_obj.get_sequence(count, window_end, iso)

    string_ex_seq = h37rv_genome[0]['DNA_seq'][count-1: window_end]

    count += step_size

    if graph_ex_seq != string_ex_seq:
        print('Fail')
        print(graph_ex_seq)
        print(string_ex_seq)

        quit()


a_seq = graph_obj.get_sequence(start_pos, stop_pos, iso)


print(a_seq)
print(len(a_seq))
print('')
print(h37rv_genome[0]['DNA_seq'][start_pos-1: stop_pos])


