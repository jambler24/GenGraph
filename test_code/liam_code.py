
# -----------------------------------------------------"SequenceHomology" functions added below

def nucleotide_sequence_alignment(pos1, pos2, path, isolate1, isolate2):

    """extracts the same gene sequences from two isolates and performs a pairwise alignment
                :param pos1: start position of gene in isolate1
                :param pos2: end position of gene in isolate1
                :param path: path to genome graph contained within XML file containing isolate1 and isolate2
                :param isolate1: the first strain from the genome graph being compared (i.e. H37Rv in this case since it is being used as a reference)
                :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
                :return:

        subseq1 - nucleotide sequence string of gene found in isolate1
        subseq2_coords - the converted coordinates to allow extraction of the same gene sequence in isolate2
        subseq2 - nucleotide sequence string of gene found in isolate2
        score - the similarity score of the two gene sequences
        nucleotide_alignment - the pairwise alignment of the two nucleotide alignments
        nucleotide_matrix - a NumPy array version of the pairwise alignment to allow iteration for mutation detection in later functions
        isolate1_gene - the alignment array version of "subseq1"
        isolate2_gene - the alignment array version of "subseq2"

        NOTES
        -------------

        this function is only able to align two similar sequences and doesn't allow for multiple sequence alignment needed for
        sequence similarity assessment. Future work should be directed toward adding this capability to the function.

        EXAMPLE
        -------------

        comparing the carB gene sequence from the H37Rv and H37Ra strains extracted from an XML file named "H37R_pangenome.xml"

        nucleotide_sequence_alignment(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')


    """

    # IMPORT .XML FILE AND EXTRACT GENE SEQUENCE FROM THE DIFFERENT STRAINS

    graph_obj = import_gg_graph(path)
    # extract subgraphs from both isolates you want to compare for a specified nucleotide range
    # isolate1 needs to be the ancestral strain and isolate2 the derived strain in order for this function to work
    subseq1 = extract_original_seq_region_fast(graph_obj, pos1, pos2, isolate1)
    # "pos1" and "pos2" are relative to "isolate1" and so need to convert coords to get gene sequence in "isolate2"
    # so that similar regions are being compared between the strains
    subseq2_coords = convert_coordinates(graph_obj, pos1, pos2, isolate1, isolate2)
    subseq2_coords = list(subseq2_coords.values())
    subseq2 = extract_original_seq_region_fast(graph_obj, subseq2_coords[0], subseq2_coords[1], isolate2)

    # mimic deletion in carB gene for "deletion_detection" function example

    # subseq2 = subseq2.replace('CCC', '', 1)

    # mimic insertion in carB gene for "insertion_detection" function example

    # subseq1 = subseq1.replace('CCC', '', 1)

    # mimic SNP in carB gene for "substitution_detection" function example

    # subseq1 = subseq1.replace('T', 'C', 1)

    # PERFORM NUCLEOTIDE ALIGNMENT using Biopython module

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -0.5
    nucleotide_alignment = aligner.align(subseq1, subseq2)
        [0  ]# produces many different alignments for the same two sequences, of which the first one will be chosen
    list_alignment = list(str
        (nucleotide_alignment).splitlines()  )# convert alignment to string to be able to be to loop through each nucleotide - 'PairwiseAlign is not iterable'
    # create an array where each character in alignment gets its own index
    isolate1_gene = np.array(list(list_alignment[0]))
    isolate2_gene = np.array(list(list_alignment[2]))
    aligned = np.array(list(list_alignment[1]))
    nucleotide_matrix = np.row_stack((isolate1_gene, aligned, isolate2_gene))
    matches = sum(np.char.count(aligned, '|'))
    score = "Similarity = %.1f:" % (matches / (len(subseq1)) * 100)

    return subseq1, subseq2_coords, subseq2, score, nucleotide_alignment, nucleotide_matrix, isolate1_gene, isolate2_gene


def protein_sequence_alignment(pos1, pos2, path, isolate1, isolate2):

    """converts extracted gene sequences in "nucleotide_alignment" to protein sequences and performs protein pairwise alignment
    :param pos1: start position of gene in isolate1
    :param pos2: end position of gene in isolate1
    :param path: path to genome graph contained within XML file containing isolate1 and isolate2
    :param isolate1: the first strain from the genome graph being compared(i.e. H37Rv in this case since it is being used as a reference)
    :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
    :return:
    protein_alignment: the pairwise alignment of the amino acid sequences converted from the nucleotide sequences
    of "subseq1" ("isolate1" protein) and "subseq2" ("isolate2" protein)
    protein_matrix: a NumPy array version of the pairwise alignment to allow iteration for amino acid change detection
    isolate1_protein: a Numpy array version of the "isolate1" protein - the first row of "protein_matrix"
    aligned_protein: a NumPy array version of the alignment of the two proteins - the second row of "protein_matrix"
    isolate2_protein: a NumPy array version of the "isolate2" protein - the third row of "protein_matrix"
    protein1_amino_acids: a list version of "isolate1_protein"
    protein2_amino_acids: a list version of "isolate2_protein"

    NOTES
    -----------------

    this function is only able to align two similar protein sequences and doesn't allow for multiple sequence alignment needed for
    sequence similarity assessment. Future work should be directed toward adding this capability to the function.

    An original aligner was developed to allow for protein pairwise alignment but may presents bug as unit testing has
    not been done

    EXAMPLE
    -----------------

    comparing the translated versions of the carB gene sequences from the H37Rv and H37Ra strains extracted from an XML file named "H37R_pangenome.xml"

    protein_sequence_alignment(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')

    """

    subseq1, subseq2_coords, subseq2, score, nucleotide_alignment, nucleotide_matrix, isolate1_gene, isolate2_gene = nucleotide_sequence_alignment \
        (pos1, pos2, path, isolate1, isolate2)
    print(isolate1_gene)
    print(isolate2_gene)
    # ALIGNMENT - Protein sequence

    # create protein sequences from nucleotide sequences "subseq1" and "subseq2"
    # "_" in table is a stop codon

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    protein1 = ""
    protein2 = ""
    subseq1_codons = [subseq1[i:i + 3] for i in range(0, len(subseq1), 3)]
    subseq2_codons = [subseq2[i:i + 3] for i in range(0, len(subseq2), 3)]

    for i in range(len(subseq1_codons)):
        if len(subseq1_codons[i]) % 3 == 0:
            protein1 += table[subseq1_codons[i]]
        else:
            protein1 += '-'

    for i in range(len(subseq2_codons)):
        if len(subseq2) % 3 == 0:
            protein2 += table[subseq2_codons[i]]
        else:
            protein2 += '-'

    protein1_amino_acids = list(protein1)
    protein2_amino_acids = list(protein2)

    # tried to perform protein pairwise alignment using Biopython but presented a lot of bugs - created an alignment tool below
    # the commented-out code

    # perform protein pairwise alignment - #since we are looking at coding regions, let's take the two protein sequences and align them to see which amino acids are different between them
    # aligner = Align.PairwiseAligner()
    # aligner.mode = 'local'
    # aligner.open_gap_score = -0.5
    # aligner.extend_gap_score = -1
    # protein_alignment = aligner.align(protein1, protein2)[0]
    # list_protein_alignment = list(str(protein_alignment).splitlines())  # convert alignment to string to be able to be to loop through each nucleotide - 'PairwiseAlign is not iterable'
    # create an array where each character in alignment gets its own index
    # ancestral_protein = np.array(list(list_protein_alignment[0]))
    # derived_protein = np.array(list(list_protein_alignment[2]))
    # aligned_protein = np.array(list(list_protein_alignment[1]))
    # protein_matrix = np.row_stack((ancestral_protein, aligned_protein, derived_protein))

    # here is the protein aligner that was developed to allow more accurate comparison of the protein sequences

    stringOne = protein1
    stringTwo = protein2

    finalStringOne = ''
    middleString = ''
    finalStringTwo = ''

    pointerOne = 0
    pointerTwo = 0
    for i in range(len(stringOne)):
        if pointerOne != len(stringOne) - 1 and pointerTwo != len(stringOne) - 1:
            if stringOne[pointerOne] == stringTwo[pointerTwo]:
                finalStringOne = finalStringOne + str(stringOne[pointerOne])
                finalStringTwo = finalStringTwo + str(stringTwo[pointerTwo])
                middleString = middleString + '|'
                pointerOne += 1
                pointerTwo += 1

        elif stringOne[pointerOne] != stringTwo[pointerTwo]:
            if pointerOne == 0 and pointerTwo == 0 and stringOne[pointerOne] != stringTwo[pointerTwo] and stringOne
                [pointerOne + 1] == stringTwo[pointerTwo + 1]:
                finalStringOne = finalStringOne + str(stringOne[pointerOne])
                finalStringTwo = finalStringTwo + str(stringTwo[pointerTwo])
                middleString = middleString + 'X'
                pointerOne += 1
                pointerTwo += 1

        elif pointerOne > 0 and pointerTwo > 0 and stringOne[pointerOne - 1] == stringTwo[pointerTwo - 1]:
            if stringOne[pointerOne + 1] == stringTwo[pointerTwo + 1] and stringOne[pointerOne + 1] != stringTwo
                [pointerTwo] and stringOne[pointerOne] != stringTwo[pointerTwo + 1]:
                finalStringOne = finalStringOne + str(stringOne[pointerOne])
                finalStringTwo = finalStringTwo + str(stringTwo[pointerTwo])
                middleString = middleString + 'X'
                pointerOne += 1
                pointerTwo += 1
            elif stringOne[pointerOne] == stringTwo[pointerTwo + 1] and stringOne[pointerOne + 1] == stringTwo
                [pointerTwo + 2]:
                check = True
            while check:
                finalStringOne = finalStringOne + '-'
                finalStringTwo = finalStringTwo + str(stringTwo[pointerTwo])
                middleString = middleString + '-'
                pointerTwo += 1
                if stringOne[pointerOne] == stringTwo[pointerTwo]:
                    check = False
                    break
        elif stringOne[pointerOne] == stringTwo[pointerTwo + 2] and stringOne[pointerOne + 1] == stringTwo
            [pointerTwo + 3] and stringOne[pointerOne + 2] == stringTwo[pointerTwo + 4]:
            finalStringOne = finalStringOne + '--'
            finalStringTwo = finalStringTwo + str(stringTwo[pointerTwo]) + str(stringTwo[pointerTwo + 1])
            middleString = middleString + '--'
            pointerTwo += 2

        elif stringOne[pointerOne + 1] == stringTwo[pointerTwo] and stringOne[pointerOne + 2] == stringTwo
            [pointerTwo + 1]:
            check = True
        while check:
            finalStringOne = finalStringOne + str(stringOne[pointerOne])
            finalStringTwo = finalStringTwo + '-'
            middleString = middleString + '-'
            pointerOne += 1
            if stringOne[pointerOne + 1] == stringTwo[pointerTwo + 1]:
                check = False
            break

        elif stringOne[pointerOne + 2] == stringTwo[pointerTwo] and stringOne[pointerOne + 3] == stringTwo
            [pointerTwo + 1] and stringOne[pointerOne + 4] == stringTwo[pointerTwo + 2]:
        finalStringOne = finalStringOne + str(stringOne[pointerOne]) + str(stringOne[pointerOne + 1])
        finalStringTwo = finalStringTwo + '--'
        middleString = middleString + '--'
        pointerOne += 2

else:
# and pointerOne == range(len(stringOne)) and pointerTwo == range(len(stringOne)):

if stringOne[-1] != stringTwo[-1]:
    if stringOne[pointerOne] != '_' or stringTwo[pointerTwo] != '_':
        finalStringOne = finalStringOne + str(stringOne[-1])
        finalStringTwo = finalStringTwo + str(stringTwo[-1])
        middleString = middleString + 'X'
else:
    if stringOne[pointerOne] == '_':
        finalStringOne = finalStringOne + '-'
        finalStringTwo = finalStringTwo + str(stringTwo[-1])
        middleString = middleString + '-'
    elif stringTwo[pointerTwo] == '_':
        finalStringOne = finalStringOne + str(stringOne[-1])
        finalStringTwo = finalStringTwo + '-'
        middleString = middleString + '-'

elif stringOne[-1] == stringTwo[-1]:
finalStringOne = finalStringOne + str(stringOne[-1])
finalStringTwo = finalStringTwo + str(stringTwo[-1])
middleString = m

protein_alignment = finalStringOne + '\n' + middleString + '\n' + finalStringTwo

isolate1_protein = np.array(list(finalStringOne))
aligned_protein = np.array(list(middleString))
isolate2_protein = np.array(list(finalStringTwo))
protein_matrix = np.row_stack((isolate1_protein, aligned_protein, isolate2_protein))

return protein_alignment, protein_matrix, isolate1_protein, aligned_protein, isolate2_protein, protein1_amino_acids, protein2_amino_acids


def substitution_detection(pos1, pos2, path, isolate1, isolate2):

    """uses the nucleotide and protein pairwise alignments to detect substitutions and any codon/amino acid changes that occur
    :param pos1: start position of gene in isolate1
    :param pos2: end position of gene in isolate1
    :param path: path to genome graph contained within XML file containing isolate1 and isolate2
    :param isolate1: the first strain from the genome graph being compared(i.e. H37Rv in this case since it is being used as a reference)
    :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
    :return:
    df_substitutions: Pandas DataFrame reporting any substitutions that occurred between the "isolate1" and "isolate2"
    versions of the gene

    NOTES
    -------------

    This function allows for detection of substitutions and the effects these substitutions have on coding sequences.
    However, in order to assess the biological significance of these substitutions, the ability for
    multiple sequence alignment to be performed in the function is needed.

    This function identifies any amino acid changes that occur due to substitutions.
    This function can't discern whether an amino acid change is conserved or not
    (i.e. a SNP that leads to an amino acid change from leucine to isoleucine represents a conservative mutation and is seen as a mismatch)

    positions where substitutions occur are reported relative to "isolate1" which will most likely be the H37Rv since it
    is used as a reference

    EXAMPLE
    -------------

    using the carB gene sequences from "isolate1" and "isolate2" as an example:

    substitution_detection(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')

    this will return "no substitutions in the coding sequence (CDS)" since these sequences are the same in both strains

    now let's mimic a SNP that occurs between the H37Rv and H37Ra carB genes:

    search for "mimic SNP in carB gene" in the "nucleotide_sequence_alignment" function and run code below comment
    to mimic SNP and to show how a SNP that occurs between two gene sequences is identified and displayed in Pandas
    DataFrame

    """

    subseq1, subseq2_coords, subseq2, score, nucleotide_alignment, nucleotide_matrix, isolate1_gene, isolate2_gene = nucleotide_sequence_alignment \
        (pos1, pos2, path, isolate1, isolate2)
    protein_alignment, protein_matrix, isolate1_protein, aligned_protein, isolate2_protein, protein1_amino_acids, protein2_amino_acids = protein_sequence_alignment \
        (pos1, pos2, path, isolate1, isolate2)

    # positions
    # the position of a nucleotide in the alignment and in the genome are different
    # need to get the position from the alignment and minus the number of gaps that have occurred during the alignment to get actual position

    substitution_positions = []
    for i in range(len(isolate1_gene)):
        if nucleotide_matrix[1][i] == 'X' or nucleotide_matrix[1][i] == '.' and nucleotide_matrix[0][i] != '.':
            alignment_symbols = list(nucleotide_matrix[1][:i])
            gaps = alignment_symbols.count('-')
            substitution_positions.append(pos1 + i - gaps)

    # nucleotides at those positions

    substitution_nucleotides_ref = [nucleotide_matrix[0][j] for j in range(len(isolate1_gene)) if nucleotide_matrix[1][j] == 'X' or nucleotide_matrix[1][j] == '.' and nucleotide_matrix[0][j] != '.']
    substitution_nucleotides_alt = [nucleotide_matrix[2][j] for j in range(len(isolate1_gene)) if nucleotide_matrix[1][j] == 'X' or nucleotide_matrix[1][j] == '.' and nucleotide_matrix[0][j] != '.']

    # codon change and amino acid change detection
    # create codons to detect codon changes from mutations

    isolate1_gene = np.ndarray.tolist(isolate1_gene)
    isolate1_gene = ''.join([x for x in isolate1_gene])
    codon_isolate1_gene = [isolate1_gene[i:i + 3] for i in range(0, len(subseq1), 3)]

    isolate2_gene = np.ndarray.tolist(isolate2_gene)
    isolate2_gene = ''.join([x for x in isolate2_gene])
    codon_isolate2_gene = [isolate2_gene[i:i + 3] for i in range(0, len(subseq2), 3)]

    codons = []
    codon_mutations_ref = []
    codon_mutations_alt = []
    amino_acids = []
    for i in range(len(codon_isolate1_gene)):
        for j in range(len(codon_isolate1_gene[i])):
            if codon_isolate1_gene[i][j] != codon_isolate2_gene[i][j]:
                codon_mutations_ref.append(codon_isolate1_gene[i][j])
            codon_mutations_alt.append(codon_isolate2_gene[i][j])
            codons.append(codon_isolate1_gene[i] + '/' + codon_isolate2_gene[i])
            if isolate1_protein[i] != isolate2_protein[i]:
                amino_acids.append(isolate1_protein[i] + '/' + isolate2_protein[i])
            else:
                amino_acids.append('synonymous_coding')

    codons = [x for x in codons if x != []]
    codon_mutations_ref = [x for x in codon_mutations_ref if x != []]
    codon_mutations_alt = [x for x in codon_mutations_alt if x != []]

    substitution_codons = \
        [codons[i] if substitution_nucleotides_ref[i] == codon_mutations_ref[i] and substitution_nucleotides_alt[i] ==
                     codon_mutations_alt[i] else 'No change' for i in range(len(substitution_nucleotides_ref))]

    # group positions and nucleotides of MNPs and append to new lists below - creates list of lists in each empty list below

    sub_pos_range = []
    sub_nucleotide_ref = []
    sub_nucleotide_alt = []
    sub_codons = []
    sub_amino_acids = []
    listPosition = 0
    if len(substitution_positions) == 0:
        print(sub_pos_range)
    elif len(substitution_positions) == 1:
        sub_pos_range.append(substitution_positions[0])
    sub_nucleotide_ref.append(substitution_nucleotides_ref[0])
    sub_nucleotide_alt.append(substitution_nucleotides_alt[0])
    sub_codons.append(substitution_codons[0])
    sub_amino_acids.append(amino_acids[0])
    elif len(substitution_positions) == 2:
    subst_pos_list = []
    subst_ref_nucl_list = []
    subst_alt_nucl_list = []
    subst_codon_list = []
    subst_amino_acid = []
    if (substitution_positions[0] + 1) == (substitution_positions[-1]):
        subst_pos_list.append([substitution_positions[0]])
        subst_pos_list.append([substitution_positions[-1]])
        sub_pos_range.append(subst_pos_list)
        subst_ref_nucl_list.append(substitution_nucleotides_ref[0])
        subst_ref_nucl_list.append(substitution_nucleotides_ref[-1])
        sub_nucleotide_ref.append(subst_ref_nucl_list)
        subst_alt_nucl_list.append(substitution_nucleotides_alt[0])
        subst_alt_nucl_list.append(substitution_nucleotides_alt[-1])
        sub_nucleotide_alt.append(subst_alt_nucl_list)
        subst_codon_list.append(substitution_codons[0])
        subst_codon_list.append(substitution_codons[-1])
        sub_codons.append(subst_codon_list)
        subst_amino_acid.append(amino_acids[0])
        subst_amino_acid.append(amino_acids[-1])
        sub_amino_acids.append(subst_amino_acid)
    else:
        sub_pos_range.append([substitution_positions[0]])
        sub_pos_range.append([substitution_positions[-1]])
        sub_nucleotide_ref.append(substitution_nucleotides_ref[0])
        sub_nucleotide_ref.append(substitution_nucleotides_ref[-1])
        sub_nucleotide_alt.append(substitution_nucleotides_alt[0])
        sub_nucleotide_alt.append(substitution_nucleotides_alt[-1])
        sub_codons.append(substitution_codons[0])
        sub_codons.append(substitution_codons[-1])
        sub_amino_acids.append(amino_acids[0])
        sub_amino_acids.append(amino_acids[-1])

    else:
    while listPosition < len(substitution_positions) - 2:
        originalPosition = listPosition
        tempList_pos = []
        tempList_nucl_ref = []
        tempList_nucl_alt = []
        tempList_codon = []
        tempList_amino_acid = []
        checkIfContinue = True
        while checkIfContinue:
            if (substitution_positions[listPosition] + 1) == substitution_positions[listPosition + 1]:
                if originalPosition == listPosition:
                    tempList_pos.append(substitution_positions[listPosition])
                tempList_pos.append(substitution_positions[listPosition + 1])
                tempList_nucl_ref.append(substitution_nucleotides_ref[listPosition])
                tempList_nucl_ref.append(substitution_nucleotides_ref[listPosition + 1])
                tempList_nucl_alt.append(substitution_nucleotides_alt[listPosition])
                tempList_nucl_alt.append(substitution_nucleotides_alt[listPosition + 1])
                tempList_codon.append(substitution_codons[listPosition])
                tempList_codon.append(substitution_codons[listPosition + 1])
                tempList_amino_acid.append(amino_acids[listPosition])
                tempList_amino_acid.append(amino_acids[listPosition + 1])
                else:
                tempList_pos.append(substitution_positions[listPosition + 1])
                tempList_nucl_ref.append(substitution_nucleotides_ref[listPosition + 1])
                tempList_nucl_alt.append(substitution_nucleotides_alt[listPosition + 1])
                tempList_codon.append(substitution_codons[listPosition + 1])
                tempList_amino_acid.append(amino_acids[listPosition + 1])
                if substitution_positions[listPosition + 1] == substitution_positions[-1]:
                    break
            else:
                if len(tempList_pos) == 0:
                    tempList_pos.append(substitution_positions[listPosition])
                tempList_nucl_ref.append(substitution_nucleotides_ref[listPosition])
                tempList_nucl_alt.append(substitution_nucleotides_alt[listPosition])
                tempList_codon.append(substitution_codons[listPosition])
                tempList_amino_acid.append(amino_acids[listPosition])
                checkIfContinue = False
        listPosition += 1

        sub_pos_range.append(tempList_pos)
        sub_nucleotide_ref.append(tempList_nucl_ref)
        sub_nucleotide_alt.append(tempList_nucl_alt)
        sub_codons.append(tempList_codon)
        sub_amino_acids.append(tempList_amino_acid)

    check = True
    lastPos = -1
    tempLast_pos = []
    tempLast_nucl_ref = []
    tempLast_nucl_alt = []
    tempLast_codons = []
    tempLast_amino_acids = []
    while check:
        if (substitution_positions[lastPos]) == (substitution_positions[lastPos - 1] + 1):
            if -1 == lastPos:
                tempLast_pos.append(substitution_positions[lastPos])
                tempLast_pos.append(substitution_positions[lastPos - 1])
                tempLast_nucl_ref.append(substitution_nucleotides_ref[lastPos])
                tempLast_nucl_ref.append(substitution_nucleotides_ref[lastPos - 1])
                tempLast_nucl_alt.append(substitution_nucleotides_alt[lastPos])
                tempLast_nucl_alt.append(substitution_nucleotides_alt[lastPos - 1])
                tempLast_codons.append(substitution_codons[lastPos])
                tempLast_codons.append(substitution_codons[lastPos - 1])
                tempLast_amino_acids.append(amino_acids[lastPos])
                tempLast_amino_acids.append(amino_acids[lastPos - 1])
            else:
                tempLast_pos.append(substitution_positions[lastPos - 1])
                tempLast_nucl_ref.append(substitution_nucleotides_ref[lastPos - 1])
                tempLast_nucl_alt.append(substitution_nucleotides_alt[lastPos - 1])
                tempLast_codons.append(substitution_codons[lastPos - 1])
                tempLast_amino_acids.append(amino_acids[lastPos - 1])
        else:
            if len(tempLast_pos) == 0:
                sub_pos_range.append([substitution_positions[lastPos - 1]])
                tempLast_pos.append(substitution_positions[lastPos])
                sub_nucleotide_ref.append(substitution_nucleotides_ref[lastPos - 1])
                tempLast_nucl_ref.append(substitution_nucleotides_ref[lastPos])
                sub_nucleotide_alt.append(substitution_nucleotides_alt[lastPos - 1])
                tempLast_nucl_alt.append(substitution_nucleotides_alt[lastPos])
                sub_codons.append([substitution_codons[lastPos - 1]])
                tempLast_codons.append(substitution_codons[lastPos])
                sub_amino_acids.append([amino_acids[lastPos - 1]])
                tempLast_amino_acids.append(amino_acids[lastPos])

        check = False
        lastPos -= 1
    tempLast_pos.reverse()
    tempLast_nucl_ref.reverse()
    tempLast_nucl_alt.reverse()
    tempLast_codons.reverse()
    tempLast_amino_acids.reverse()

    sub_pos_range.append(tempLast_pos)
    sub_nucleotide_ref.append(tempLast_nucl_ref)
    # sub_nucleotide_ref.append(substitution_nucleotides_ref[-1])
    sub_nucleotide_alt.append(tempLast_nucl_alt)
    # sub_nucleotide_alt.append(substitution_nucleotides_alt[-1])
    sub_codons.append(tempLast_codons)
    sub_amino_acids.append(tempLast_amino_acids)

    if len(sub_pos_range) == 0 or len(sub_pos_range) == 1:
        print(sub_pos_range)
    else:
        if sub_pos_range[-2] == sub_pos_range[-1]:
            sub_pos_range.pop(-1)

    for i in range(len(sub_codons)):
        sub_codons[i] = list(set(sub_codons[i]))

    for i in range(len(sub_amino_acids)):
        sub_amino_acids[i] = list(set(sub_amino_acids[i]))

    # joining together positions and the nucleotides associated with those positions
    # the list comprehensions below joins the reference and alternate nucleotides for the positions ranges together to report them as a single nucleotide string

    sub_nucleotide_ref = [''.join(sub_nucleotide_ref[i][0:]) if len(sub_nucleotide_ref[i]) != 1 else ", ".join(
        map(str, sub_nucleotide_ref[i])) for i in range(len(sub_nucleotide_ref))]
    sub_nucleotide_alt = [''.join(sub_nucleotide_alt[i][0:]) if len(sub_nucleotide_alt[i]) != 1 else ", ".join(
        map(str, sub_nucleotide_alt[i])) for i in range(len(sub_nucleotide_alt))]

    if len(sub_pos_range) > 0:
        substitutions_data = {'chromosome': 1,
                              'positions (isolate1)': sub_pos_range,
                              'reference allele': sub_nucleotide_ref,
                              'alternate allele': sub_nucleotide_alt,
                              'mutation type': 'substitution',
                              'frameshift': '-',
                              'old codon/new codon': sub_codons,
                              'old AA/new AA': sub_amino_acids}
    df_substitutions = pd.DataFrame(substitutions_data)
    return df_substitutions
    else:
    return 'no substitutions in this CDS'


def insertion_detection(pos1, pos2, path, isolate1, isolate2):
    """uses the nucleotide alignment to detect insertions and reports them in a Pandas dataframe in VCF format
    :param pos1: start position of gene in isolate1
    :param pos2: end position of gene in isolate1
    :param path: path to genome graph contained within XML file containing isolate1 and isolate2
    :param isolate1: the first strain from the genome graph being compared(i.e. H37Rv in this case since it is being used as a reference)
    :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
    :return:
    df_insertions: Pandas DataFrame reporting any insertions that occurred between the "isolate1" and "isolate2"
    versions of the gene

    NOTES
    ------------------

    This function allows for the detection of insertions that occur between similar gene sequences extracted from a genome graph
    However, this function can only compare two similar sequences and thus in order to understand the biological significance
    of the insertions that may occur, the ability to perform a multiple sequence alignment within the function will be required

    positions of insertions are reported relative to "isolate2" in the Pandas DataFrame since insertions will occur across
    a position range in "isolate2" and not "isolate1"
    It is important to note that InDels are identified relative to "isolate1"
    (i.e. a deletion is where a nucleotide from "isolate1" is removed in "isolate2" and an insertion is where a nucleotide is
    added to "isolate2" that is not present in "isolate1" - this is determined by looking at the nucleotide pairwise alignment)

    EXAMPLE
    ------------------

    Using the carB gene sequences from "isolate1" and "isolate2" as an example:

    insertion_detection(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')

    this will return "no insertions in the coding sequence (CDS)" since these sequences are the same in both strains

    now let's mimic a three-nucleotide insertion that occurs between the H37Rv and H37Ra carB genes:

    search for "mimic insertion in carB gene" in the "nucleotide_sequence_alignment" function and run code below comment
    to mimic the insertion to show how an insertion that occurs between two gene sequences is identified and
    displayed in Pandas DataFrame
    """

    subseq1, subseq2_coords, subseq2, score, nucleotide_alignment, nucleotide_matrix, isolate1_gene, isolate2_gene = nucleotide_sequence_alignment(
        pos1, pos2, path, isolate1, isolate2)

    # determine positions where insertions occurred

    insertion_positions = []
    for i in range(len(isolate1_gene)):
        if nucleotide_matrix[0][i] == '-' or nucleotide_matrix[0][i] == '.':
            alignment_symbols = list(nucleotide_matrix[0][:i])
            gaps = alignment_symbols.count('-')
            insertion_positions.append(subseq2_coords[0] + i)  # get position of insertion in "isolate2"

    # nucleotides that have been inserted

    insertion_nucleotides = [nucleotide_matrix[2][j] for j in range(len(isolate1_gene)) if
                             nucleotide_matrix[0][j] == '-' or nucleotide_matrix[0][j] == '.']

    # group consecutive positions and their respective nucleotides for insertions > 1 nucleotide

    insert_pos_range = []
    insert_nucleotide = []
    listPosition = 0
    if len(insertion_positions) == 0:
        print(insert_pos_range)
    elif len(insertion_positions) == 1:
        insert_pos_range.append(insertion_positions[0])
    insert_nucleotide.append(insertion_nucleotides[0])
    elif len(insertion_positions) == 2:
    ins_list = []
    if (insertion_positions[0] + 1) == (insertion_positions[-1]):
        ins_list.append([insertion_positions[0]])
        ins_list.append([insertion_positions[-1]])
        insert_pos_range.append(ins_list)
    else:
        insert_pos_range.append([insertion_positions[0]])
        insert_pos_range.append([insertion_positions[-1]])
    else:
    while listPosition < len(insertion_positions) - 2:
        originalPosition = listPosition
        tempList_pos = []
        tempList_nucl = []
        checkIfContinue = True
        while checkIfContinue:
            if (insertion_positions[listPosition] + 1) == insertion_positions[listPosition + 1]:
                if originalPosition == listPosition:
                    tempList_pos.append(insertion_positions[listPosition])
                tempList_pos.append(insertion_positions[listPosition + 1])
                tempList_nucl.append(insertion_nucleotides[listPosition])
                tempList_nucl.append(insertion_nucleotides[listPosition + 1])
                else:
                tempList_pos.append(insertion_positions[listPosition + 1])
                tempList_nucl.append(insertion_nucleotides[listPosition + 1])
                if insertion_positions[listPosition + 1] == insertion_positions[-1]:
                    break
            else:
                if len(tempList_pos) == 0:
                    tempList_pos.append(insertion_positions[listPosition])
                tempList_nucl.append(insertion_nucleotides[listPosition])
                checkIfContinue = False
        listPosition += 1

        insert_pos_range.append(tempList_pos)
        insert_nucleotide.append(tempList_nucl)

    check = True
    lastPos = -1
    tempLast_pos = []
    tempLast_nucl = []
    if len(insert_pos_range[0]) == len(insertion_positions):
        check = False
    while check:
        if (insertion_positions[lastPos]) == (insertion_positions[lastPos - 1] + 1):
            if -1 == lastPos:
                tempLast_pos.append(insertion_positions[lastPos])
                tempLast_pos.append(insertion_positions[lastPos - 1])
                tempLast_nucl.append(insertion_nucleotides[lastPos])
                tempLast_nucl.append(insertion_nucleotides[lastPos - 1])
            else:
                tempLast_pos.append(insertion_positions[lastPos - 1])
                tempLast_nucl.append(insertion_nucleotides[lastPos - 1])
        else:
            if len(tempLast_pos) == 0:
                insert_pos_range.append(insertion_positions[lastPos - 1])
                tempLast_pos.append(insertion_positions[lastPos])
                insert_nucleotide.append(insertion_nucleotides[lastPos - 1])
                tempLast_nucl.append(insertion_nucleotides[lastPos])
        check = False
        lastPos -= 1
    tempLast_pos.reverse()
    tempLast_nucl.reverse()

    if len(tempLast_pos) > 0:
        insert_pos_range.append(tempLast_pos)
        insert_nucleotide.append(tempLast_nucl)

    insert_nucleotide = [''.join(insert_nucleotide[i][0:]) if len(insert_nucleotide[i]) != 1 else ", ".join(
        map(str, insert_nucleotide[i])) for i in range(len(insert_nucleotide))]

    print(insert_pos_range)
    print(insert_nucleotide)

    # detect whether insertion is a frameshift mutation or not

    frameshift_insertions = []
    for i in range(len(insert_nucleotide)):
        if len(insert_nucleotide[i]) % 3 != 0:
            frameshift_insertions.append('frameshift!')
        else:
            frameshift_insertions.append('-')

    # represent results in a pandas dataframe

    if len(insert_pos_range) > 0:
        insertion_data = {'chromosome': 1,
                          'positions (isolate2)': insert_pos_range,
                          'reference allele': '-',
                          'alternate allele': insert_nucleotide,
                          'mutation type': 'insertion',
                          'frameshift': frameshift_insertions,
                          }
    df_insertions = pd.DataFrame(insertion_data)
    return df_insertions
    else:
    return 'there are no insertions in the CDS'


def deletion_detection(pos1, pos2, path, isolate1, isolate2):
    """uses the nucleotide alignment to detect deletions and reports them in a Pandas dataframe in VCF format
    :param pos1: start position of gene in isolate1
    :param pos2: end position of gene in isolate1
    :param path: path to genome graph contained within XML file containing isolate1 and isolate2
    :param isolate1: the first strain from the genome graph being compared(i.e. H37Rv in this case since it is being used as a reference)
    :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
    :return:
    df_deletions: Pandas dataframe reporting any deletions that occurred between the "isolate1" and "isolate2"
    versions of the gene

    NOTES
    -------------------

    This function allows for the detection of deletions that occur between similar gene sequences extracted from a genome graph
    However, this function can only compare two similar sequences and thus in order to understand the biological significance
    of the deletions that may occur, the ability to perform a multiple sequence alignment within the function will be required

    positions of deletions are reported relative to "isolate1" in the Pandas DataFrame
    (i.e. a deletion is where a nucleotide from "isolate1" is removed in "isolate2" and an insertion is where a nucleotide is
    added to "isolate2" that is not present in "isolate1" - this is determined by looking at the nucleotide pairwise alignment)

    It is important to note that InDels are identified relative to "isolate1" but that deletions are reported relative to "isolate1"
    and insertions are reported relative "isolate2"


    EXAMPLE
    -------------------

    Using the carB gene sequences from "isolate1" and "isolate2" as an example:

    insertion_detection(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')

    this will return "no insertions in the coding sequence (CDS)" since these sequences are the same in both strains

    now let's mimic a three-nucleotide deletion that occurs in the H37Ra carB gene:

    search for "mimic deletion in carB gene" in the "nucleotide_sequence_alignment" function and use code below comment
    to mimic the deletion. Run code to show how an deletion that occurs between two gene sequences is identified and
    displayed in Pandas DataFrame

    """

    subseq1, subseq2_coords, subseq2, score, nucleotide_alignment, nucleotide_matrix, isolate1_gene, isolate2_gene = nucleotide_sequence_alignment(
        pos1, pos2, path, isolate1, isolate2)

    # positions in CDS where deletions occur

    deletion_positions = []
    for i in range(len(isolate1_gene)):
        if nucleotide_matrix[2][i] == '-' or nucleotide_matrix[2][i] == '.':
            alignment_symbols = list(nucleotide_matrix[0][:i])
            gaps = alignment_symbols.count('-')
            deletion_positions.append(pos1 + i - gaps)

    # deleted nucleotides

    deletion_nucleotides = [nucleotide_matrix[0][j] for j in range(len(isolate1_gene)) if
                            nucleotide_matrix[2][j] == '-' or nucleotide_matrix[2][j] == '.']

    # group MN deletion positions and the nucleotides associated with them (MN = multiple nucleotide)

    del_pos_range = []
    del_nucleotide = []
    listPosition = 0
    if len(deletion_positions) == 0:
        print(del_pos_range)
    elif len(deletion_positions) == 1:
        del_pos_range.append(deletion_positions[0])
    del_nucleotide.append(deletion_nucleotides[0])
    elif len(deletion_positions) == 2:
    del_list = []
    if (deletion_positions[0] + 1) == (deletion_positions[-1]):
        del_list.append([deletion_positions[0]])
        del_list.append([deletion_positions[-1]])
        del_pos_range.append(del_list)
    else:
        del_pos_range.append([deletion_positions[0]])
        del_pos_range.append([deletion_positions[-1]])
    else:
    while listPosition < len(deletion_positions) - 2:
        originalPosition = listPosition
        tempList_pos = []
        tempList_nucl = []
        checkIfContinue = True
        while checkIfContinue:
            if (deletion_positions[listPosition] + 1) == deletion_positions[listPosition + 1]:
                if originalPosition == listPosition:
                    tempList_pos.append(deletion_positions[listPosition])
                tempList_pos.append(deletion_positions[listPosition + 1])
                tempList_nucl.append(deletion_nucleotides[listPosition])
                tempList_nucl.append(deletion_nucleotides[listPosition + 1])
                else:
                tempList_pos.append(deletion_positions[listPosition + 1])
                tempList_nucl.append(deletion_nucleotides[listPosition + 1])
                if deletion_positions[listPosition + 1] == deletion_positions[-1]:
                    break
            else:
                if len(tempList_pos) == 0:
                    tempList_pos.append(deletion_positions[listPosition])
                tempList_nucl.append(deletion_nucleotides[listPosition])
                checkIfContinue = False
        listPosition += 1

        del_pos_range.append(tempList_pos)
        del_nucleotide.append(tempList_nucl)

    check = True
    lastPos = -1
    tempLast_pos = []
    tempLast_nucl = []
    if len(del_pos_range[0]) == len(deletion_positions):
        check = False
    else:
        while check:
            if (deletion_positions[lastPos]) == (deletion_positions[lastPos - 1] + 1):
                if -1 == lastPos:
                    tempLast_pos.append(deletion_positions[lastPos])
                tempLast_pos.append(deletion_positions[lastPos - 1])
                tempLast_nucl.append(deletion_nucleotides[lastPos])
                tempLast_nucl.append(deletion_nucleotides[lastPos - 1])
                else:
                tempLast_pos.append(deletion_positions[lastPos - 1])
                tempLast_nucl.append(deletion_nucleotides[lastPos - 1])
            else:
                if len(tempLast_pos) == 0:
                    del_pos_range.append(deletion_positions[lastPos - 1])
                tempLast_pos.append(deletion_positions[lastPos])
                del_nucleotide.append(deletion_nucleotides[lastPos - 1])
                tempLast_nucl.append(deletion_nucleotides[lastPos])
                check = False
        lastPos -= 1
        tempLast_pos.reverse()
        tempLast_nucl.reverse()

        if len(tempLast_pos) > 0:
            del_pos_range.append(tempLast_pos)
        del_nucleotide.append(tempLast_nucl)

        del_pos_range.append(tempLast_pos)
        del_nucleotide.append(tempLast_nucl)

    if len(del_pos_range) == 0 or len(del_pos_range) == 1:
        print(del_pos_range)
    else:
        if del_pos_range[-2] == del_pos_range[-1]:
            del_pos_range.pop(-1)

    # concatenate nucleotides together for MN deletions

    deletion_nucleotides = [
        ''.join(del_nucleotide[i][0:]) if len(del_pos_range[i]) != 1 else ", ".join(map(str, del_nucleotide[i])) for i
        in range(len(del_pos_range))]

    # determine whether the small/large deletion causes a frameshift or not

    frameshift_deletions = []
    for i in range(len(del_pos_range)):
        if len(del_pos_range[i]) % 3 != 0:
            frameshift_deletions.append('frameshift!')
        else:
            frameshift_deletions.append('-')

    # code to join together nucleotides of MN deletions
    for i in range(len(del_pos_range)):
        if len(del_pos_range[i]) != 1:
            ",".join(map(str, del_pos_range[i]))
            ''.join(del_nucleotide[i][0:])
        else:
            ", ".join(map(str, del_pos_range[i]))
            ", ".join(map(str, del_nucleotide[i]))

    # represent results in a pandas dataframe

    if len(del_pos_range) > 0:
        deletion_data = {'chromosome': 1,
                         'positions (isolate1)': del_pos_range,
                         'reference allele': deletion_nucleotides,
                         'alternate allele': '-',
                         'mutation type': 'deletion',
                         'frameshift': frameshift_deletions,
                         }
    df_deletions = pd.DataFrame(deletion_data)
    return df_deletions
    else:
    return 'there are no deletions in the CDS'


# SCORING MATRIX FOR CORE, ACCESSORY AND UNIQUE GENES

# take the two proteins and compares amino acid sequences of the two

# 'ancestral_protein' and 'derived_protein' previously assigned so use these two variables for matrix

def scoring_matrix(pos1, pos2, path, isolate1, isolate2):
    """gene classification function that uses protein sequence similarity of the same gene from two different isolates
    to determine whether a gene is core, accessory or unique
    :param pos1: start position of gene in isolate1
    :param pos2: end position of gene in isolate1
    :param path: path to genome graph contained within XML file containing isolate1 and isolate2
    :param isolate1: the first strain from the genome graph being compared(i.e. H37Rv in this case since it is being used as a reference)
    :param isolate2: the second strain from the genome graph being compared (i.e. H37Ra strain)
    :return:
    score: the amino acid sequence similarity score of the two protein sequences

    NOTES
    --------------------

    Gene classification is based on protein sequence similarity score - this function will return whether the gene being assessed
    is likely to be characterised as either a "core gene", "accessory gene" or "unique gene" using sequence similarity cut-off values

    Since only two sequences can be compared at once within the function, gene classification will not be accurate since
    the definitions of core genes (present in all strains), accessory genes (present in more than two strains) and unique
    genes (specific to a certain strain) require that multiple similar genes from different isolates be compared within a
    multiple sequence alignment. Adding this capability to the function will allow for more accurate classification of genes

    The scoring system of this function is simple (match = 1, mismatch/gap = -1) and can be changed to suit needs of user.
    The scoring system has a limitation in that but is not able to discern conservative amino acid changes that occur
    which would lead to a mismatch occurring instead of a match. Gaps and mismatches are also given the same score which may
    be lead to more inaccuracy of the scoring matrix and therefore these could be improved upon.


    EXAMPLE
    --------------------

    Using the carB gene sequences from "isolate1" and "isolate2" as an example:

    scoring_matrix(1557101, 1560448, './H37R_pangenome.xml', 'H37Rv', 'H37Ra')

    since these genes are very similar, this gene may be essential to the survival of these M. tuberculosis strains
    and may therefore be considered a core gene - this prediction would become more accurate once multiple gene sequences
    can be compared (may be classified an accessory gene after multiple sequence alignment)

    """

    protein_alignment, protein_matrix, isolate1_protein, aligned_protein, isolate2_protein, protein1_amino_acids, protein2_amino_acids = protein_sequence_alignment(
        pos1, pos2, path, isolate1, isolate2)

    score = []
    for i in range(len(isolate1_protein)):
        if isolate1_protein[i] == isolate2_protein[i]:
            score.append(1)
        else:
            score.append(-1)
    score = 100 * (sum(score) / len(isolate1_protein))
    if 95 <= score <= 100:
        print(score)
    print('this gene is possibly a core gene')
    elif 90 <= score < 95:
    print(score)
    print('this gene is possibly an accessory gene')
    else:
    print(score)
    print('this gene is possibly a unique gene')

# -----------------------------------------------------Sequece_Homology Functions end here
