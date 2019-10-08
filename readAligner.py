from gengraph import *
from gengraphTool import *
import networkx as nx
import matplotlib.pyplot as plt
import difflib


class Aligner(nx.DiGraph):

    kmerDict = {}
    global firstPath
    firstPath = True

    def fast_kmer_create(self, kmerLength):
        '''

        :param kmerLength: A number that indicated the length of the kmers that you want as an output
        :return: A dictionary with all the kmers in the graph, each entry in the dictionary has a kmer name, the nodes that the kmer covers, the sequence of the kmer and
                the starting position of the kmer on the graph
        '''
        # Start Positions stored in the dictionary are character list positions, so the first character in the node sequence is position 0 not position 1
        graphNodes = list(self.nodes)
        kmerNumber = 1

        #Runs through the graph nodes creating kmers and checking if each of the nodes has an inversion in it
        for i in graphNodes:
            currentNodeName = str(i)
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
            if inv == True:
                for j in range(len(nodeSequence)):
                    if kmerStartPos + kmerLength <= len(nodeSequence):
                        kmerSeq = nodeSequence[kmerStartPos:kmerStartPos + kmerLength]
                        kmerRevSeq = revNodeSequence[kmerStartPos:kmerStartPos + kmerLength]
                        dictStartPos = kmerStartPos
                        kmerStartPos += 1
                        tempDict = {'nodes': [currentNodeName], 'sequence': kmerSeq, 'inversesequence':kmerRevSeq, 'startposition': dictStartPos}
                        newkmer = {'kmer_' + str(kmerNumber): tempDict}
                        kmerNumber += 1
                        Aligner.kmerDict.update(newkmer)
                    else:
                        kmerPartseq = nodeSequence[kmerStartPos:]
                        dictStartPos = kmerStartPos
                        kmerStartPos += 1
                        currentNode = [currentNodeName]
                        self.kmer_over_multiple_nodes(kmerPartseq, currentNode, kmerLength, kmerNumber, dictStartPos,True)
                        kmerNumber += 1
            else:
                for j in range(len(nodeSequence)):
                    if kmerStartPos + kmerLength <= len(nodeSequence):
                        kmerSeq = nodeSequence[kmerStartPos:kmerStartPos+kmerLength]
                        dictStartPos = kmerStartPos
                        kmerStartPos += 1
                        tempDict = {'nodes': [currentNodeName], 'sequence': kmerSeq, 'startposition': dictStartPos}
                        newkmer = {'kmer_' + str(kmerNumber): tempDict}
                        kmerNumber += 1
                        Aligner.kmerDict.update(newkmer)
                    else:
                        kmerPartseq = nodeSequence[kmerStartPos:]
                        dictStartPos = kmerStartPos
                        kmerStartPos += 1
                        currentNode = [currentNodeName]
                        self.kmer_over_multiple_nodes(kmerPartseq, currentNode, kmerLength, kmerNumber,dictStartPos,False)
                        kmerNumber += 1
        return Aligner.kmerDict

    def kmer_over_multiple_nodes(self, kmerPart, currentNode, kmerLength, kmerNumber,kmerStartPos,hasInv):
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
                    tempDict = {'nodes': allNodes, 'sequence': kmerSeq, 'inversesequence':kmerRevSeq, 'startposition': kmerStartPos}
                    kmerKey = 'kmer_' + str(kmerNumber) + '_' + str(k + 1)
                    if kmerKey in Aligner.kmerDict.keys():
                        lastDigit = kmerKey[-1]
                        lastDigit = int(lastDigit) + 1
                        newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
                        kmerKey = newkmerKey
                    newkmer = {kmerKey: tempDict}
                    Aligner.kmerDict.update(newkmer)
                else:
                    #If the node is not longenough to make a full kmer you will need to recursively move onto its neighbours until
                    #all of the kmers have been successfully created
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
                        tempDict = {'nodes': allNodes, 'sequence': kmerSeq, 'inversesequence':revKmerSeq, 'startposition': kmerStartPos}
                        kmerKey = 'kmer_' + str(kmerNumber) + '_' + str(k + 1)
                        if kmerKey in Aligner.kmerDict.keys():
                            lastDigit = kmerKey[-1]
                            lastDigit = int(lastDigit) + 1
                            newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
                            kmerKey = newkmerKey
                        newkmer = {kmerKey: tempDict}
                        Aligner.kmerDict.update(newkmer)
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
                        if kmerKey in Aligner.kmerDict.keys():
                            lastDigit = kmerKey[-1]
                            lastDigit = int(lastDigit) + 1
                            newkmerKey = kmerKey[:len(kmerKey) - 1] + str(lastDigit)
                            kmerKey = newkmerKey
                        newkmer = {kmerKey: tempDict}
                        Aligner.kmerDict.update(newkmer)
            else:
                kmerTemp = kmerPart + nextNodeSeq
                nodesInvolved = []
                posx = 0
                for x in range(len(currentNode)):
                    nodesInvolved.append(currentNode[posx])
                    posx += 1
                nodesInvolved.append(connectedNodes[k])
                #print(kmerTemp + " " +str(nodesInvolved))
                self.kmer_over_multiple_nodes(kmerTemp, nodesInvolved,kmerLength,kmerNumber,kmerStartPos,False)
                # recursively goes over multiple nodes
                # will need the currentNode that is parsed at this step so instead you have the two previous nodes the kmer goes over so you have all three nodes saved into the dictionary

    def kmers_of_node(self, nodeName):
        '''

        :param nodeName: the Node that you want to use create kmers from
        :return: returns a dictionary of all of the kmers involved with the given node.
        '''
        kmersOfNode = {}
        if len(Aligner.kmerDict) == 0:
            #default kmer length is 4
            Aligner.fast_kmer_create(self, 4)

        keyList = list(Aligner.kmerDict.keys())
        #creates the kmers of the current node if the dictionary of all of the kmers is not already made
        for i in keyList:
            currentKmer = Aligner.kmerDict.get(i)
            kmerNodeList = currentKmer['nodes']
            for j in kmerNodeList:
                if j == nodeName:
                    kmer = {i: currentKmer}
                    kmersOfNode.update(kmer)
        return kmersOfNode

    def create_query_kmers(self, querySequence, kmerLength):
        '''

        :param querySequence: the sequence of the query you want to align to the graph
        :param kmerLength: The length of the kmers that you want to make from the query sequence.
                These should be the same as the graph kmer length
        :return:Returns a dictionary with all of the kmers from the query sequence. The dictionary contains the query
                kmer name and query kmer sequence.
        '''
        startPos = 0
        queryKmers = {}
        kmerNumber = 1
        for i in range(len(querySequence)):
            if startPos + kmerLength <=len(querySequence):
                kmerSeq = querySequence[startPos:startPos+kmerLength]
                kmer = {'qkmer_' + str(kmerNumber): kmerSeq}
                kmerNumber += 1
                queryKmers.update(kmer)
                startPos += 1
        return queryKmers

    # need to fix readgaps. It is sometimes still not giving the correct values for gaps in between
    def readGaps(self,nodeNeighbours, query, firstReadInfo, secondReadInfo, distanceBetween):
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
                    self.dist = distanceBetween + ((len(self.nodes[firstReadInfo['nodescoveredbyread'][-1]]['sequence'])) - (firstReadInfo['alignementendpositioninlastnode'] + 1)) + secondReadInfo['alignmentstartpositioninfirstnode']
                    self.firstPath = False
                    #print(self.dist)
                    return self.dist
                else:
                    return self.dist
        for i in nodeNeighbours:
            if len(self.nodes[i]['sequence']) + int(distanceBetween) < len(query):
                distanceBetween = int(distanceBetween) + len(self.nodes[i]['sequence'])
                newNeighbours = list(self[i])
                self.readGaps(newNeighbours, query,firstReadInfo, secondReadInfo, distanceBetween)
            else:
                return distanceBetween
        return distanceBetween


    def debruin_read_alignment(self , query, kmerLength):
        '''
        Running this Function will print out blocks of information on the aligned reads. Each block of 7 lines will correspond to a single aligned read.
        Line 1: The Aligned read sequence
        Line 2: The Nodes that the aligned read covers
        Line 3: All nodes that the read aligns to inversely(Nodes with an aligned inversion)
        Line 4: The starting alignment position on the first node involved in the alignment
        Line 5: The ending alignment position on the last node involved in the alignment
        Line 6: Percentage of the initial query that has aligned to the graph
        Line 7:Percentage of the aligned read that has successfully aligned to the graph
        :param query: This is the read, as a string, that you want to try and align to the graph.
        :param kmerLength: The length of the kmers that will be created
        :return: returns a dictionary with the sequences that align, x's represent sections of the sequence that doesnt correctly align
                It also returns the nodes that the sequence aligns to and the position on the first node the sequence starts aligning to and the position on the last node that the sequence ends aligning to.
                These positions are in string positions, so the first base pair is position 0.
        '''
        queryKmerDict = self.create_query_kmers(query, kmerLength)
        referenceKmerDict = self.fast_kmer_create(kmerLength)
        # creates the query and graph kmers with the desired kmer length
        finalkmerGroups = []
        groupedListOftwos = []
        matchedKmers = []
        finalAlignedReads = {}
        inversionNodes = []
        #list of lists of matched kmers
        #first pos is query kmer name, second pos is reference nodes that the kmer covers, third is the start position of reference kmer, fourth is reference kmer name, 5th is qkmerseq

        queryKmerKeys = list(queryKmerDict.keys())
        referenceKmerKeys = list(referenceKmerDict.keys())

        #mathces the each query kmer to refernce graph kmers with the same sequence
        for i in queryKmerKeys:
            currentQueryKmer = queryKmerDict.get(i)
            for j in referenceKmerKeys:
                currentReferenceKmer = referenceKmerDict.get(j)
                if currentReferenceKmer['sequence'] == currentQueryKmer:
                    tempList = [i, currentReferenceKmer['nodes'], currentReferenceKmer['startposition'], j, currentQueryKmer,'normal']
                    matchedKmers.append(tempList)
                elif 'inversesequence' in currentReferenceKmer:
                    if currentReferenceKmer['inversesequence'] == currentQueryKmer:
                        tempList = [i, currentReferenceKmer['nodes'], currentReferenceKmer['startposition'], j, currentQueryKmer, 'inverse']
                        matchedKmers.append(tempList)
        #results in a list of list with the matched kmers
        #has each match between query kmer and reference graph kmer

        #runs through the matched kmers and groups the kmers into groups of twos by the kmers that line up next to each other
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
                            node = referenceKmerDict[addOne[addOne.find('-')+1:]]['nodes'][0]
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
                                if(len(k[1]) == 1 and k[5] == 'inverse'):
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

        #takes the groups of two lists and groups all o the kmers that are in a line and next to each other on the graph
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
        #result in final groups with all kmers that are next to each other on the graphs grouped up
        AlignNumber = 1

        #creates a dictionary with all the final aligned kmers that align 100% perfectly to the reference graph
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
                    #shouldnt be finalkmerName, got to get each kmerName
                    if r[0] == qkmerName and kmerName == r[3]:
                        kmerseq = r[4]
                kmernode = referenceKmerDict[kmerName]['nodes']
                nodesCovered.extend(kmernode)
                positionsLeft = 0
                nodesCovered = list(dict.fromkeys(nodesCovered))
                if m == n[0]:
                    startpos = referenceKmerDict[kmerName]['startposition']
                if m == n[-1]:
                    if len(referenceKmerDict[kmerName]['nodes']) == 1 and referenceKmerDict[kmerName]['nodes'][-1] == finalNode:
                        endpos = referenceKmerDict[kmerName]['startposition'] + kmerLength -1
                    else:
                        for b in (referenceKmerDict[kmerName]['nodes']):
                            if b == referenceKmerDict[kmerName]['nodes'][0]:
                                tempnumber = len(self.nodes[b]['sequence']) - referenceKmerDict[kmerName]['startposition']
                                positionsLeft += tempnumber
                            elif b != referenceKmerDict[kmerName]['nodes'][0] and b != referenceKmerDict[kmerName]['nodes'][-1]:
                                positionsLeft += len(self.nodes[b]['sequence'])
                        endpos = len(referenceKmerDict[kmerName]['sequence']) - positionsLeft -1
                if len(matchedSquence) == 0:
                    matchedSquence = matchedSquence + kmerseq
                else:
                    matchedSquence = matchedSquence + kmerseq[-1]
            match = False
            for s in finalAlignedReads:
                if matchedSquence in str(finalAlignedReads[s]['sequence']) and set(nodesCovered).issubset(finalAlignedReads[s]['nodescoveredbyread']):
                    match = True
            if match == False:
                tempValues = {'sequence': matchedSquence, 'nodescoveredbyread': nodesCovered, 'alignmentstartpositioninfirstnode': startpos, 'alignementendpositioninlastnode':endpos}
                tempAlign = {'Aligned_Read_' + str(AlignNumber): tempValues}
                AlignNumber += 1
                finalAlignedReads.update(tempAlign)
        # At the moment this returns reads that matched that are any size larger than one kmer length.
        # finalAlignedReads Are all reads that map 100 percent accurately
        # Start grouping the final aligned read dictionary to show bigger reads and not only 100% alignment
        finalGroupedReads = {}
        readGroups = []
        unlinkedReads = []

        #then finds the distances between the 100 percent mapped reads and if those distances are small enough it groups those 100 percent mapped reads together
        #this represents the final aligned reads which may not be 100 percent accurately mapped to each other.
        for z in finalAlignedReads:
            endNode = finalAlignedReads[z]['nodescoveredbyread'][-1]
            links = 0
            for c in finalAlignedReads:
                distanceBetween = 0
                startNode = finalAlignedReads[c]['nodescoveredbyread'][0]
                if z!=c and list(finalAlignedReads.keys()).index(z) < list(finalAlignedReads.keys()).index(c):
                    if int(z[13:]) + 1 == int(c[13:]):
                        if endNode == startNode:
                            distanceBetween = finalAlignedReads[c]['alignmentstartpositioninfirstnode'] - 1 - finalAlignedReads[z]['alignementendpositioninlastnode']
                        else:
                            if nx.has_path(self,endNode,startNode):
                                path = nx.all_shortest_paths(self,endNode,startNode)
                                pathList = list(path)
                                for v in range(len(pathList)):
                                    tdistanceBetween = 0
                                    pointer = v
                                    startDist = 0
                                    for x in pathList[v]:
                                        if x == pathList[v][0]:
                                            startDist = len(self.nodes[x]['sequence']) - finalAlignedReads[z]['alignementendpositioninlastnode'] -1
                                            tdistanceBetween = tdistanceBetween + startDist
                                        elif x == pathList[v][-1]:
                                            tdistanceBetween = tdistanceBetween + finalAlignedReads[c]['alignmentstartpositioninfirstnode']
                                        else:
                                            tdistanceBetween = tdistanceBetween + len(self.nodes[x]['sequence'])
                                    if distanceBetween == 0:
                                        distanceBetween = tdistanceBetween
                                    elif tdistanceBetween < distanceBetween and distanceBetween != 0:
                                        distanceBetween = tdistanceBetween
                        if distanceBetween < len(query) and distanceBetween != 0:
                            readGroups.append([z, distanceBetween, c])
                            links += 1
            if links == 0:
                contains = False
                for o in readGroups:
                    if o[0] == z or o[2] == z:
                        contains = True
                if contains == False:
                    unlinkedReads.append(z)

        #Alignement read groups will shop the start and end point of the reads but will hva egaps in between that can be larger then gaps between the sections in the reads
        finalReadGroups = []
        intfinalReadGroups = []
        tempReadGroups = []

        for u in readGroups:
            temp =[u[0], u[-1]]
            tempReadGroups.append(temp)
        #print(tempReadGroups)
        out = []
        while len(tempReadGroups) > 0:
            first, *rest = tempReadGroups
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
        #print(intfinalReadGroups)

        for x in intfinalReadGroups:
            readsWithDist = []
            for c in range(len(x)-1):
                for z in readGroups:
                    if z[0] == x[c] and z[-1] == x[c+1]:
                        if c == 0:
                            readsWithDist.append(x[c])
                            readsWithDist.append(z[1])
                            readsWithDist.append(x[c+1])
                        else:
                            readsWithDist.append(z[1])
                            readsWithDist.append(x[c+1])
            finalReadGroups.append(readsWithDist)

        #all the final reads are now grouped together with the distances between the reads represented in the list of reads
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
                elif i%2 == 0:
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
                nodeNos.append(str(r[str(r).find('_') +1:str(r).rfind('_')]) + str(r[str(r).rfind('_') + 1:]))
            for a in nodeNos:
                for b in nodesCoveredByRead:
                    if a == str(b[str(b).find('_') +1:str(b).rfind('_')]) + str(b[str(b).rfind('_') + 1:]):
                        newNodesCoveredByRead.append(b)
            startAlignPos +=1
            endAlignPos += 1
            if newNodesCoveredByRead[-1] in inversionNodes:
                endAlignPos = endAlignPos * -1
            if newNodesCoveredByRead[0] in inversionNodes:
                startAlignPos = startAlignPos * -1
            percentOfQuery = difflib.SequenceMatcher(None, query, checkPercentSeq).ratio() * 100
            percentAligned = (len(sequence)-xCounter)/len(sequence)*100
            tempVal = {'sequence':sequence, 'nodescoveredbyread':newNodesCoveredByRead, 'alignmentstartpositioninfirstnode':startAlignPos, 'alignementendpositioninlastnode':endAlignPos,
                       'percentageofqueryaligned':percentOfQuery,'percentageofalignedreadtograph':percentAligned}
            key = {'Alignment_'+str(AlignNo):tempVal}
            AlignNo += 1
            allAlignedReads.update(key)

        for w in unlinkedReads:
            read = finalAlignedReads.get(w)
            percentOfQuery = difflib.SequenceMatcher(None, query, read['sequence']).ratio() * 100
            percentAligned = (len(read['sequence'])) / len(read['sequence']) * 100
            s = read['alignmentstartpositioninfirstnode']
            e = read['alignementendpositioninlastnode']
            s += 1
            e += 1
            if read['nodescoveredbyread'][-1] in inversionNodes:
                e = e * -1
            if read['nodescoveredbyread'][0] in inversionNodes:
                s = s * -1
            tempVals = {'sequence':read['sequence'], 'nodescoveredbyread':read['nodescoveredbyread'], 'alignmentstartpositioninfirstnode':s, 'alignementendpositioninlastnode':e,
                       'percentageofqueryaligned':percentOfQuery,'percentageofalignedreadtograph':percentAligned}
            key = {'Alignment_' + str(AlignNo): tempVals}
            AlignNo += 1
            allAlignedReads.update(key)

        #print out the blocks of information for each aligned read
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
                print(inversionNodes)
            print(allAlignedReads[t]['alignmentstartpositioninfirstnode'])
            print(allAlignedReads[t]['alignementendpositioninlastnode'])
            print(allAlignedReads[t]['percentageofqueryaligned'])
            print(allAlignedReads[t]['percentageofalignedreadtograph'])
            print()

        return allAlignedReads


if __name__ == '__main__':
    G = nx.DiGraph()
    G.add_node('Aln_1_1', ids='H37Rv,H37Rv1,H37Rv2', H37Rv_leftend=1, H37Rv_rightend=100, H37Rv1_leftend=1,
               H37Rv1_rightend=100, H37Rv2_leftend=1, H37Rv2_rightend=100,
               sequence='TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATG')
    G.add_node('Aln_1_2', ids='H37Rv,H37Rv1,H37Rv2', H37Rv_leftend=101, H37Rv_rightend=200, H37Rv1_leftend=101,
               H37Rv1_rightend=200, H37Rv2_leftend=-101, H37Rv2_rightend=-200,
               sequence='CTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAG')
    G.add_node('Aln_1_3', ids='H37Rv,H37Rv1,H37Rv2', H37Rv_leftend=201, H37Rv_rightend=100, H37Rv1_leftend=201,
               H37Rv1_rightend=1000, H37Rv2_leftend=201, H37Rv2_rightend=1000,
               sequence='CAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCCCCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGATCCAACTCGGGGTCCGCATCGCTCCGCCGGCGACCGACGAAGCCGACGACACTACCGTGCCGCCTTCCGAAAATCCTGCTACCACATCGCCAGACACCACAACCGACAACGACGAGATTGATGACAGCGCTGCGGCACGGGGCGATAACCAGCACAGTTGGCCAAGTTACTTCACCGAGCGCCCGCACAATACCGATTCCGCTACCGCTGGCGTAACCAGCCTTAACCGTCGCTACACCTTTGATACGTTCGTTATCGGCGCCTCCAACCGGTTCGCGCACGCCGCCGCCTTGGCGATCGCAGAAGCACCCGCCCGCGCTTACAACCCCCTGTTCATCTGGGGCGAGTCCGGTCTCGGCAAGACACACCTGCTACACGCGGCAGGCAACTATGCCCAACGGTTGTTCCCGGGAATGCGGGTCAAATATGTCTCCACCGAGGAATTCACCAACGACTTCATTAACTCGCTCCGCGATGACCGCAAGGTCGCATTCAAACGCAGCTACCGCGACGTAGACGTGCTGTTGGTCGACGACATCCAATTCATTGAAGGCAAAGAGGGTATTCAAGAGGAGTTCTTCCACACCTTCAACACCTTGCACAATGCCAACAAGCAAATCGTCATCTCATCTGACCGCCCACCCAAGCAGCTCGCCACCCTCGAGGACCGGCTGAGAACCCGCTTTGAGTGGGGGCTGATCACTGACGTACAACCACCCG')
    G.add_edge('Aln_1_1', 'Aln_1_2')
    G.add_edge('Aln_1_2', 'Aln_1_3')
    graph_obj = import_gg_graph('./TestGraphs/mix_of_snps.xml')
    alignerGraph = Aligner(G)
    alignerGraph.debruin_read_alignment('CCCCAGTCGCCTCGCGACTCTAATCCAGCTTTGCCAAACGAAATCGAGCG',7)
    #nx.draw(graph_obj, with_labels=True)
    #plt.show()


