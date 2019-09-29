from gengraph import *
from gengraphTool import *
import networkx as nx
import matplotlib.pyplot as plt
import difflib


class Aligner(nx.DiGraph):

    kmerDict = {}
    global firstPath
    firstPath = True

    def get_number_isolates(self):
        x = []

        return self.nodes


    def fast_kmer_create(self, kmerLength):
        # Start Positions stored in the dictionary are character list positions, so the first character in the node sequence is position 0 not position 1
        graphNodes = list(self.nodes)
        kmerNumber = 1


        for i in graphNodes:
            currentNodeName = str(i)
            nodeSequence = self.nodes[currentNodeName]['sequence']
            kmerStartPos = 0
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
                    self.kmer_over_multiple_nodes(kmerPartseq, currentNode, kmerLength, kmerNumber,dictStartPos)
                    kmerNumber += 1
        return Aligner.kmerDict

    def kmer_over_multiple_nodes(self, kmerPart, currentNode, kmerLength, kmerNumber,kmerStartPos):
        remainingNumber = kmerLength - len(kmerPart)
        connectedNodes = list(self[currentNode[-1]])
        connectedNodeNumber = len(connectedNodes)

        for k in range(connectedNodeNumber):
            nextNode = self.nodes[connectedNodes[k]]
            nextNodeSeq = nextNode['sequence']
            if len(nextNodeSeq) >= remainingNumber:
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
                self.kmer_over_multiple_nodes(kmerTemp, nodesInvolved,kmerLength,kmerNumber,kmerStartPos)
                # recursively goes over multiple nodes
                # will need the currentNode that is parsed at this step so instead you have the two previous nodes the kmer goes over so you have all three nodes saved into the dictionary

    def kmers_of_node(self, nodeName):
        kmersOfNode = {}
        if len(Aligner.kmerDict) == 0:
            #default kmer length is 4
            Aligner.fast_kmer_create(self, 4)

        keyList = list(Aligner.kmerDict.keys())

        for i in keyList:
            currentKmer = Aligner.kmerDict.get(i)
            kmerNodeList = currentKmer['nodes']
            for j in kmerNodeList:
                if j == nodeName:
                    kmer = {i: currentKmer}
                    kmersOfNode.update(kmer)
        return kmersOfNode

    def create_query_kmers(self, querySequence, kmerLength):
        # returns a dictionary in order of the kmers of the query
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
        #print(distanceBetween)
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
        Running this Function will also print out blocks of information on the aligned read. Each block of 6 lines will correspond to a single aligned read.
        :param query: This is the read, as a string, that you want to try and align to the graph.
        :param kmerLength:
        :return: returns a dictionary with the sequences that align, x's represent sections of the sequence that doesnt correctly align
                It also returns the nodes that the sequence aligns to and the position on the first node the sequence starts aligning to and the position on the last node that the sequence ends aligning to.
                These positions are in string positions, so the first base pair is position 0.
        '''
        queryKmerDict = self.create_query_kmers(query, kmerLength)
        referenceKmerDict = self.fast_kmer_create(kmerLength)
        finalkmerGroups = []
        groupedListOftwos = []
        matchedKmers = []
        finalAlignedReads = {}
        #list of lists of matched kmers
        #first pos is query kmer name, second pos is reference nodes that the kmer covers, third is the start position of reference kmer, fourth is reference kmer name, 5th is qkmerseq

        queryKmerKeys = list(queryKmerDict.keys())
        referenceKmerKeys = list(referenceKmerDict.keys())


        for i in queryKmerKeys:
            currentQueryKmer = queryKmerDict.get(i)
            for j in referenceKmerKeys:
                currentReferenceKmer = referenceKmerDict.get(j)
                if currentReferenceKmer['sequence'] == currentQueryKmer:
                    tempList = [i, currentReferenceKmer['nodes'], currentReferenceKmer['startposition'], j, currentQueryKmer]
                    matchedKmers.append(tempList)


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
                        tempGroupedList.append(addOne)
                        tempGroupedList.append(addTwo)
                        groupedListOftwos.append(tempGroupedList)
                    elif len(self.nodes[k[1][0]]['sequence']) - 1 == k[2]:
                        if len(k[1]) > 1:
                            if k[1][1] == x[1][0] and x[2] == 0:
                                addOne = str(k[0]) + '-' + str(k[3])
                                addTwo = str(x[0]) + '-' + str(x[3])
                                tempGroupedList.append(addOne)
                                tempGroupedList.append(addTwo)
                                groupedListOftwos.append(tempGroupedList)

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
        AlignNumber = 1
        #print(finalkmerGroups)
        for n in finalkmerGroups:
            matchedSquence = ''
            nodesCovered = []
            startpos = 0
            endpos = 0
            finalKmerName = n[-1][n[-1].find('-') + 1:]
            finalNode = referenceKmerDict[finalKmerName]['nodes'][-1]
            for m in n:
                kmerName = m[m.find('-') + 1:]
                kmerseq = referenceKmerDict[kmerName]['sequence']
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
        print(finalAlignedReads)
        # finalAlignedReads Are all reads that map 100 percent accurately
        # Start grouping the final aligned read dictionary to show bigger reads and not only 100% alignment
        finalGroupedReads = {}
        readGroups = []

        for z in finalAlignedReads:
            endNode = finalAlignedReads[z]['nodescoveredbyread'][-1]
            for c in finalAlignedReads:
                distanceBetween = 0
                startNode = finalAlignedReads[c]['nodescoveredbyread'][0]
                if z!=c and list(finalAlignedReads.keys()).index(z) < list(finalAlignedReads.keys()).index(c):
                    if int(z[13:]) + 1 == int(c[13:]):
                        if endNode == startNode:
                            distanceBetween = finalAlignedReads[c]['alignmentstartpositioninfirstnode'] - finalAlignedReads[z]['alignementendpositioninlastnode']
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
        #print(readGroups)


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

        print(finalReadGroups)

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
                    for v in range(int(z[i])):
                        xCounter += 1
                        sequence = sequence + 'x'
            nodeNos = []
            newNodesCoveredByRead = []
            for r in nodesCoveredByRead:
                nodeNos.append(r[str(r).rfind('_') + 1:])
            for a in nodeNos:
                for b in nodesCoveredByRead:
                    if a == b[str(b).rfind('_') + 1:]:
                        newNodesCoveredByRead.append(b)

            percentOfQuery = difflib.SequenceMatcher(None, query, checkPercentSeq).ratio() * 100
            percentAligned = (len(sequence)-xCounter)/len(sequence)*100
            tempVal = {'sequence':sequence, 'nodescoveredbyread':newNodesCoveredByRead, 'alignmentstartpositioninfirstnode':startAlignPos, 'alignementendpositioninlastnode':endAlignPos,
                       'percentageofqueryaligned':percentOfQuery,'percentageofalignedreadtograph':percentAligned}
            key = {'Alignment_'+str(AlignNo):tempVal}
            AlignNo += 1
            allAlignedReads.update(key)
        print(allAlignedReads)
        print('6 line block format for each aligned read')
        print('Line 1: The Aligned read sequence ')
        print('Line 2: The Nodes that the aligned read covers')
        print('Line 3: The starting alignment position on the first node involved in the alignment')
        print('Line 4: The ending alignment position on the last node involved in the alignment')
        print('Line 5: Percentage of the initial query that has aligned to the graph')
        print('Line 6:Percentage of the aligned read that has successfully aligned to the graph')
        print()

        for t in list(allAlignedReads.keys()):
            print(allAlignedReads[t]['sequence'])
            print(allAlignedReads[t]['nodescoveredbyread'])
            print(allAlignedReads[t]['alignmentstartpositioninfirstnode'])
            print(allAlignedReads[t]['alignementendpositioninlastnode'])
            print(allAlignedReads[t]['percentageofqueryaligned'])
            print(allAlignedReads[t]['percentageofalignedreadtograph'])
            print()

        return allAlignedReads


if __name__ == '__main__':
    graph_obj = import_gg_graph('./TestGraphs/test_kmer.xml')
    alignerGraph = Aligner(graph_obj)
    #print(graph_obj.nodes['Aln_1_1']['sequence'])
    #print(graph_obj.nodes['Aln_1_2']['sequence'])
    #print(graph_obj.nodes['Aln_1_4']['sequence'])
    #print(graph_obj.nodes['Aln_1_5']['sequence'])
    #print(graph_obj.nodes['Aln_1_7']['sequence'])
    #print(graph_obj.nodes['Aln_1_8']['sequence'])
    #print(graph_obj.nodes['Aln_1_10']['sequence'])
    #print(alignerGraph.fast_kmer_create(4))
    #print(alignerGraph.kmers_of_node('Aln_1_1'))
    #print(alignerGraph.create_query_kmers('TTGACCGATGACCCCGGTT', 3))
    alignerGraph.debruin_read_alignment('TTGACCTGACCCGCTTCACCAGTGGAAC',4)
    #nx.draw(graph_obj, with_labels=True)
    #plt.show()


