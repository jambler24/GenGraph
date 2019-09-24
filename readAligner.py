from gengraph import *
from gengraphTool import *
import networkx as nx
import matplotlib.pyplot as plt



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
        return self.dist


    def debruin_read_alignment(self , query, kmerLength):
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
            connectedToEndNode = list(self[endNode].keys())
            for c in finalAlignedReads:
                distanceBetween = 0
                gapNodes = [endNode]
                continueTest = True
                while continueTest == True:
                    if z != c and list(finalAlignedReads.keys()).index(z) < list(finalAlignedReads.keys()).index(c):
                        if gapNodes[-1] in finalAlignedReads[c]['nodescoveredbyread'][0]:
                            distanceBetween = distanceBetween + (finalAlignedReads[c]['alignmentstartpositioninfirstnode'] - 1 - finalAlignedReads[z]['alignementendpositioninlastnode'])
                            continueTest = False
                        else:
                            self.firstPath = True
                            dist = self.readGaps(list(self[gapNodes[-1]]), query, finalAlignedReads[z], finalAlignedReads[c], distanceBetween)
                            #print(dist)
                            distanceBetween = dist
                            continueTest = False
                    else:
                        continueTest = False
                if distanceBetween < len(query) and distanceBetween != 0:
                    if int(z[13:]) + 1 == int(c[13:]):
                        readGroups.append([z,distanceBetween, c])
        # May have issues with completing the code at this stage. Double check program exits correctly
        print(readGroups)

        #Alignement read groups will shop the start and end point of the reads but will hva egaps in between that can be larger then gaps between the sections in the reads
        finalReadGroups = []
        tempReadGroups = []

        for u in readGroups:
            temp =[u[0], u[-1]]
            tempReadGroups.append(temp)

        print(tempReadGroups)
        print(finalReadGroups)






if __name__ == '__main__':
    graph_obj = import_gg_graph('./test_kmer.xml')
    alignerGraph = Aligner(graph_obj)
    print(graph_obj.nodes['Aln_1_1']['sequence'])
    print(graph_obj.nodes['Aln_1_2']['sequence'])
    print(graph_obj.nodes['Aln_1_4']['sequence'])
    print(graph_obj.nodes['Aln_1_5']['sequence'])
    print(graph_obj.nodes['Aln_1_7']['sequence'])
    print(graph_obj.nodes['Aln_1_8']['sequence'])
    print(graph_obj.nodes['Aln_1_10']['sequence'])
    #print(alignerGraph.fast_kmer_create(4))
    #print(alignerGraph.kmers_of_node('Aln_1_1'))
    #print(alignerGraph.create_query_kmers('TTGACCGATGACCCCGGTT', 3))
    print(alignerGraph.debruin_read_alignment('TTGACCTGACCCGCTTCACCAGTGGAAC',5))
    #nx.draw(graph_obj, with_labels=True)
    #plt.show()


