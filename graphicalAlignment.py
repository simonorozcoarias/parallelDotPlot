import sys
from Bio import SeqIO
from Bio import Align
import matplotlib.pyplot as plt
from Bio import pairwise2
import numpy as np


def graphicalAlignment(seqs1, seqs2, matchScore, mismatchScore, gapOpen, gapExtend):
    sequenceObjects1 = SeqIO.parse(seqs1, "fasta")
    mergedSecuences1 = [str(x.seq) for x in sequenceObjects1]
    horizontalSeq = ''.join(mergedSecuences1)

    sequenceObjects2 = SeqIO.parse(seqs2, "fasta")
    mergedSecuences2 = [str(x.seq) for x in sequenceObjects2]
    verticalSeq = ''.join(mergedSecuences2)

    n = len(horizontalSeq)
    m = len(verticalSeq)
    alphabet = ['A', 'C', 'G', 'T', 'N']
    alpha = len(alphabet)
    windowSize = 24 # el por defecto de dotter es 25

    dotmatrix = [[0 for y in range(m)] for x in range(n)]

    scoreMatrix = [[0 for y in range(n)] for x in range(alpha)]
    newSum = [0 for x in range(n)]
    oldSum = [0 for x in range(n)]
    symVec = verticalSeq

    tmp = [0 for x in range(n)]
    addvec = [0 for x in range(n)]
    delvec = [0 for x in range(n)]

    for i in range(len(horizontalSeq)):
        charac = horizontalSeq[i]
        for j in range(alpha):
            scoreMatrix[j][i] = pairwise2.align.globalms(charac, alphabet[j], matchScore, mismatchScore, gapOpen, gapExtend, score_only=True)

    for i in range(n):
        print("it %d of %d" % (i, n))
        tmp = oldSum
        oldSum = newSum
        newSum = tmp
        index = alphabet.index(str(symVec[i]).upper())
        addvec = scoreMatrix[index]
        if i > windowSize:
            delvec = scoreMatrix[alphabet.index(str(symVec[i-windowSize]).upper())]
        else:
            delvec = [0 for x in range(n)]
        newSum[0] = addvec[0]
        for j in range(1, windowSize):
            newSum[j] = oldSum[j-1]+addvec[j]
        for j in range(windowSize+1, m):
            newSum[j] = oldSum[j-1]+addvec[j]-delvec[j-windowSize]
            if newSum[j] > 0 and i > windowSize:
                dotmatrix[i-int(windowSize/2)][j-int(windowSize/2)] = newSum[j]/windowSize

    averagedMatrix = [[0 for y in range(int(m/windowSize)+1)] for x in range(int(n/windowSize)+1)]

    scoresum = 0
    index = 0
    indey = 0
    i = 0
    j = 0
    while i < len(dotmatrix):
        for x in range(i, i+windowSize):
            for y in range(j, j+windowSize):
                scoresum += dotmatrix[i][j]
        print("len %d, index %d, indey %d" % (len(averagedMatrix), index, indey))
        averagedMatrix[index][indey] = scoresum/(windowSize*windowSize)
        scoresum = 0
        j += windowSize
        indey += 1
        if j > len(dotmatrix[i]):
            j = 0
            indey = 0
            i += windowSize
            index += 1

    imageresult = np.array(averagedMatrix)
    plt.imshow(imageresult, interpolation='nearest', cmap='Greys',
      extent=(0.5,np.shape(imageresult)[0]+0.5,0.5,np.shape(imageresult)[1]+0.5))
    plt.show()
    plt.show()


if __name__ == '__main__':
    # file = sys.argv[1]
    fileseq = "2sequences.fasta"
    matchScore = 1
    mismatchScore = -1
    gapOpen = -5
    gapExtend = -3
    graphicalAlignment(fileseq, fileseq, matchScore, mismatchScore, gapOpen, gapExtend)
