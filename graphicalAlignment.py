import sys
from Bio import SeqIO
from Bio import Align
import matplotlib.pyplot as plt
from Bio import pairwise2
import numpy as np


def graphicalAlignment(seqs1, seqs2, matchScore, mismatchScore):
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
           # scoreMatrix[j][i] = pairwise2.align.globalmx(charac, alphabet[j], matchScore, mismatchScore,score_only=True)
            if charac == alphabet[j]:
               scoreMatrix[j][i] = 5
            else:
                scoreMatrix[j][i] = -4
            # print("%s - %s => %d"%(charac, alphabet[j], scoreMatrix[j][i]))

    for i in range(n):
        #print("it %d of %d" % (i, n))
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
        for j in range(windowSize, m):
            newSum[j] = oldSum[j-1]+addvec[j]-delvec[j-windowSize]
            if newSum[j] > 0 and i > windowSize:
                dotmatrix[i-int(windowSize/2)][j-int(windowSize/2)] = newSum[j]/windowSize

    averagedMatrix = [[0 for y in range(int(m/windowSize)+1)] for x in range(int(n/windowSize)+1)]

    # greymaptool
    minThre = 20
    maxThre = 80
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
        averagedMatrix[index][indey] = (scoresum/(windowSize)) * 256/5*windowSize
        # debemos invertir los valores debido a que entre mas puntaje, mas cercano a negro (cero) deberia ser
        # averagedMatrix[index][indey] = 255 - averagedMatrix[index][indey]
        """if averagedMatrix[index][indey] < minThre:
            averagedMatrix[index][indey] = 0
        elif averagedMatrix[index][indey] > maxThre:
            averagedMatrix[index][indey] = 255"""

        scoresum = 0
        j += windowSize
        indey += 1
        if j > len(dotmatrix[i]):
            j = 0
            indey = 0
            i += windowSize
            index += 1

    # normalizacion de datos (entre 0-255) y aplicacion del filtro (greymaptool)
    maxValue = max(max(averagedMatrix))
    for i in range(len(averagedMatrix)):
        row=""
        for j in range(len(averagedMatrix[i])):
            averagedMatrix[i][j] = int(averagedMatrix[i][j]*255/maxValue)
            if averagedMatrix[i][j] < minThre:
                averagedMatrix[i][j] = 0
            elif averagedMatrix[i][j] > maxThre:
                averagedMatrix[i][j] = 255
            row += str(averagedMatrix[i][j])+" "
        print(row)
    imageresult = np.array(averagedMatrix)
    plt.imshow(imageresult, interpolation='nearest', cmap='Greys',
      extent=(0.5,np.shape(imageresult)[0]+0.5,0.5,np.shape(imageresult)[1]+0.5))
    plt.show()
    plt.show()


if __name__ == '__main__':
    # file = sys.argv[1]
    # fileseq = "2sequences.fasta"
    fileseq = "TE.fasta"
    matchScore = 5 #en el paper de dotter ponen esto
    mismatchScore = -4 #en el paper de dotter ponen esto
    graphicalAlignment(fileseq, fileseq, matchScore, mismatchScore)
