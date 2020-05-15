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
    windowSize = 25 # el por defecto de dotter es 25

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
            if charac == alphabet[j]:
               scoreMatrix[j][i] = 5
            else:
                scoreMatrix[j][i] = -4

    for i in range(n):
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
    #print(dotmatrix)
    

    # greymaptool
    minThre = 40
    maxThre = 100
    scoresum = []
    index = 0
    indey = 0
    i = 0
    j = 0
    mem = 0.5 * 8000
    T = 13 #int(np.sqrt((n*m)/mem))
    averagedMatrix = [[0 for y in range(int(m/T)+1)] for x in range(int(n/T)+1)]
    print("T value: ",T)
    while i < len(dotmatrix):
        for x in range(i, i+T):
            for y in range(j, j+T):
                scoresum.append(dotmatrix[i][j])
        #print("len %d, index %d, indey %d" % (len(averagedMatrix), index, indey))
        averagedMatrix[index][indey] = max(scoresum)
        scoresum = []
        j += T
        indey += 1
        if j > len(dotmatrix[i]):
            j = 0
            indey = 0
            i += T
            index += 1
    
    imageresult = np.array(averagedMatrix)
    #imageresult = np.array(dotmatrix)

    plt.figure()

    plt.imshow(imageresult, cmap='Greys')
    plt.savefig("dotter_res.png")


#def winsizeFromLamdak(matriz, tob, alpha)

'''
def karlin():
	##Calculo de lambda

	MAXIT = 20
	SUMLIMIT = 0.01
    up = 0.5
    bandera= 0
    while bandera != 1:
    	up *= 2
    	ptr1 = p
    	suma = 0
    	for i in range(low,high+1,1):
    		ptr1 += 1
    		suma += ptr1 * np.exp(up *i)
    	if suma >= 1.0:
    		bandera = 1

    #Resolviendo por el método de la bisección
    lamda = 0
    for j in range(0,25):
    	new_val = (lamda+up)/2.0;
    	ptr1 = p
    	suma =0
    	for i in range(low,high+1):
    		
    	    suma+= ptr1 * np.exp(new_val*i)
    	    ptr1+=1    	
        if suma > 1.0:
        	up=new_val
        else:
        	lamda = new_val

    beta = np.exp(lamda)

    #Calculando la entropia relativa, parametro H
    ptr1 = p
    av = 0
    go = 0
    for i in range(low,high+1):
    	ptr1 += 1
    	av += ptr1*i*np.exp(lamda*i)
    	H = lamda*av

    if low == -1 or high == 1:
    	if high == 1:
    		K = av
    	else:
    		K= Sum*Sum/av
    	K *= 1.0 -1/beta
    	go = 1

    	#GOTO
    if go == 0:

	    Sum = 0
	    lo, hi = 0, 0
	    P = 0 # Que no conozco
	    P = suma = oldsum = olsum2 = 1
	    j=0
        while j<MAXIT and suma> SUMLIMIT:
	    #for j in range(): # linea que no entiendo
		    first = last = rang
		    hi += high
		    lo += low
		    ptrP = P + hi -lo

		    while ptrP>=P:
		    #for ptrP in range(): #  no se como serían los límites
		        ptr1 = ptrP - first
		        ptr2 = p + first
		        suma = 0
		        for i in range(first,last+1):
		        	ptr1 -= 1 
		        	ptr2 += 1
		        	suma += ptr1 * ptr2
		    	if first != 0:
		    		first -= 1
		        if ptrP - P <= rang:
		        	last = last -1
		        ptrP = suma-1
		    new_val = beta**(lo-1)
		    suma = 0
		    #Aqui va un for pero no sse como hacer esa comparacion
		    i = lo
		    while i!=0:
		    	ptrP += 1
		    	new_val *= beta
		    	suma += ptrP * new_val
		    	i+=1
	    	for j in range(i,hi+1):
	    		ptrP+=1
	    		suma += ptrP
	    	oldsum2 = oldsum 
	    	
    # Progresión geometrica
    go=0
    ratio = olsum/oldsum
    if ratio >= (1- SUMLIMIT*0.001):
    	K=0.1
    	go = 1
    if go ==0:
	    while suma > SUMLIMIT*0.01:
	    	oldsum *= ratio
	    	j+=1
	    	suma = oldsum/j
	    	Suma += suma
	    for i in range(): # este for no lo entiendo
	        for j in range(-i,high): # este otro nada que lo entiendo 
	        	i += 1
	           if p[i-low]:
	           	j = fct_gcd(j,i)

def fct_gcd(a, b):
	b = abs(b)
	if(b>a):
		c=a
		a=b
		b=c
	while(b!=0):
		c=a%b
		a=b
		b=c
	return a
'''


def scale(original_image, width, height):
    resize_image = np.zeros(shape=(width, height))
    for W in range(width):
        for H in range(height):
            new_width = int(W * original_image.shape[0] / width)
            new_height = int(H * original_image.shape[1] / height)
            resize_image[W][H] = original_image[new_width][new_height]
    print("Resized image size : ", resize_image.shape)
    return resize_image

if __name__ == '__main__':
    # file = sys.argv[1]
    # fileseq = "2sequences.fasta"
    fileseq = "TE.fasta.txt"
    matchScore = 5 #en el paper de dotter ponen esto
    mismatchScore = -4 #en el paper de dotter ponen esto
    graphicalAlignment(fileseq, fileseq, matchScore, mismatchScore)