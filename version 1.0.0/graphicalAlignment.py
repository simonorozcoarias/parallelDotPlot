import sys
from Bio import SeqIO
from Bio import Align
#import PyQt5
#import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

from Bio import pairwise2
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
#import opencv-python as cv2
def graphicalAlignment(file2, xSize, ySize):
  f = open(file, 'r')
  tes = SeqIO.parse(file, "fasta")
  result = [[0 for y in range(ySize)] for x in range(xSize)]
  for te in tes:
    aligner = Align.PairwiseAligner()
    lengthX = int(len(te.seq)/xSize)
    lengthY = int(len(te.seq)/ySize)
    startY = 0
    seq1Array = [te.seq[lengthX*x:lengthX*(x+1)] for x in range(xSize)]
    print ("finished creating seq Array")
    result = [[pairwise2.align.localmx(seq1Array[y],seq1Array[x],5, -4, score_only=True) for y in range(ySize)] for x in range(xSize)]
    #result = [[aligner.score(seq1Array[y],seq1Array[x]) for y in range(ySize)] for x in range(xSize)]
    print ("finished creating matrix")
    """for i in range(ySize):
      print "cylce %i of %i"%(i,ySize)
      subseq = te.seq[startY:startY+lengthY]
      startX = 0
      for j in range(xSize):
        subseqX = te.seq[startX:startX+lengthX]
        result[i][j] = aligner.score(subseq,subseqX)
        #result[i][j] = pairwise2.align.globalxx(subseq,subseqX)
        startX += lengthX
      startY += lengthY"""
    maxi = max(max(result))
    #resultImage = [[(result[y][x]*255)/maxi for y in range(len(result))] for x in range(len(result[y]))]
    resultImage = list(np.zeros([len(result),len(result[0])]))
    minTh = 180
    maxTh = 200
    for i in range(len(result)):
      for j in range(len(result[i])):
        resultImage[i][j] = (result[i][j]*255)/maxi
        if resultImage[i][j] < minTh:
          resultImage[i][j] = 0
        elif resultImage[i][j] > maxTh:
          resultImage[i][j] = 255
    print(resultImage)
    """for i in range(len(resultImage)):
      for j in range(len(resultImage[i])):
        if resultImage[i][j] < 180:
          resultImage[i][j] = 0"""
    #print resultImage
    print (te.id)
    plt.imshow(resultImage, interpolation='nearest', cmap='Greys',
      extent=(0.5,np.shape(resultImage)[0]+0.5,0.5,np.shape(resultImage)[1]+0.5))
    plt.show()
    plt.show()

def graphicalAlignment2(file2, xSize, ySize, windowSize):
  f = open(file2, 'r')
  tes = SeqIO.parse(file2, "fasta")
  for te in tes:

    lengthX = int(len(te.seq)/windowSize)
    lengthY = int(len(te.seq)/windowSize)
    lengthX = int((len(te.seq)-lengthX)/windowSize)+1
    lengthY = int((len(te.seq)-lengthY)/windowSize)+1

    startY = 0
    seq1Array = [str(te.seq[windowSize*x+x:windowSize*(x)+x+windowSize]) for x in range(0,lengthX)]

    print(seq1Array)
    print(len(seq1Array))
    print ("finished creating seq Array")
    result = [[pairwise2.align.localmx(seq1Array[y],seq1Array[x],5, -4, score_only=True) for y in range(lengthY)] for x in range(lengthX)]
    print ("finished creating matrix")

    resultAveraged = [[(result[x][y]/windowSize) * 256/5*windowSize for y in range(lengthY)] for x in range(lengthX)]
    maxi=max(max(resultAveraged))

    resultAveraged = [[(resultAveraged[x][y]*255) /maxi for y in range(lengthY)] for x in range(lengthX)]

    #print(resultAveraged)
    #print(te.id)

    Size=np.shape(resultAveraged)
    #print(Size)

    #GRAY MAP TOOL
    opc = True
    while(opc):

      perc=int(input('Ingrese el porcentaje de eliminacion de ruido:  '))

      print ("El porcentaje ingresado es: ",perc)

      thres=(255*perc)/100
      print ("El umbral ingresado es: ",thres)
      posi = np.where(np.asarray(resultAveraged) <= thres)
      umbra = np.asarray(resultAveraged)
      umbra[posi] = 0

      # posi = np.where(np.asarray(resultAveraged) > 255-thres)
      # umbra[posi] = 255


      print("Original image size : ", umbra.shape)

      original_image=umbra
      width , height = 800,600
      resize_image = np.zeros(shape=(width,height))

      for W in range(width):
          for H in range(height):
              new_width = int( W * original_image.shape[0] / width )
              new_height = int( H * original_image.shape[1] / height )
              resize_image[W][H] = original_image[new_width][new_height]

      print("Resized image size : " , resize_image.shape)


      plt.figure()
      plt.imshow(resize_image, cmap='Greys')
      plt.show()

      opcInt = int(input("Desea calcular nuevamente la matriz con otro valor 1/0??:  "))
      if opcInt == 0:
          opc = False
      #print(type(opc))
    return 0


if __name__ == '__main__':
  #file2 = sys.argv[1]
  window = 25 #int(sys.argv[2])
  #res2 = sys.argv[3]
  #file = "coffeeClasifiedTE.fa"
  file2 = "TE.fasta"
  graphicalAlignment2(file2, 800, 600, window)
