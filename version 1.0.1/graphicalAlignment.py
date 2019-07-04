
import sys
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from Bio import pairwise2
import getopt

#from matplotlib.lines import Line2D
#from matplotlib.patches import Rectangle

def printHelp():
    print("\nParallel dot-plot allows to align 2 DNA fasta files in multiple CPU's ")
    print("changing various parameters in an interactive way \n")
    print("================================================================================\n")
    print("		Usage: %s" % sys.argv[0])
    print("      			-i <Sequences input file1>")
    print("      			-a <Sequences input file2>")
    print("      			-w WindowSize of paired bases <int> ")
    print("      			-d Image width <int>")
    print("      			-e Image height <int>")
    print("      			-g Graphic mode <bool>")
    print("      			-h <Extended information>")
    print("      			-m <>")
    print("      			-j <>")
    print("      			-p <>")
    print("\n")
    #print("If you require more info about the Parallel Schemes available use \n \n     		")+colored("%s -r" % sys.argv[0],'green')
    #print("\nFor support contact:")+colored("Leonardo Camargo Forero, M.Sc @ leonardo.camargo@bios.co\n",'green')
        ##print("\n")


def seq2array(file, windowSize):
    te = list(SeqIO.parse(file, "fasta"))
    length = int(len(te[0].seq)/windowSize)
    length = int((len(te[0].seq)-length)/windowSize)+1
    seqArray = [str(te[0].seq[windowSize*x+x:windowSize*(x)+x+windowSize]) for x in range(0,length)]
    return(length, seqArray)

def graphicalAlignment2(file1, file2, width, height, windowSize, name1, name2, graph_flag):
    lengthX, seq1Array = seq2array(file1, windowSize)
    lengthY, seq2Array = seq2array(file2, windowSize)
    print(seq1Array)
    print(len(seq1Array))
    print ("finished creating seq Array")
    result = [[pairwise2.align.localmx(seq2Array[y],seq1Array[x],5, -4, score_only=True) for y in range(lengthY)] for x in range(lengthX)]
    print ("finished creating matrix")
    resultAveraged = [[(result[x][y]/windowSize) * 256/5*windowSize for y in range(lengthY)] for x in range(lengthX)]
    maxi=max(max(resultAveraged))
    resultAveraged = [[(resultAveraged[x][y]*255) /maxi for y in range(lengthY)] for x in range(lengthX)]
    inicial = np.mean(resultAveraged)

    #GRAY MAP TOOL
    opc = True
    contador = False
    while(opc):
        umbra, barSize = graytool(resultAveraged, inicial, contador)
        original_image = umbra
        resize_image = scale(original_image,width,height)
        plt.figure()
        plt.imshow(resize_image, cmap='Greys')
        plt.title("Dot-Plot alignment -"+str(name1)+" Vs "+str(name2)+"-")
        plt.ylabel(str(name2))
        plt.xlabel(str(name1))
        imagen = "Dot-plot_"+str(name1)+" Vs "+str(name2)+"_"+str(window)+str(barSize)+str(width)+".png"
        plt.savefig(imagen)
        if graph_flag:
            plt.show()
            #Si el usuario ingresa una bandera que quiere solo el resultad arrojado por el dotter entonces que no se haga esta parte
            opcInt = int(input("Do you want to calculate the matrix again with another value? Yes/No <1/0> :  "))
            if opcInt == 0:
                opc = False
                print("image saved as: "+imagen)
            else:
                contador = True
        else:
            print("image saved as: "+imagen)
            opc = False
    return 0


def perc_graytool(resultAveraged):
    resultAveraged = np.asarray(resultAveraged)
    perc=int(input('Insert de noise elimination percentage:  '))
    print("The inserted percentage is: ", perc)
    thres=(255*perc)/100
    print("The attenuated values are below of : ", thres)
    posi = np.where(np.asarray(resultAveraged) <= thres)
    umbra = np.asarray(resultAveraged)
    umbra[posi] = 0
    return(umbra)

def graytool(resultAveraged, prom, cont):
    resultAveraged = np.asarray(resultAveraged)
    if cont:
        Min_cutoff = int(input('Insert the minimum cutoff:  '))
        Bar_size = int(input('Insert the bar size:  '))
    else:
        Min_cutoff = prom.round(2) +20
        Bar_size = 40
    Max_cutoff = Min_cutoff + Bar_size
    print("Cutoff min , max y Bar Size : ", Min_cutoff,Max_cutoff, Bar_size)
    posi = np.where(resultAveraged <= Min_cutoff)
    umbra = np.asarray(resultAveraged)
    umbra[posi] = 0
    posi1 = np.where(resultAveraged >= Max_cutoff)
    umbra[posi1] = 255
    print("Original image size : ", umbra.shape)
    return(umbra, Bar_size)

def scale(original_image, width, height):
    resize_image = np.zeros(shape=(width, height))
    for W in range(width):
        for H in range(height):
            new_width = int(W * original_image.shape[0] / width)
            new_height = int(H * original_image.shape[1] / height)
            resize_image[W][H] = original_image[new_width][new_height]
    print("Resized image size : ", resize_image.shape)
    return(resize_image)

def join_seq(file):
    name = extract_name(file)
    archivo='joined_'+name+'.fasta'
    f = open(archivo, 'w')
    te = list(SeqIO.parse(file, "fasta"))
    sec, inicio, fin, id, pos_fin = '', [0], [], [], 0
    for tes in te:
        id.append(str(tes.id))
        pos_fin = pos_fin + (len(tes.seq))
        sec = sec + str(tes.seq)
        fin.append(pos_fin-1)
        inicio.append(pos_fin)
    fin.append(0)
    #Cuando se vayan a identificar las secuncias se debe tener en cuenta que no se debe tomar la ultima posicion del
    #vector inicio ya que esta sería el final de las secuencias y por eso se le pone 0 a lo ultimo en el vector fin
    f.write(">"+name+":"+str(len(sec))+"\n")
    f.write(sec)
    f.close()
    return(archivo, name, inicio, fin, id)

def extract_name(file):
    base =os.path.basename(file)
    nombre=os.path.splitext(base)[0]
    return(nombre)

if __name__ == '__main__':
    #Declaramos las variables inicales
    window, width, height, graph = 25, 1024, 1024, False
    joined_file1, joined_file2, name1, name2 = '', '', '', ''
    try:
        flags, params = getopt.getopt(sys.argv[1:], "h:i:a:w:d:e:g:")
    except getopt.GetoptError as error:
        print(str(error))
        print("Usage: %s -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -d <width> -e <height> -o <bool interactive mode> -help <More info>" % sys.argv[0])
        sys.exit(2)
    for o, a in flags:
        if o == '-h':
            printHelp()
            #print("Ayuda")
            sys.exit()
        elif o == '-i':
            file1=a
            joined_file1, name1, inic1, fin1, id1 = join_seq(file1)
        elif o == '-a':
            file2 = a
            joined_file2, name2, inic2, fin2, id2 = join_seq(file2)
        elif o == '-w':
            window = int(a)
        elif o == '-d':
            width = int(a)
        elif o == '-e':
            height = int(a)
        elif o == '-g':
            graph = a
            if (a == 'True') or (a == '1'):
                graph = True
            elif (a == 'False') or (a == '0'):
                graph = False
            else:
                print("Insert True or 1 if you wish run with interactive mode")
                sys.exit(3)
        elif o == '-t':
            typePar = a
        elif o == '-m':
            machinesFile = a
        elif o == '-c':
            print("c")
        elif o == '-j':
            jobSpace = a
        elif o == '-p':
            splitter = a
    if (joined_file1 != '' and joined_file2 != ''):
        graphicalAlignment2(joined_file1, joined_file2,  width, height, window, name1, name2,graph) #800 y 600
    else:
        print("Please insert the sequences to align")
        print("Usage: -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -d <width> -e <height> -o <bool interactive mode> -help <More info>" )
        sys.exit(3)
