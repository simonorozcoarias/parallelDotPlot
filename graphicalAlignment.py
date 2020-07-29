import sys
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from Bio import pairwise2
import getopt
import time
import multiprocessing


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
    print("      			-t number of threads <int>")
    print("\n")

def seq2array(file, windowSize):
    te = list(SeqIO.parse(file, "fasta"))
    length = int(len(te[0].seq)/windowSize)
    length = int((len(te[0].seq)-length)/windowSize)+1
    seqArray = [str(te[0].seq[windowSize*x+x:windowSize*(x)+x+windowSize]) for x in range(0,length)]
    return(length, seqArray)


def parallelAligmentProt(lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, id):
    init = id*seqs_per_procs + id
    end = (id*seqs_per_procs + id) + seqs_per_procs
    if strand == 0:  # horizontal
        if end > lengthX:
            end = lengthX
        resultLocal = [[pairwise2.align.localmx(seq1Array[i], seq2Array[j], 5, -4, score_only=True, one_alignment_only = True) for j in range(lengthY)] for i in range(init, end)]
    else:  # vertical
        if end > lengthY:
            end = lengthY
        resultLocal = [[pairwise2.align.localmx(seq1Array[i], seq2Array[j], 5, -4, score_only=True) for j in range(init, end)] for i in range(lengthX)]
    print("process %i done" % id)
    return resultLocal

def parallelAligmentDNA(lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, id):
    init = id*seqs_per_procs + id
    end = (id*seqs_per_procs + id) + seqs_per_procs
    if strand == 0:  # horizontal
        if end > lengthX:
            end = lengthX
        resultLocal = [[hammingDistance(seq1Array[i], seq2Array[j], 5, -4) for j in range(lengthY)] for i in range(init, end)]
    else:  # vertical
        if end > lengthY:
            end = lengthY
        resultLocal = [[hammingDistance(seq1Array[i], seq2Array[j], 5, -4) for j in range(init, end)] for i in range(lengthX)]
    print("process %i done" % id)
    return resultLocal

def hammingDistance(seq1, seq2, matScore, misScore):
    if len(seq1) < len(seq2):
       length = len(seq1)
    else:
       length = len(seq2)
    mismatches = [1 for x in range(length) if seq1[x] != seq2[x]]
    numMis = sum(mismatches) + abs(len(seq1) - len(seq2)) 
    score = numMis*misScore + (length-sum(mismatches))*matScore
    if score < 0:
       score = 0
    return score

def graphicalAlignment2(width, height, windowSize, name1, name2, graph_flag, result, lengthX, lengthY, threads):
    resultAveraged = [[(result[x][y]/windowSize) * 256/5*windowSize for y in range(lengthY)] for x in range(lengthX)]
    maxi=max(max(resultAveraged))
    resultAveraged = [[(resultAveraged[x][y]*255) /maxi for y in range(lengthY)] for x in range(lengthX)]
    inicial = np.mean(resultAveraged)

    # GRAY MAP TOOL
    opc = True
    contador = False
    while(opc):
        #umbra, barSize = graytool(resultAveraged, inicial, contador)
        #original_image = umbra
        original_image = np.asarray(resultAveraged)
        resize_image = scale(original_image,width,height)
        plt.figure()
        plt.imshow(resize_image, cmap='Greys')
        plt.title("Dot-Plot alignment -"+str(name1)+" Vs "+str(name2)+"-")
        plt.ylabel(str(name2))
        plt.xlabel(str(name1))
        imagen = "Dot-plot_"+str(name1)+"_Vs_"+str(name2)+"_"+str(window)+"_"+str(threads)+"_"+str(width)+".png"
        plt.savefig(imagen)
        if graph_flag:
            plt.show()
            
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
    return umbra

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
    f.write(">"+name+":"+str(len(sec))+"\n")
    f.write(sec)
    f.close()
    return(archivo, name, inicio, fin, id)

def extract_name(file):
    base =os.path.basename(file)
    nombre=os.path.splitext(base)[0]
    return(nombre)

if __name__ == '__main__':
    start_time = time.time()
    # Declaramos las variables inicales
    window, width, height, graph, threads = 25, 1024, 1024, False, 1
    joined_file1, joined_file2, name1, name2 = '', '', '', ''
    try:
        flags, params = getopt.getopt(sys.argv[1:], "h:i:a:w:d:e:g:t:")
    except getopt.GetoptError as error:
        print(str(error))
        print("Usage: %s -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -d <width> -e <height> -o <bool interactive mode> -t num_threads <int> -help <More info>" % sys.argv[0])
        sys.exit(2)
    for o, a in flags:
        if o == '-h':
            printHelp()
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
            threads = int(a)
 
    if (joined_file1 != '' and joined_file2 != ''):
        lengthX, seq1Array = seq2array(joined_file1, window)
        lengthY, seq2Array = seq2array(joined_file2, window)
        print ("finished creating seq Array")
	# creating the multiprocessing environment
        if lengthX > lengthY:
            n = lengthX
            strand = 0
        else:
            n = lengthY
            strand = 1
        seqs_per_procs = int(n/threads)+1
        end_time = time.time()
        print("module 1 time=", end_time - start_time)

        start_time = time.time()
        # execute sequence alignment in multiprocess mode
        pool = multiprocessing.Pool(processes=threads)
        localresults = [pool.apply_async(parallelAligmentDNA, args=(lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, x)) for x in range(threads)]
        # join all partial results in one
        results = [p.get() for p in localresults]
        result = [[0 for j in range(lengthY)] for i in range(lengthX)]
        x = 0
        y = 0
        if strand == 0:  # horizontal
            for matrix in results:
                for i in range(len(matrix)):
                    for j in range(len(matrix[i])):
                        result[x][y] = matrix[i][j]
                        y += 1
                    x += 1
                    y = 0
                x += 1
                y = 0
        else:  # vertical
            for matrix in results:
                for i in range(len(matrix[0])):
                    for j in range(len(matrix)):
                        result[x][y] = matrix[j][i]
                        x += 1
                    y += 1
                    x = 0
                y += 1
                x = 0
        print ("finished creating matrix")
        end_time = time.time()
        print("module 2 time=", end_time - start_time)

        # create images with created score matrix
        start_time = time.time()
        graphicalAlignment2(width, height, window, name1, name2, graph, result, lengthX, lengthY,threads)
        end_time = time.time()
        print("module 3 time=", end_time - start_time)
    else:
        print("Please insert the sequences to align")
        print("Usage: -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -d <width> -e <height> -g <bool interactive mode> -t num_threads <int> -help <More info>" )
        sys.exit(3)
