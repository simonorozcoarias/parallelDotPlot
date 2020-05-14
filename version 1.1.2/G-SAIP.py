import sys
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from Bio import pairwise2
import time
import multiprocessing
from optparse import OptionParser

def seq2array(file, windowSize):
    te = list(SeqIO.parse(file, "fasta"))
    length = int(len(te[0].seq) / windowSize)
    length = int((len(te[0].seq) - length) / windowSize)
    seqArray = [str(te[0].seq[windowSize * x + x:windowSize * x + x + windowSize]).upper() for x in range(0, length)]

    return length, seqArray

def parallelAligmentDNA(lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, id, kmer, mode, n, threads):
    remain = n % threads
    if id < remain:
        init = id * (seqs_per_procs + 1)
        end = init + seqs_per_procs + 1
    else:
        init = id * seqs_per_procs + remain
        end = init + seqs_per_procs

    if strand == 0:  # horizontal
        if end > lengthX:
            end = lengthX
        resultLocal = [[Distance(seq1Array[i], seq2Array[j], 5, -4, kmer, mode, i, j) for j in range(lengthY)] for i in
                       range(init, end)]
    else:  # vertical
        if end > lengthY:
            end = lengthY
        resultLocal = [[Distance(seq1Array[i], seq2Array[j], 5, -4, kmer, mode, i, j) for j in range(init, end)] for i in
                       range(lengthX)]
    print("process %i done" % id)
    return resultLocal

def Distance(seq1, seq2, matScore, misScore, k, mode, x, y):
    if mode == 0:  # SLOW (Pairwise2 Biopython) -- HIGH detail
        score = pairwise2.align.globalmx(seq1, seq2, matScore, misScore, score_only=1)
    elif mode == 1:  # FAST (subsequences with kmer) -- MEDIUM detail
        score = 0
        for i in range(len(seq1) - k):
            subseq = seq1[i:i + k]
            if subseq in seq2:
                score += 1
    elif mode == 2:  # VERY - FAST (Hamming Distance) -- LOW detail
        if len(seq1) < len(seq2):
            length = len(seq1)
        else:
            length = len(seq2)
        mismatches = [1 for x in range(length) if seq1[x] != seq2[x] or seq1[x] == 'N' or seq2[x] == 'N']
        numMis = sum(mismatches) + abs(len(seq1) - len(seq2))
        score = numMis * misScore + (length - sum(mismatches)) * matScore
        '''if score < 0:
            score = 0'''
        print("Score for slice ",x,",",y,": ",score)
    return score


def graphicalAlignment2(width, height, windowSize, name1, name2, graph_flag, result, lengthX, lengthY, threads, kmer, mode, option, Filter, Name, image_format):

    if mode == 1:
        maxi = windowSize - kmer
    else:
        print("La matriz queda de score : ")
        print(result)
        minimo = result.min()
        maximo = result.max()
        resultAveraged = [[(((result[x][y]-minimo)*(255))/(maximo-minimo)) for y in range(lengthY)]for x in range(lengthX)]
        #resultAveraged = [[(result[x][y]/int(windowSize))*256/5*int(windowSize) for y in range(lengthY)] for x in range(lengthX)]
        maxi = max(max(resultAveraged))
        print("Max num:", maxi, "Orig_max:",maximo,"Orig_min:",minimo)

    resultAveraged = [[(result[x][y] * 255) / maxi for y in range(lengthY)] for x in range(lengthX)]
    inicial = np.mean(resultAveraged) + np.std(resultAveraged)

    # GRAY MAP TOOL
    opc = True
    contador = False
    while (opc):
        if Filter:
            umbra, barSize = graytool(resultAveraged, inicial, contador)
            original_image = umbra
        else:
            original_image = np.asarray(resultAveraged)
        resize_image = scale(original_image, width, height)
        plt.figure()
        plt.imshow(resize_image, cmap='Greys')
        plt.title("G-SAIP -" + str(name1) + " Vs " + str(name2) + "-")
        plt.ylabel(str(name2))
        plt.xlabel(str(name1))
        ##PonerParametro
        if Name is None:
            date = time.localtime( time.time() )
            imagen = "G-SAIP" + "%04d"%date.tm_year + "%02d"%date.tm_mon + "%02d"%date.tm_mday + "%02d"%date.tm_hour + "%02d"%date.tm_min + "%02d"%date.tm_sec + "."+ str(image_format)
        else:
            imagen = str(Name)+"."+str(image_format)
        plt.savefig(imagen)
        if graph_flag:
            plt.show()
            
            opcInt = int(input("Do you want to calculate the matrix again with another value? Yes/No <1/0> :  "))
            if opcInt == 0:
                opc = False
                print("image saved as: " + imagen)
            else:
                contador = True
        else:
            print("image saved as: " + imagen)
            opc = False
    return 0

def graytool(resultAveraged, prom, cont):
    resultAveraged = np.asarray(resultAveraged)
    if cont:
        Min_cutoff = int(input('Insert the minimum cutoff:  '))
        Bar_size = int(input('Insert the bar size:  '))
    else:
        Min_cutoff = prom.round(2) + 20
        Bar_size = 40
    Max_cutoff = Min_cutoff + Bar_size
    print("Cutoff min , max y Bar Size : ", Min_cutoff, Max_cutoff, Bar_size)
    posi = np.where(resultAveraged <= Min_cutoff)
    umbra = np.asarray(resultAveraged)
    umbra[posi] = 0
    posi1 = np.where(resultAveraged >= Max_cutoff)
    umbra[posi1] = 255
    print("Original image size : ", umbra.shape)
    return umbra, Bar_size

def scale(original_image, width, height):
    resize_image = np.zeros(shape=(width, height))
    for W in range(width):
        for H in range(height):
            new_width = int(W * original_image.shape[0] / width)
            new_height = int(H * original_image.shape[1] / height)
            resize_image[W][H] = original_image[new_width][new_height]
    print("Resized image size : ", resize_image.shape)
    return resize_image

def join_seq(file):
    name = extract_name(file)
    archivo = 'joined_' + name + '.fasta'
    f = open(archivo, 'w')
    te = list(SeqIO.parse(file, "fasta"))
    sec, inicio, fin, id, pos_fin = '', [0], [], [], 0
    for tes in te:
        id.append(str(tes.id))
        pos_fin = pos_fin + (len(tes.seq))
        sec = sec + str(tes.seq)
        fin.append(pos_fin - 1)
        inicio.append(pos_fin)
    fin.append(0)
    size = len(sec)
    # Cuando se vayan a identificar las secuncias se debe tener en cuenta que no se debe tomar la ultima posicion del
    # vector inicio ya que esta seria el final de las secuencias y por eso se le pone 0 a lo ultimo en el vector fin
    f.write(">" + name + ":" + str(len(sec)) + "\n")
    f.write(sec)
    f.close()
    return archivo, name, inicio, fin, id, size

def extract_name(file):
    base = os.path.basename(file)
    nombre = os.path.splitext(base)[0]
    return nombre

if __name__ == '__main__':
    start_time = time.time()
    # Declaramos las variables inicales
    usage = "usage: python G-SAIP.py -q file.fasta ... [options]"
    parser = OptionParser(usage=usage)

    parser.add_option('-q','--query',dest='file1',type=str,default=None,help='Query sequence in FASTA format')
    parser.add_option('-s','--subject',dest='file2',type=str,default=None,help='Subject sequence in FASTA format')
    parser.add_option('-w','--window',dest='window',type=int,default=None,help='Window size')
    parser.add_option('-t','--threads',dest='threads',type=int,default=1,help='Number of threads to use')
    parser.add_option('-k','--kmer',dest='kmer',type=int,default=3,help='kmer value')
    parser.add_option('-m','--mode',dest='option',type=str,default=None,help='Calculation score mode')
    parser.add_option('-i','--width',dest='width',type=int,default=1024,help='Output image width')
    parser.add_option('-e','--height',dest='height',type=int,default=1024,help='Output image height')
    parser.add_option('-g','--graphic',dest='graph',type=str,default='FALSE',help='Interactive mode')
    parser.add_option('-f','--filter',dest='Filter',type=str,default='FALSE',help='Apply default filter')
    parser.add_option('-n','--name',dest='Name',type=str,default=None,help='Output image name')
    parser.add_option('-o','--format',dest='image_format',type=str,default='png',help='Specify the output image format')


    (options,arguments) = parser.parse_args()

    file1=options.file1
    file2=options.file2
    window=options.window
    threads=options.threads
    kmer=options.kmer
    option=options.option
    width=options.width
    height=options.height
    graph=options.graph
    Filter=options.Filter
    Name=options.Name
    image_format=options.image_format


    size1, size2= 0, 0
    joined_file1, joined_file2, name1, name2, inic1, inic2, fin1, fin2, id1, id2 = '', '', '', '', None, None, None, None, None, None

    if file1 is not None:
    	joined_file1, name1, inic1, fin1, id1, size1 = join_seq(file1)
    else:
    	print("Please insert at least a Query file in FASTA format")
    	sys.exit(1)

    if file2 is not None:
    	joined_file2, name2, inic2, fin2, id2, size2 = join_seq(file2)

    if file1 != '' and file2 is None:
        joined_file2, name2, inic2, fin2, id2, size2 = joined_file1, name1, inic1, fin1, id1, size1

    if option is not None:
        if option.upper() == 'SLOW':
            mode = 0
        elif option.upper() == 'FAST':
            mode = 1
        elif option.upper() == 'VERY-FAST':
            mode = 2
        else:
            print("Insert a valid option for calculation mode")
            sys.exit(4)
        print("Calculation mode: ", option)

    if image_format.lower() != 'png' and image_format.lower() != 'pdf' and image_format.lower() != 'svg':
        print("Insert a valid image format: svg,png or pdf")
        sys.exit(4)
    else:
        image_format = image_format.lower()

    if (graph.upper() == 'TRUE') or (str(graph) == '1'):
        graph = True
    elif (graph.upper() == 'FALSE') or (str(graph) == '0'):
        graph = False
    else:
        print("Insert True or 1 if you wish run with interactive mode")
        sys.exit(3)
    
    if (Filter.upper() == 'TRUE') or (str(Filter) == '1'):
        Filter = True
    elif (Filter.upper() == 'FALSE') or (str(Filter) == '0'):
        Filter = False
    else:
        print("Insert True or 1 if you wish run with filter mode")
        sys.exit(3)
    
    # Here, we define the score calculation mode according to the secuences size. OJO For user manual
    if option is None:
        if size1 >= 1e6 or size2 >= 1e6:  # Fast mode with medium detail
            mode = 2
            kmer = 0
            option = 'VERY-FAST'
            print("Default calculation mode: VERY-FAST")
        elif size1 > 50000 or size2 > 50000:
            mode = 1
            option = 'FAST'
            print("Default calculation mode: FAST")
        else:  # Slow mode is choosen
            mode = 0
            kmer = 0
            option = 'SLOW'
            print("Default calculation mode: SLOW")
    
    sizes = [size1, size2]
    print("secuencias len: ", sizes)
    if min(sizes) <= 1000 and window is None:  # Secuencias de 1000 bases
        mode = 0  # No se puede hacer con el fast
        kmer = 0
        window = int(min(sizes)*0.1) + 1 # 1 % de la longitud de la secuencia el +1 es por si el 0.01 % da por debajo de 1
    elif window is None:
        window = int(min(sizes)/512)
    if mode == 1 and window <= kmer: # El kmer no puede ser mayor a la ventana, se debe usar otro método.
        print("window: ", window, " Kmer: ", kmer)
        print("Window value smaller tha k-mer value, please, select other Score mode or modify window value")
        sys.exit(5)


    lengthX, seq1Array = seq2array(joined_file1, window)
    lengthY, seq2Array = seq2array(joined_file2, window)
    print("finished creating seq Array")
    # creating the multiprocessing environment
    if lengthX > lengthY:
        n = lengthX
        strand = 0
    else:
        n = lengthY
        strand = 1
    print("lenX: ", lengthX, "lenY: ", lengthY, "n: ", n, "strand: ", strand, "Window", window)
    seqs_per_procs = int(n / threads)
    print("seqs_per_procs", seqs_per_procs)
    end_time = time.time()
    print("module 1 time=", end_time - start_time)
    start_time = time.time()
    # execute sequence alignment in multiprocess mode
    pool = multiprocessing.Pool(processes=threads)
    localresults = [pool.apply_async(parallelAligmentDNA, args=(
        lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, x, kmer, mode, n, threads))
                    for x in range(threads)]
    # join all partial results in one
    results = [p.get() for p in localresults]
    end_time = time.time()
    #print("Result_len", len(results))
    print("module 2 time=", end_time - start_time)
    start_time = time.time()
    if strand == 0:  # horizontal
        result = np.row_stack(results)
    else:  # vertical
        result = np.column_stack(results)
    lengthY = result.shape[1]
    lengthX = result.shape[0]
    print("finished creating matrix")
    end_time = time.time()
    print("module 3 time=", end_time - start_time)
    # create images with created score matrix
    start_time = time.time()
    graphicalAlignment2(width, height, window, name1, name2, graph, result, lengthX, lengthY, threads, kmer, mode, option, Filter, Name, image_format)
    end_time = time.time()
    print("module 4 time=", end_time - start_time)
    #deleting joined_files
    if joined_file1 == joined_file2:
        os.remove(joined_file1)
    else:
        os.remove(joined_file1)
        os.remove(joined_file2)