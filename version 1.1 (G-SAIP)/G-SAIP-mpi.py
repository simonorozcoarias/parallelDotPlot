import sys
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from Bio import pairwise2
import getopt
import time
import multiprocessing
from mpi4py import MPI

def printHelp():
    print("\nParallel dot-plot allows to align 2 DNA fasta files in multiple CPU's ")
    print("changing various parameters in an interactive way \n")
    print("================================================================================\n")
    print("		Usage: %s" % sys.argv[0])
    print("      			-i <Sequences input file1>")
    print("      			-a <Sequences input file2>")
    print("      			-w WindowSize of paired bases <int> ")
    print("      			-t Number of threads <int>")
    print("      			-m Score calculation mode <string>")
    print("      			-k K-mer value <int>")
    print("      			-d Image width <int>")
    print("      			-e Image height <int>")
    print("      			-g Graphic mode <bool>")
    print("      			-f Image filter <bool>")
    print("      			-N Output image name <string>")
    print("      			-F Output image format <string>")
    print("      			-help <Extended information>")
    print("\n")

def seq2array(file, windowSize):
    te = list(SeqIO.parse(file, "fasta"))
    length = int(len(te[0].seq) / windowSize)
    length = int((len(te[0].seq) - length) / windowSize) #Pilas que aqui se sumaba un 1
    seqArray = [str(te[0].seq[windowSize * x + x:windowSize * x + x + windowSize]) for x in range(0, length)]

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
        resultLocal = [[Distance(seq1Array[i], seq2Array[j], 5, -4, kmer, mode) for j in range(lengthY)] for i in
                       range(init, end)]
    else:  # vertical
        if end > lengthY:
            end = lengthY
        resultLocal = [[Distance(seq1Array[i], seq2Array[j], 5, -4, kmer, mode) for j in range(init, end)] for i in
                       range(lengthX)]
    print("process %i done" % id)
    return resultLocal

def Distance(seq1, seq2, matScore, misScore, k, mode):
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
        if score < 0:
            score = 0
    return score

def graphicalAlignment2(width, height, windowSize, name1, name2, graph_flag, result, lengthX, lengthY, kmer, mode, option, Filter, Name, image_format):

    if mode == 1:
        maxi = windowSize - kmer
    else:
        resultAveraged = [[(result[x][y]/int(windowSize))*256/5*int(windowSize) for y in range(lengthY)] for x in range(lengthX)]
        maxi = max(max(resultAveraged))

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
            # Si el usuario ingresa una bandera que quiere solo el resultad arrojado por el dotter entonces que no se haga esta parte
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

    """
    Joining sequence file
    """
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

    """
    Extracting name from file
    """
    base = os.path.basename(file)
    name = os.path.splitext(base)[0]
    return name

def main():

    """
    Main function
    """
    start_time = time.time()

    ### Parallel process rank assignment
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    rank_msg = '[rank '+str(rank)+' msg]'
    #######################################################################

    # Declaramos las variables inicales
    window, width, height, graph, threads, kmer, mode, size1, size2, option, Filter = None, 1024, 1024, False, 1, 3, -1, 0, 0, None, False
    joined_file1, joined_file2, name1, name2, image_format, Name, file1, file2, inic1, inic2, fin1, fin2, id1, id2 = '', '', '', '', "png", None, None, None, None, None, None, None, None, None
    try:
        flags, params = getopt.getopt(sys.argv[1:], "h:i:a:w:d:e:g:t:k:m:f:N:F:")
    except getopt.GetoptError as error:
        print(str(error))
        print(
            "Usage: %s -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -d <width> -e <height> -o <bool interactive mode> -t num_threads <int> -help <More info>" %
            sys.argv[0])
        sys.exit(2)
    for o, a in flags:
        if o == '-h':
            if rank ==0: printHelp()
            sys.exit()
        elif o == '-i':
            file1 = a
            joined_file1, name1, inic1, fin1, id1, size1 = join_seq(file1)
        elif o == '-a':
            file2 = a
            joined_file2, name2, inic2, fin2, id2, size2 = join_seq(file2)
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
                if rank ==0: print(rank_msg+" Insert True or 1 if you wish run with interactive mode")
                sys.exit(3)
        elif o == '-t':
            threads = int(a)
        elif o == '-m':
            oparallelAligmentDNAption = a.upper()
            if option == 'SLOW':
                mode = 0
            elif option == 'FAST':
                mode = 1
            elif option == 'VERY-FAST':
                mode = 2
            else:
                if rank ==0: print(rank_msg+" Insert a valid option for calculation mode")
                sys.exit(4)
            if rank ==0: print(rank_msg+" Calculation mode: ", option)
        elif o == '-f':
            Filter = a
            if (a == 'True') or (a == '1'):
                Filter = True
            elif (a == 'False') or (a == '0'):
                Filter = False
            else:
                if rank ==0: print(rank_msg+" Insert True or 1 if you wish run with filter mode")
                sys.exit(3)
        elif o == '-j':
            jobSpace = a
        elif o == '-N':
            Name = a
        elif o == '-F':
            image_format = a.lower()
            if image_format != 'png' and image_format != 'pdf' and image_format != 'svg':
                if rank ==0: print(rank_msg+" Insert a valid image format: svg,png or pdf")
                sys.exit(4)
            else:
                image_format = a.lower()
        elif o == '-p':
            splitter = a
        elif o == '-k':
            kmer = int(a)

    if file1 != '' and file2 is None:
        joined_file2, name2, inic2, fin2, id2, size2 = joined_file1, name1, inic1, fin1, id1, size1

    # Here, we define the score calculation mode according to the secuences size. OJO PARA el MANUAL de usuario
    if option is None:
        if size1 >= 1e6 or size2 >= 1e6:  # Fast mode with medium detail
            mode = 2
            kmer = 0
            option = 'VERY-FAST'
            if rank ==0: print(rank_msg+" Default calculation mode: VERY-FAST")
        elif size1 > 50000 or size2 > 50000:
            mode = 1
            option = 'FAST'
            if rank ==0: print(rank_msg+" Default calculation mode: FAST")
        else:  # Slow mode is choosen
            mode = 0
            kmer = 0
            option = 'SLOW'
            if rank ==0: print(rank_msg+" Default calculation mode: SLOW")

    sizes = [size1, size2]
    #print("secuencias len: ", sizes)
    if min(sizes) <= 1000 and window is None:  # Secuencias de 1000 bases
        mode = 0  # No se puede hacer con el fast
        kmer = 0
        window = int(min(sizes)*0.1) + 1 # 1 % de la longitud de la secuencia el +1 es por si el 0.01 % da por debajo de 1
    elif window is None:
        window = int(min(sizes)/512)

    if mode == 1 and window <= kmer: # Condicionar que el kmer no puede ser mayor a la ventana, y decir que use otro método.
        if rank ==0:
            print(rank_msg+" window: ", window, " Kmer: ", kmer)
            print(rank_msg+" Window value smaller tha k-mer value, please, select other Score mode or modify window value")
        sys.exit(5)

    if joined_file1 != '' and joined_file2 != '':
        lengthX, seq1Array = seq2array(joined_file1, window)
        lengthY, seq2Array = seq2array(joined_file2, window)
        if rank ==0: print(rank_msg+" finished creating seq Array")
        # creating the multiprocessing environment
        if lengthX > lengthY:
            n = lengthX
            strand = 0
        else:
            n = lengthY
            strand = 1


        #### MPI parallel reguib

        seqs_per_procs = int(n / threads)
        end_time = time.time()
        if rank ==0: print(rank_msg+" module 1 time=", end_time - start_time)
        start_time = time.time()
        # execute sequence alignment in multiprocess mode
        pool = multiprocessing.Pool(processes=threads)
        localresults = [pool.apply_async(parallelAligmentDNA, args=(
            lengthX, seq1Array, lengthY, seq2Array, seqs_per_procs, strand, x, kmer, mode, n, threads))
                        for x in range(threads)]
        # join all partial results in one
        results = [p.get() for p in localresults]
        end_time = time.time()
        #######################################################################



        #print("Result_len", len(results))
        if rank ==0: print(rank_msg+" module 2 time=", end_time - start_time)
        start_time = time.time()
        if strand == 0:  # horizontal
            result = np.row_stack(results)
        else:  # vertical
            result = np.column_stack(results)
        lengthY = result.shape[1]
        lengthX = result.shape[0]
        if rank ==0: print(rank_msg+" finished creating matrix")
        end_time = time.time()
        if rank ==0: print(rank_msg+" module 3 time=", end_time - start_time)
        # create images with created score matrix
        start_time = time.time()
        graphicalAlignment2(width, height, window, name1, name2, graph, result, lengthX, lengthY, kmer, mode, option, Filter, Name, image_format)
        end_time = time.time()
        if rank ==0: print(rank_msg+" module 4 time=", end_time - start_time)
        #deleting joined_files
        if joined_file1 == joined_file2:
            os.remove(joined_file1)
        else:
            os.remove(joined_file1)
            os.remove(joined_file2)
    else:
        if rank ==0:
            print(rank_msg+" Please insert the sequences to align")
            print(
                "Usage: -i <Sequence input file1> -a <Sequence input file2> -w WindowSize <int> -t num_threads <int> -k kmer_value <int> -m score_calculation_mode <str> -d <width> -e <height> -g interactive mode <bool> -f apply_filter <bool> -N output_name <str> -F output_format <str> -help <More info>")
        sys.exit(3)

if __name__ == '__main__':

    main()
