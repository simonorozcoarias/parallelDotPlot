import sys
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os,glob 
from os import path
from Bio import pairwise2
import time
from optparse import OptionParser
from mpi4py import MPI
import pickle
import subprocess


def seq2array(file, windowSize):
    te = list(SeqIO.parse(file, "fasta"))
    length = int(len(te[0].seq) / windowSize)
    length = int((len(te[0].seq) - length) / windowSize)
    seqArray = [str(te[0].seq[windowSize * x + x:windowSize * x + x + windowSize]).upper() for x in range(0, length)]

    return length, seqArray


def dna_parallel_alignment(namesX,namesY,strand,threads,sim,slen,kmer):
    """
    Performing DNA alignment in parallel
    """

    if strand==0:
        #
        
        local_result = Distance('vertical', 'horizontal', namesY, namesX, strand, sim, slen, kmer)
        #local_result = np.column_stack(local_result)
        #print(local_result)

    else: 
        local_result = Distance('horizontal', 'vertical', namesX, namesY, strand, sim, slen, kmer)
        #local_result=np.row_stack(local_result)
        #print(local_result)
        #
    return local_result


def Distance(ref,subj,sizes_ref,sizes_subj,strand,sim,slen,kmer):
    
    ''' This function run MashMap and extract the alignment score   '''
    #print("strand:",strand)
    #print(sizes_ref,sizes_subj)
    #for i in range(1,len(sizes_subj)+1): #este for itera el archivo del .rank sobre la secuencia más larga.
    #print(i,rank-1)
    command= "mashmap -r "+ref+"."+str(rank)+" -q "+str(subj)+".unique -t 1 --pi "+str(sim)+ " -s "+str(slen)+" -k "+str(kmer)+" -f none -o mashmap.temp."+str(rank)
    start_time=time.time()
    subprocess.run(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #subprocess.run(command,shell=True)
    # me falta definir el tamaño de la matriz que voy devolver acá.
    if strand==0:
        size = (sizes_ref[rank-1],sizes_subj)
    else:
        size = (sizes_subj,sizes_ref[rank-1])
    matrix = np.zeros(size)
    end_time=time.time()
    #acá debo sacar el score de los pedazos y devolver la matriz que se ha alineando
    #start_time=time.time()
    temp_name="mashmap.temp."+str(rank)
    exist = path.exists(temp_name)
    if exist:
        file = open(temp_name,'r')
        for line in file:
            score = float(line.split()[-1])
            row = line.split()
            if row[4]=='+':
                posj= int(row[0].split('_')[-1])
                posi= int(row[5].split('_')[-1])
                if strand ==0:
                    #print(posi,posj)
                    matrix[posi][posj] = (score*255)/100
                else:
                    #print(posj,posi)
                    matrix[posj][posi] = (score*255)/100
        file.close()
        os.remove(temp_name)
    
    
    #plt.imshow(matrix)
    #plt.savefig("imagen"+str(rank)+"_Proc:"+str(i))
    return matrix


def graphicalAlignment2(width, height, windowSize, name1, name2, graph_flag, result, lengthX, lengthY, kmer, mode, option, Filter, Name, image_format):


    resultAveraged = result            
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
        #plt.title("G-SAIP -" + str(name1) + " Vs " + str(name2) + "-")
        plt.ylabel(str(name2))
        plt.xlabel(str(name1))
        ##PonerParametro
        if Name is None:
            date = time.localtime(time.time())
            imagen = "G-SAIP" + "%04d" % date.tm_year + "%02d" % date.tm_mon + "%02d" % date.tm_mday + "%02d" % date.tm_hour + "%02d" % date.tm_min + "%02d" % date.tm_sec + "." + str(
                image_format)
        else:
            imagen = str(Name) + "." + str(image_format)
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
    Each process creates its own files
    """
    name = extract_name(file)
    #fasta_file = 'joined_' + name + '_' + str(rank) + '.fasta'
    fasta_file = 'joined_' + name + '.fasta'
    f = open(fasta_file, 'w')
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
    return fasta_file, name, inicio, fin, id, size

def extract_name(file):
    """
    Extracting name from file
    """
    base = os.path.basename(file)
    name = os.path.splitext(base)[0]
    return name

def receive_mpi_msg(src=MPI.ANY_SOURCE, t=MPI.ANY_TAG, deserialize=False):
    """
    Sending MPI messages
    """
    data_dict = {}
    status = MPI.Status()
    data = comm.recv(source=src, tag=t, status=status)
    if deserialize: data = pickle.loads(data)
    data_dict['data'] = data
    data_dict['sender'] = int(status.Get_source())
    return data_dict


def send_mpi_msg(destination, data, serialize=False):
    """
    Sending MPI messages
    """
    if serialize: data = pickle.dumps(data)
    comm.send(data, dest=destination)

def CreateFiles(seqs_per_procs,ArraySplit,ArrayUnique,lengthSplit,lengthUnique,typeSplit,typeUnique,threads,n):
    # ESe divide en Rank pedazos el vertical
    procs = threads -1
    init = 0
    end = 0
    namesS = []
    remain = n % (threads-1)
    for i in range(0,threads-1):  
        if i < remain:
            init = i* (seqs_per_procs + 1)
            end = init + seqs_per_procs + 1
        else:
            init = i * seqs_per_procs + remain
            end = init + seqs_per_procs
        if end > lengthSplit:
            end = lengthSplit
        fasta_name = str(typeSplit)+"."+str(i+1)
        file = open(fasta_name,"w")
        num_slices=0
        for j in range(init,end):
            file.write(">pos_"+str(typeSplit)+"_"+str(j)+"_"+str(num_slices)+"\n"+ArraySplit[j]+"\n")
            num_slices+=1
        file.close()
        namesS.append(num_slices)
    # Ahora solo 1 horizontal donde están todas las secuencias
    namesU = 0
    num_slices = 0
    fasta_name = str(typeUnique)+".unique"
    file = open(fasta_name,"w")
    for j,seq in enumerate(ArrayUnique):
        file.write(">pos"+str(typeUnique)+"_"+str(j)+"_"+str(j)+"\n"+seq+"\n")
        num_slices+=1
    file.close()
    namesU=num_slices
    return namesS,namesU

def main():
    """
    Main function
    """
    start_time = time.time()

    ### Parallel process rank assignment
    global comm, threads, rank, rank_msg
    global bandera
    comm = MPI.COMM_WORLD
    threads = comm.Get_size()
    rank = comm.Get_rank()
    rank_msg = '[rank ' + str(rank) + ' msg]'
    #######################################################################
    #print(threads)
    # Declaramos las variables inicales
    usage = "usage: mpirun -np <threads> python3 G-SAIP.py -q file.fasta ... [options]"
    parser = OptionParser(usage=usage)

    parser.add_option('-q', '--query', dest='file1', type=str, default=None, help='Query sequence in FASTA format')
    parser.add_option('-s', '--subject', dest='file2', type=str, default=None, help='Subject sequence in FASTA format')
    parser.add_option('-w', '--window', dest='window', type=int, default=None, help='Window size')
    parser.add_option('-k', '--kmer', dest='kmer', type=int, default=16, help='Mashmap kmer value')
    parser.add_option('-m', '--mode', dest='option', type=str, default=None, help='Calculation score mode')
    parser.add_option('-i', '--width', dest='width', type=int, default=1024, help='Output image width')
    parser.add_option('-e', '--height', dest='height', type=int, default=1024, help='Output image height')
    parser.add_option('-g', '--graphic', dest='graph', type=str, default='FALSE', help='Interactive mode')
    parser.add_option('-f', '--filter', dest='Filter', type=str, default='FALSE', help='Apply default filter')
    parser.add_option('-n', '--name', dest='Name', type=str, default=None, help='Output image name')
    parser.add_option('-o', '--format', dest='image_format', type=str, default='png', help='Specify the output image format')
    parser.add_option('-p', '--pidentity', dest='sim', type=int, default=95, help='Mashmap similarity percentage')
    parser.add_option('-S', '--segmentlen', dest='slen', type=int, default=500, help='Mashmap segment length')
    
    (options, arguments) = parser.parse_args()

    file1 = options.file1
    file2 = options.file2
    window = options.window
    kmer = options.kmer
    option = options.option
    width = options.width
    height = options.height
    graph = options.graph
    Filter = options.Filter
    Name = options.Name
    image_format = options.image_format
    sim = options.sim
    slen = options.slen

    size1, size2 = 0, 0
    joined_file1, joined_file2, name1, name2, inic1, inic2, fin1, fin2, id1, id2 = '', '', '', '', None, None, None, None, None, None

    if file1 is not None:
        if rank == 0:  
            joined_file1, name1, inic1, fin1, id1, size1 = join_seq(file1)
            for m in range(1,threads):
                send_mpi_msg(m,[joined_file1, name1, inic1, fin1, id1, size1],serialize=True)
        else:
            a = receive_mpi_msg(deserialize=True)
            [joined_file1, name1, inic1, fin1, id1, size1] = a['data']
        
    else:
        if rank == 0:
            print("Please insert at least a Query file in FASTA format")
        sys.exit(1)

    if file2 is not None:
        if rank == 0:  
            joined_file2, name2, inic2, fin2, id2, size2 = join_seq(file2)
            for m in range(1,threads):
                send_mpi_msg(m,[joined_file2, name2, inic2, fin2, id2, size2],serialize=True)
        else:
            a = receive_mpi_msg(deserialize=True)
            [joined_file2, name2, inic2, fin2, id2, size2] = a['data']

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
            if rank == 0:
                print("Insert a valid option for calculation mode")
            sys.exit(4)
        print("Calculation mode: ", option)

    if image_format.lower() != 'png' and image_format.lower() != 'pdf' and image_format.lower() != 'svg':
        if rank == 0:
            print("Insert a valid image format: svg,png or pdf")
        sys.exit(4)
    else:
        image_format = image_format.lower()

    if (graph.upper() == 'TRUE') or (str(graph) == '1'):
        graph = True
    elif (graph.upper() == 'FALSE') or (str(graph) == '0'):
        graph = False
    else:
        if rank == 0:
            print("Insert True or 1 if you wish run with interactive mode")
        sys.exit(3)

    if (Filter.upper() == 'TRUE') or (str(Filter) == '1'):
        Filter = True
    elif (Filter.upper() == 'FALSE') or (str(Filter) == '0'):
        Filter = False
    else:
        if rank == 0:
            print("Insert True or 1 if you wish run with filter mode")
        sys.exit(3)

    # Here, we define the score calculation mode according to the secuences size. OJO PARA el MANUAL de usuario
    if option is None:
        if size1 >= 1e6 or size2 >= 1e6:  # Fast mode with medium detail
            mode = 2
            #kmer = 0
            option = 'VERY-FAST'
            if rank == 0: print(rank_msg + " Default calculation mode: VERY-FAST")
        elif size1 > 50000 or size2 > 50000:
            mode = 1
            kmer = 5
            option = 'FAST'
            if rank == 0: print(rank_msg + " Default calculation mode: FAST")
        else:  # Slow mode is choosen
            mode = 0
            kmer = 5
            option = 'SLOW'
            if rank == 0: print(rank_msg + " Default calculation mode: SLOW")

    sizes = [size1, size2]
    # print("secuencias len: ", sizes)
    if min(sizes) <= 1000 and window is None:  # Secuencias de 1000 bases
        mode = 0  # No se puede hacer con el fast
        #kmer = 0
        window = int(
            min(sizes) * 0.1) + 1  # 1 % de la longitud de la secuencia el +1 es por si el 0.01 % da por debajo de 1
    elif window is None:
        window = int(min(sizes) / 512)

    if mode == 1 and window <= kmer:  # Condicionar que el kmer no puede ser mayor a la ventana, y decir que use otro método.
        if rank == 0:
            print(rank_msg + " window: ", window, " Kmer: ", kmer)
            print(
                rank_msg + " Window value smaller tha k-mer value, please, select other Score mode or modify window value")
        sys.exit(5)

    if joined_file1 != '' and joined_file2 != '':
        #######################################################################
        #### MPI parallel region
        if rank == 0:         
            lengthX, seq1Array = seq2array(joined_file1, window)
            lengthY, seq2Array = seq2array(joined_file2, window)
            print(rank_msg + " finished creating seq Array")
            if lengthX >= lengthY:
                n = lengthX
                strand = 0
            else:
                n = lengthY
                strand = 1
            end_time = time.time()
            print(rank_msg + " module 1 time=", end_time - start_time)
            for m in range(1, threads):
                send_mpi_msg(m,[lengthX, lengthY, strand, n],serialize=True)
        else:
            a = receive_mpi_msg(deserialize=True)
            [lengthX, lengthY, strand, n] = a['data']
        

        #######################################################################
        #### MPI parallel region
        start_time = time.time()
        seqs_per_procs = int(n / (threads-1))
        
        if rank == 0: 
            if strand==0:  # Esto quiere decir que la más larga es la horizontal, por lo tanto vamos a partir en Rank pedazos la vertical y la horizontal queda toda en un archivo.
                namesY,namesX=CreateFiles(seqs_per_procs,seq2Array,seq1Array,lengthY,lengthX,"vertical","horizontal",threads,n)
            else:
                namesX,namesY=CreateFiles(seqs_per_procs,seq1Array,seq2Array,lengthX,lengthY,"horizontal","vertical",threads,n)
            for m in range(1,threads):
                send_mpi_msg(m,[namesX,namesY],serialize=True)

            #En este punto será que ya puedo borrar el joined?################### 
            local_results_dict = {}
            #local_results_dict[0] = local_result
            for r in range(1, threads):
                data_dict = receive_mpi_msg(deserialize=True)
                local_results_dict[data_dict['sender']] = data_dict['data']
            # join all partial results in one
            results = [local_results_dict[r] for r in sorted(local_results_dict.keys())]
        else:
            a=receive_mpi_msg(deserialize=True)
            [namesX,namesY]=a['data']
            # execute sequence alignment using MPI
            local_result = dna_parallel_alignment(namesX,namesY,strand,threads,sim,slen,kmer)
            send_mpi_msg(0, local_result, serialize=True)
            
        end_time = time.time()
        #######################################################################

        if rank == 0:
            
            print(rank_msg + " module 2 time=", end_time - start_time)
            start_time = time.time()
            if strand == 0:
                result = np.row_stack(results)  # horizontal
            else:
                result = np.column_stack(results)  # vertical
            lengthY = result.shape[1]
            lengthX = result.shape[0]
            #print(lengthX,lengthY)
            print(rank_msg + " finished creating matrix")
            end_time = time.time()
            print(rank_msg + " module 3 time=", end_time - start_time)
            # create images with created score matrix
            start_time = time.time()
            graphicalAlignment2(width, height, window, name1, name2, graph, result, lengthX, lengthY, kmer, mode,
                                option, Filter, Name, image_format)
            end_time = time.time()
            print(rank_msg + " module 4 time=", end_time - start_time)
            # deleting joined_files
            start_time=time.time()
            if joined_file1 == joined_file2:
                os.remove(joined_file1)
            else:
                os.remove(joined_file1)
                os.remove(joined_file2)
            
            # deleting vertical and horizontal files

            files = glob.glob("horizontal*") + glob.glob("vertical*") + glob.glob("core.*") 
            for file in files:
                os.remove(file)
            end_time=time.time()
            print("Removing files time= ",end_time - start_time)
        comm.Barrier()
        
    else:
        if rank == 0:
            print(rank_msg + " Please insert the sequences to align")
            print(
                "Usage: -i <Sequence input file1> -a <Sequence input file2>  -h <More info>")
        sys.exit(3)


if __name__ == '__main__':
    main()
