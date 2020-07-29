## Graphical alignment of sequences through parallel programming: an approach from the post-genomic era

A graphical alignment or "dot plot" is a method of visual representation of genomic data analysis, commonly used to compare the similarity of two biological sequences. The DOTTER program, developed in 1995, is the most widely used tool for this type of task. The biggest problem with this software is the high runtime for large scale genomic data. GEPARD (2007), performs faster alignments for larger sequences than DOTTER, but reducing the execution time of the alignment of a chromosome against itself, from 382 years with DOTTER to 61 minutes with GEPARD, although with a low level of detail because it uses an approximation method. This article proposes a strategy that works on multiple processors to perform genomic-level alignments in a shorter run time than GEPARD, achieving accelerations up to 27.9 times using 64 processors from the nominal value. The strategy allows the identification of chromosomal rearrangements, repetitive elements, comparison between genomes of different species and the graphic measurement of the assembly quality of genomic sequences quickly. 

### Requirements.

This algorithm requieres Python3 and the following packages. 

1. Numpy
2. Matplotlib
3. Sys
4. Time
5. Multiprocessing
6. Getopt

### How to install requirements

Requirements can be installed under conda environments as follows:

```sh
$ conda create -n G_Align python=3.7
$ conda activate G_Align
$ conda install -c anaconda numpy
$ conda install -c conda-forge matplotlib
$ pip install multiprocessing
```

Other requirements are installed in python 3.7 distribution.

### Usage

To execute the graphical aligner you have to run:

```sh
python3 graphicalAlignment.py -i file1.fasta -a file2.fasta -t 20
```

if you want to define a specific window size, type:
```sh
python3 graphicalAlignment.py -i file1.fasta -a file2.fasta -t 20 -w 500
```
### Example

To create a dot-plot of *Homo sapiens* 21 chromosome against itself write:

```sh
python3 graphicalAlignment.py -i Chr21_HomoSapiens.fa -a Chr21_HomoSapiens.fa -t 20 -w 41990
```
Input file and output image of this line is in Sample_data

### Help:

for extended information please execute `python3 graphicalAlignment.py -h`

## Contact
For more information please write to:
 
- johan.pinad@autonoma.edu.co
- simon.orozco.arias@gmail.com

