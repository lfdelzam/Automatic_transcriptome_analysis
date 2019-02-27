#!/usr/bin/env python3
import argparse
import sys

parser = argparse.ArgumentParser(description="the program reads a file containing the name of the genes more expressed, and a gff file. Thier correspondig CDSs positions in the genome\
are extracted and search the positions in the genome file and produce a fasta file with the sequences")
parser.add_argument ('-g', dest= 'gff', help='gff file', required = True)
parser.add_argument ('-i', dest= 'input', help='input file containing the genome, type .fna', required = True)
parser.add_argument ('-t', dest= 'top', help='input file containing the top genes more expressed', required = True)
parser.add_argument ('-o', dest= 'out', help='output fasta file containing the list of extracted sequences', required = True)
args = parser.parse_args()

cds_dict = {}
genome = ''
gene_list = []

#reading the file containing the name of the genes
with open(args.top, 'r') as fin, open(args.gff, 'r') as fcds:
    for line in fin:
        line = line.rstrip()
        if not line.startswith("baseMean"):
            line = line.split("\t")
            gene_list.append(line[0]) #list of gene names
#reading the file conataining thier correspondig CDS positions in the genome
    for line3 in fcds:
        line3 = line3.rstrip().split()
        if line3[2] == 'CDS':
            id = line3[9][1:-2] #name of the gene
            if id in gene_list:
      # Storing CDSs positions in the genome for each gene (start and stop)
                if id in cds_dict:
                    cds_dict[id] += ","+line3[3]
                    cds_dict[id] += ","+line3[4]
                else:
                    cds_dict[id] = line3[3]
                    cds_dict[id] += ","+line3[4]

# Reading the genome file
with open(args.input, 'r') as fin2:
    for line2 in fin2:
        line2 = line2.rstrip()
        if not line2.startswith('>'):
            genome += line2 #Genome stored in a single line
#searching the positions in the genome file and producing a fasta file with the sequences
with open (args.out, 'w') as fout:
    for key in cds_dict:
        print('>'+key, file=fout)
        value = cds_dict[key]
        value = value.split(",") #to have access to each CDS position value
        for i in range(0,len(value),2):
            start = int(value[i])-1 #
            stop = int(value[i+1])+1
            print(genome[start:stop], end='', file=fout)
        print('', file=fout) #ading new line to the last line of the gene sequence
