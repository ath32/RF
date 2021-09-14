### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### CHOOSE SOURCE FOLDER ###

source = 'eukaryotes/FASTA'

### MAIN ###

def main():

    csv_total = []

    stops = ['TAA', 'TGA', 'TAG']

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get accession
            accession = f.split('.')[0]
            print (accession)

            #Get raw file and access sequences
            path = os.path.join(source, f)
            raw = open(path).read()
            split = [i for i in raw.split('>') if i != '']
            sequences = ["".join(i.split('\n')[1:]).upper() for i in split]
            clean_seqs = [i for i in sequences if i[len(i)-3:] in stops and 'N' not in i]
            if len(clean_seqs) > 0:

                #Calculate stats
                taa, tga, tag = get_stop_usage(clean_seqs)
                mean_gc, mean_gc3 = get_gc(clean_seqs)
                csv_total.append([accession, taa, tga, tag, mean_gc, mean_gc3])
                print (accession, taa, tga, tag, mean_gc, mean_gc3)

    headers = ['Accession', 'TAA', 'TGA', 'TAG', 'Mean_GC', 'Mean_GC3']
    create_csv(headers, csv_total, 'eukaryote_usage_20.csv', 'eukaryotes/CSV')

### FUNCTIONS ###

def get_stop_usage(sequences):

    taa = 0
    tga = 0
    tag = 0

    for seq in sequences:

        stop = seq[len(seq)-3:]
        if stop == 'TAA':
            taa += 1
        elif stop == 'TGA':
            tga += 1
        elif stop == 'TAG':
            tag += 1

    return taa / len(sequences), tga / len(sequences), tag / len(sequences)

def get_gc(sequences):

    gc_list = []
    gc3_list = []

    for seq in sequences:
        gc_list.append(GC123(seq[:len(seq)-3])[0])
        gc3_list.append(GC123(seq[:len(seq)-3])[3])

    return np.mean(gc_list), np.mean(gc3_list)

def create_csv(headers, csv_total, filename, subdir):

    filepath = os.path.join(subdir, filename)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
