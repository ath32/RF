### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### CHOOSE SOURCE FOLDER ###

source = 'eukaryotes/UTR'

### MAIN ###

def main():

    csv_total = []

    stops = ['TAA', 'TGA', 'TAG']

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            if 'utr' in f:

                #Get accession
                accession = f.split('.')[0]

                #Get raw file and access sequences
                path = os.path.join(source, f)
                raw = open(path).read()
                split = [i for i in raw.split('>') if i != '']
                sequences = ["".join(i.split('\n')[1:]).upper() for i in split]
                clean_seqs = [i[len(i)-30:] for i in sequences]
                print (accession, clean_seqs[:3])
                # clean_seqs = [i for i in sequences if i[:3] == 'ATG' and i[len(i)-3:] in stops]

                #Calculate stats
                taa, tga, tag = get_stop_usage(clean_seqs)
                mean_gc, mean_gc3 = get_gc(clean_seqs)
                csv_total.append([accession, taa, tga, tag, mean_gc, mean_gc3])
                print (accession, taa, tga, tag, mean_gc, mean_gc3)

    headers = ['Accession', 'TAA', 'TGA', 'TAG', 'Mean_GC', 'Mean_GC3']
    create_csv(headers, csv_total, 'eukaryote_utr_usage.csv', 'eukaryotes/CSV')

### FUNCTIONS ###

def get_stop_usage(sequences):

    taa = 0
    tga = 0
    tag = 0

    trinucleotides = []

    for sequence in sequences:
      for i in range(len(sequence)-5):
          trinucleotide = sequence[i:i+3]
          trinucleotides.append(trinucleotide)

    for trinuc in trinucleotides:
        if trinuc == 'TAA':
            taa += 1
        elif trinuc == 'TGA':
            tga += 1
        elif trinuc == 'TAG':
            tag += 1

    return taa / (taa+tga+tag), tga / (taa+tga+tag), tag / (taa+tga+tag)

def get_gc(sequences):

    gc_list = []
    gc3_list = []

    for seq in sequences:
        gc_list.append(GC123(seq)[0])
        gc3_list.append(GC123(seq)[3])

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
