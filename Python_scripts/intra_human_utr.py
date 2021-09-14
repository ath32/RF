### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### SOURCE ###

fasta = os.path.abspath("./eukaryotes/UTR_full/Homo_sapiens_utr.fa")

### MAIN ###

def main():

    stops = ['TAA', 'TGA', 'TAG']

    #Define raw data
    raw = open(fasta).read()
    sequences = ["".join(i.split('\n')[1:]) for i in raw.split('>')]
    clean_seqs = [i for i in sequences if i[:3] in stops]

    #Bin genes by GC3
    bins = [[], [], [], [], [], [], [], [], [], []]
    gc3_list = []

    for seq in clean_seqs:
        gc3_list.append(float(GC123(seq)[0]))

    q1 = np.quantile(gc3_list, 0.1)
    q2 = np.quantile(gc3_list, 0.2)
    q3 = np.quantile(gc3_list, 0.3)
    q4 = np.quantile(gc3_list, 0.4)
    q5 = np.quantile(gc3_list, 0.5)
    q6 = np.quantile(gc3_list, 0.6)
    q7 = np.quantile(gc3_list, 0.7)
    q8 = np.quantile(gc3_list, 0.8)
    q9 = np.quantile(gc3_list, 0.9)

    for seq in clean_seqs:
        gc3 = float(GC123(seq)[0])
        if gc3 < q1:
            bins[0].append(seq)
        elif q1 < gc3 < q2:
            bins[1].append(seq)
        elif q2 < gc3 < q3:
            bins[2].append(seq)
        elif q3 < gc3 < q4:
            bins[3].append(seq)
        elif q4 < gc3 < q5:
            bins[4].append(seq)
        elif q5 < gc3 < q6:
            bins[5].append(seq)
        elif q6 < gc3 < q7:
            bins[6].append(seq)
        elif q7 < gc3 < q8:
            bins[7].append(seq)
        elif q8 < gc3 < q9:
            bins[8].append(seq)
        elif q9 < gc3:
            bins[9].append(seq)

    #Analyse stops in each group

    for ind, bin in enumerate(bins):

        mean_gc = get_gc(bin)
        taa, tga, tag = get_stop_usage(bin)

        print ('Group:', int(ind)+1, "3' UTR GC:", mean_gc, 'TAA:', taa, 'TGA:', tga, 'TAG:', tag)


### FUNCTIONS ###

def get_stop_usage(sequences):

    taa = 0
    tga = 0
    tag = 0

    for seq in sequences:

        stop = seq[:3]
        if stop == 'TAA':
            taa += 1
        elif stop == 'TGA':
            tga += 1
        elif stop == 'TAG':
            tag += 1

    return taa / len(sequences), tga / len(sequences), tag / len(sequences)

def get_gc(sequences):

    gc_list = []

    for seq in sequences:
        gc_list.append(GC123(seq)[0])

    return np.mean(gc_list)

### RUN ###
if __name__ == '__main__':
    main()
