### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### SOURCE ###

a_source = 'archaea/FASTA'
b_source = 'bacteria/UTRs'

### MAIN ###

def main():

    csv_total = []

    stops = ['TAA', 'TGA', 'TAG']

    bacteria_gc_dict = {}
    bacteria_tag_dict = {}

    for root, dirs, filenames in os.walk(b_source):
        for f in filenames:

            #Get accession
            bacteria_accession = f.split('.')[0]

            #Get raw file and access sequences
            b_path = os.path.join(b_source, f)
            raw_b = open(b_path).read()
            split_b = [i for i in raw_b.split('>') if i != '']
            b_sequences = ["".join(i.split('\n')[1:]).upper() for i in split_b]
            clean_seqs_b = [i for i in b_sequences if i[:3] in stops and 'N' not in i]

            if len(clean_seqs_b) > 0:

                #Calculate stats
                taa_b, tga_b, tag_b = get_stop_usage(clean_seqs_b)
                mean_gc_b, mean_gc3_b = get_gc(clean_seqs_b)

                if bacteria_accession not in bacteria_gc_dict:
                    bacteria_gc_dict[bacteria_accession] = mean_gc_b
                    print ('bacteria gc added...')
                else:
                    continue

                if bacteria_accession not in bacteria_tag_dict:
                    bacteria_tag_dict[bacteria_accession] = tag_b
                    print ('bacteria tag added...')
                else:
                    continue

    print ('length of bacterial dict = ', len(bacteria_gc_dict))

    for root, dirs, filenames in os.walk(a_source):
        for f in filenames:

            #Get accession
            archaea_accession = f.split('.')[0]

            #Get raw file and access sequences
            a_path = os.path.join(a_source, f)
            raw_a = open(a_path).read()
            split_a = [i for i in raw_a.split('>') if i != '']
            a_sequences = ["".join(i.split('\n')[1:]).upper() for i in split_a]
            clean_seqs_a = [i for i in a_sequences if i[:3] in stops and 'N' not in i]

            if len(clean_seqs_a) > 0:

                #Calculate stats
                taa_a, tga_a, tag_a = get_stop_usage(clean_seqs_a)
                mean_gc_a, mean_gc3_a = get_gc(clean_seqs_a)
                print ('GC to match = ', mean_gc_a)

                #Select bacterial genome
                res_key, res_val = min(bacteria_gc_dict.items(), key=lambda x: abs(mean_gc_a - x[1]))
                print ('match to', archaea_accession, '=', res_key, res_val)

                csv_total.append([archaea_accession, res_key, mean_gc_a, res_val, tag_a, bacteria_tag_dict[res_key]])

    headers = ['Archaea_accession', 'Bacteria_accession', 'Archaea_GC', 'Bacteria_GC', 'Archaea_TAG', 'Bacteria_TAG']
    create_csv(headers, csv_total, 'archaea_test.csv', 'bacteria/CSV')

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
    gc3_list = []

    for seq in sequences:
        gc_list.append(GC123(seq[3:])[0])
        gc3_list.append(GC123(seq[3:])[3])

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
