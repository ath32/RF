### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### SOURCE ###

source = 'eukaryotes/UTR'

### MAIN ###

def main():

    csv_total = []

    stops = ['TAA', 'TGA', 'TAG']

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            if 'utr' in f and 'brucei' not in f:

                #Get accession
                accession = f.split('.')[0]
                print (accession)

                #Get raw file and access sequences
                path = os.path.join(source, f)
                raw = open(path).read()
                split = [i for i in raw.split('>') if i != '']
                sequences = ["".join(i.split('\n')[1:]).upper() for i in split]
                clean_seqs = [i[len(i)-30:] for i in sequences]

                #Split into codons
                filtered_codons = []
                filtered_sequences = []

                for sequence in clean_seqs:

                    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
                    codons = codons[1:]

                    #Get sequences after any ASC
                    if 'TAA' in codons:
                        taa_index = codons.index('TAA')
                    else:
                        taa_index = 0

                    if 'TGA' in codons:
                        tga_index = codons.index('TGA')
                    else:
                        tga_index = 0

                    if 'TAG' in codons:
                        tag_index = codons.index('TAG')
                    else:
                        tag_index = 0

                    if taa_index + tga_index + tag_index > 0:

                        indexes = [taa_index, tga_index, tag_index]
                        print (indexes)
                        indexes = [i for i in indexes if i != 0]
                        smallest = min(indexes)
                        new_codons = codons[smallest+1:]
                        if len(new_codons) > 1:
                            for codon in new_codons:
                                filtered_codons.append(codon)
                            new_sequences = "".join(new_codons)
                            filtered_sequences.append(new_sequences)

                #Calculate stats
                taa, tga, tag = get_stop_usage(filtered_codons)
                mean_gc, mean_gc3 = get_gc(filtered_sequences)
                csv_total.append([accession, taa, tga, tag, mean_gc, mean_gc3])
                print (accession, taa, tga, tag, mean_gc, mean_gc3)

    headers = ['Accession', 'TAA', 'TGA', 'TAG', 'Mean_GC', 'Mean_GC3']
    create_csv(headers, csv_total, 'eukaryote_utr_filtered_usage.csv', 'eukaryotes/CSV')

### FUNCTIONS ###

def get_stop_usage(codons):

    taa = 0
    tga = 0
    tag = 0

    for codon in codons:

        if codon == 'TAA':
            taa += 1
        elif codon == 'TGA':
            tga += 1
        elif codon == 'TAG':
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
