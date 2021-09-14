### IMPORTS ###

import os
import numpy as np
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### CHOOSE SOURCE FOLDER ###

cds_source = 'bacteria/UTRs'

### MAIN ###

def main():

    stops = ['TAA', 'TGA', 'TAG']

    #Parse folder
    for root, dirs, filenames in os.walk(cds_source):
        for f in filenames:

            csv_total = []

            #Get accession
            # accession = f.split('_')[0] + '_' + f.split('_')[1]
            accession = f.split('.')[0]

            #Get raw file and access sequences
            path = os.path.join(cds_source, f)
            raw = open(path).read()
            split = [i for i in raw.split('>') if i != '']

            #Parse genes
            for i in split:

                info = i.split('\n')[0]
                gene_id = info.split(';')[-2].replace(' ', '').replace(';', '')
                sequence = "".join(i.split('\n')[1:]).upper()
                if len(sequence) > 0:
                    stop = sequence[0:3]
                    gc3 = GC123(sequence)[0]
                    if stop in stops:
                        if stop == 'TAA':
                            csv_total.append([gene_id, gc3, 1, 0, 0])
                        elif stop == 'TGA':
                            csv_total.append([gene_id, gc3, 0, 1, 0])
                        elif stop == 'TAG':
                            csv_total.append([gene_id, gc3, 0, 0, 1])

            headers = ['Gene_id', 'GC3', 'TAA', 'TGA', 'TAG']
            create_csv(headers, csv_total, accession + '_logreg.csv', 'bacteria/Intragenome/CSV_UTR')
            print (accession, 'done!')


### FUNCTONS ###

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
