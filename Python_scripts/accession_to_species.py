### IMPORTS ###

import os

### SOURCE ###

source = 'archaea/genome_assemblies_genome_fasta/ncbi-genomes-2021-07-01'
cds_source = 'archaea/genome_assemblies_genome_gff/ncbi-genomes-2021-07-01'
out_source = 'archaea/genome_assemblies_genome_gff/one-per-genus'

### FUNCTIONS ###

def main():

    dict = {}
    accessions_list = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            path = os.path.join(source, f)
            raw = open(path).read()
            first_line = raw.split('>')[1].split('\n')[0]
            split = first_line.split(' ')
            accession = split[0]
            name = split[1] + '_' + split[2]
            genus = split[1]
            alt_accession = f.split('.')[0]

            #Build dictionary
            if genus not in dict:
                dict[genus] = [alt_accession]
                accessions_list.append(alt_accession)
            else:
                dict[genus].append(alt_accession)

    print (len(dict), 'species when filtered one per genus...')
    print (dict)

    for root, dirs, filenames in os.walk(cds_source):
        for f in filenames:

            acc_cds = f.split('.')[0]
            path_to_write = os.path.join(cds_source, f)
            path_out = os.path.join(out_source, f)
            raw_out = open(path_to_write).read()
            # first_line = raw.split('>')[1].split('\n')[0]
            # acc_cds = first_line.split('_cds_')[0].replace('lcl|', '').replace(' ', '')

            if acc_cds in accessions_list:
                f = open(path_out, "a")
                f.write(raw_out)
                f.close()

### RUN ###

if __name__ == '__main__':
    main()
