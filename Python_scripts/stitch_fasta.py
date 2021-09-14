### IMPORTS ###
import os
import re

### FUNCTIONS ###

fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_exon.fa")
new_fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_exon_built.fa")


### FUNCTIONS ###

def main():

    total = ""

    #Define raw data
    raw = open(fasta).read()

    #Build genes and return dictionary
    gene_dict = build_genes(raw)

    #Convert dictionary to fasta
    total = dict_to_fasta(gene_dict)

    #Output fasta
    create_fasta(new_fasta, total)


def build_genes(raw):

    chunks = raw.split('>')
    genes = {}

    for i in chunks:
        if i != '':
            split = i.split('\n')
            fasta_line = split[0]
            line_split = fasta_line.split(':')
            gene_id = line_split[0]
            sequence = split[1]

            if gene_id not in genes:
                genes[gene_id] = sequence
            elif gene_id in genes:
                genes[gene_id] += sequence

    return genes


def dict_to_fasta(dict):

    total = ''

    for key, value in dict.items():
        id = key
        seq = value

        fasta = ">" + id + "\n" + seq + "\n"
        total += fasta

    return total


def create_fasta(new_fasta, total):

    f = open(new_fasta, 'a')
    f.write(total)
    f.close()


### RUN ###
if __name__ == '__main__':
    main()
