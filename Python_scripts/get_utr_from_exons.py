### IMPORTS ###
import os
import re
import time
import generic as gen

### SOURCE ###

cds_fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_cds_built.fa")
exon_fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_exon_built.fa")
utr_fasta = os.path.abspath("./eukaryotes/UTR_full/Paramecium_tetraurelia_utr.fa")

### MAIN ###

def main():

    total = ""

    #Define raw data
    raw_cds = open(cds_fasta).read()
    raw_exon = open(exon_fasta).read()

    #Find cds in exon fasta
    chunks = raw_cds.split('>')
    exon_chunks = raw_exon.split('>')

    #Get list of exon sequences to search for our cds
    exon_seqs = []

    for i in exon_chunks:
        if i != '':
            split = i.split('\n')
            fasta_line = split[0]
            sequence = split[1].upper()
            exon_seqs.append(sequence)

    #Set empty fasta string
    fasta = ''
    count = 0

    #Search each cds against exon list
    for j in chunks:
        start = time.time()

        if j != '':
            cds_split = j.split('\n')
            cds_fasta_line = cds_split[0]
            cds_sequence = cds_split[1].upper()

            for k in exon_seqs:
                if cds_sequence in k:
                    match = k
                    start_coord = match.find(cds_sequence)
                    end_coord = int(start_coord) + len(cds_sequence)

                    three_prime_start = end_coord
                    three_prime_end = len(k)
                    three_prime_utr = k[three_prime_start: three_prime_end]

                    #Stop codon could be either of these
                    stop_codon = cds_sequence[len(cds_sequence)-3:]
                    # stop_codon = three_prime_utr[0:3]
                    print (stop_codon)

                    if stop_codon == 'TAA' or stop_codon == 'TGA' or stop_codon == 'TAG':
                        if nucleotide_check(three_prime_utr) == True:
                            new_section = ">" + cds_fasta_line + '\n' + stop_codon + three_prime_utr + '\n'
                            fasta += new_section
                            count += 1
                            stop = time.time()
                            total_time = stop - start
                            print (count, total_time)

    f = open(utr_fasta, 'a')
    f.write(fasta)
    f.close()

### FUNCTIONS ###

def nucleotide_check(sequence):
    filtered = re.findall('[ACTG]', sequence)
    f_string = "".join(filtered)

    if f_string == sequence:
        nucleotides = True
    else:
        nucleotides = False
        print ('failed nucleotide check')

    return nucleotides

### RUN ###
if __name__ == '__main__':
    main()
