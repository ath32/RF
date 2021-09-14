### CREDITS ROSINA SAVISAAR AND LIAM ABRAHAMS ###

### IMPORTS ###
import generic as gen
import re
import collections
import copy
import numpy as np
import os
from pathlib import Path
import random
import shutil

### SOURCE ###

intron_bed_source = 'archaea/BED'
intron_fasta_source = 'archaea/FASTA'
genome_fasta = 'archaea/genome_assemblies_genome_fasta/ncbi-genomes-2021-07-01'

### MAIN ###

def main():

    for root, dirs, filenames in os.walk(intron_bed_source):
        for f in filenames:

            accession = f.split('.')[0]
            bed_path = os.path.join(intron_bed_source, f)

            for root, dirs, filenames in os.walk(genome_fasta):
                for x in filenames:

                        if accession in x:

                            fasta_path = os.path.join(genome_fasta, x)
                            utr_fasta_path = os.path.join(intron_fasta_source, accession + '.fa')

                            fasta_from_intervals(bed_path, utr_fasta_path, fasta_path, names = True)

### FUNCTIONS ###

def fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = False):
    '''
    Takes a bed file and creates a fasta file with the corresponding sequences.
    If names == False, the fasta record names will be generated from the sequence coordinates.
    If names == True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    '''

    #if the index file exists, check whether the expected features are present
    genome_fasta_index = genome_fasta + '.fai'
    if(os.path.exists(genome_fasta_index)):
        bed_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(bed_file, "\t")])))
        index_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(genome_fasta_index, "\t")])))
        if(not set(bed_chrs).issubset(set(index_chrs))):
            gen.remove_file(genome_fasta_index)

    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    gen.run_process(bedtools_args)
    names, seqs = gen.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    gen.write_to_fasta(names, seqs, fasta_file)


### RUN ###
if __name__ == '__main__':
    main()
