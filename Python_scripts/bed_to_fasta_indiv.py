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

cds_bed = os.path.abspath("./eukaryotes/BED/Paramecium_tetraurelia.ASM16542v1.51_cds.bed")
cds_fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_cds.fa")
exon_bed = os.path.abspath("./eukaryotes/BED/Paramecium_tetraurelia.ASM16542v1.51_exon.bed")
exon_fasta = os.path.abspath("./eukaryotes/UTR/Paramecium_tetraurelia_exon.fa")
genome_fasta = os.path.abspath("./eukaryotes/WGS/Paramecium_tetraurelia.ASM16542v1.dna.toplevel.fa")

### MAIN ###

def main():

    fasta_from_intervals(exon_bed, exon_fasta, genome_fasta, names = True)
    fasta_from_intervals(cds_bed, cds_fasta, genome_fasta, names = True)
    
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
