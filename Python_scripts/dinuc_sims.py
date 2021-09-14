### IMPORTS ###

import os
import numpy as np
import csv
import time
from Bio.Seq import Seq
from Bio.SeqUtils import GC123
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER ###

# source = 'archaea/genome_assemblies_cds_fasta/one-per-genus'
# source = 'bacteria/FASTA'
source = 'eukaryotes/FASTA'

### MAIN ###

stops = ['TAA', 'TGA', 'TAG']

def main():

    filenames = get_files(source)
    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(filenames, ['foo', source], folder_parse, workers=workers)

    csv_total = []

    for process in processes:
        output = process.get()
        csv_total.extend(output)

    write_outputs(csv_total)

### FUNCTIONS ###

def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):

    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example, a list of gene names)
    args: a list of arguments to the function. Put in "foo" in place of the argument you are parallelizing over.
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''

    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for elem in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = elem
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)


def get_files(source):

    ''' Obtain all file names from the target folder '''

    files = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            files.append(f)

    return files


def folder_parse(filenames, source):

    ''' Function to parallelise '''

    csv_total = []

    for i,f in enumerate(filenames):

        print ('Doing genome {0}/{1}'.format(i+1, len(filenames)))

        start_time = time.time()

        #Get accession
        accession = f.split('.')[0]

        #Get raw file and access sequences
        path = os.path.join(source, f)
        raw = open(path).read()
        split = [i for i in raw.split('>') if i != '']
        sequences = ["".join(i.split('\n')[1:]).upper() for i in split]
        clean_seqs = [i for i in sequences if i[:3] == 'ATG' and i[len(i)-3:] in stops]

        #Calculate observed
        taa, tga, tag = get_stop_usage(clean_seqs)
        mean_gc, mean_gc3 = get_gc(clean_seqs)

        #Calculate expected
        total_sequence = ''
        for seq in clean_seqs:
            total_sequence += seq

        taa_e, tga_e, tag_e = get_simulations(clean_seqs, total_sequence, start_time)

        #Calculate O-E/E
        if taa_e > 0 and tga_e > 0 and tag_e > 0:
            print (taa, taa_e, tga, tga_e, tag, tag_e)
            taa_dev = (taa - taa_e) / taa_e
            tga_dev = (tga - tga_e) / tga_e
            tag_dev = (tag - tag_e) / tag_e

            csv_total.append([accession, taa, taa_e, taa_dev, tga, tga_e, tga_dev, tag, tag_e, tag_dev, mean_gc, mean_gc3])
        else:
            print (accession)

    return csv_total


def write_outputs(csv_total):

    ''' Define headers '''

    #Set headers for the CSV
    headers = ["Accession",
    "Observed_TAA", "Expected_TAA", "TAA_Dev",
    "Observed_TGA", "Expected_TGA", "TGA_Dev",
    "Observed_TAG", "Expected_TAG", "TAG_Dev",
    "Mean_GC", "Mean_GC3"
    ]

    create_csv(headers, csv_total)

def get_stop_usage(sequences):

    taa = 0
    tga = 0
    tag = 0

    for seq in sequences:

        stop = seq[len(seq)-3:]
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
        gc_list.append(GC123(seq)[0])
        gc3_list.append(GC123(seq)[3])

    return np.mean(gc_list), np.mean(gc3_list)

def get_stop_trinuc_usage(sequences):

    taa = 0
    tga = 0
    tag = 0

    trinucleotides = []

    for sequence in sequences:
        trinucleotides.append(sequence[0])

    for trinuc in trinucleotides:
        if trinuc == 'TAA':
            taa += 1
        elif trinuc == 'TGA':
            tga += 1
        elif trinuc == 'TAG':
            tag += 1

    print (taa, tga, tag)

    return taa / (taa+tga+tag), tga / (taa+tga+tag), tag / (taa+tga+tag)


def get_simulations(sequence_list, total_sequence, start_time):

    ''' Simulate 10,000 null sequences for each gene based upon dinucleotide content '''

    #Probability of first base - dictionary
    total_sequence = total_sequence.replace('N', '')
    nt_dict = {}
    length = len(total_sequence)

    for i in range(length):
        base = total_sequence[i]
        if base not in nt_dict:
            nt_dict[base] = 1
        else:
            nt_dict[base] += 1

    nt_freq = {k: v / length for k, v in nt_dict.items()}

    print ('nucleotide dictionary done...')

    #Second base / Next base
    trans = {}
    for i in range(len(total_sequence)-1):

        dinuc = total_sequence[i:i+2]
        first_dinuc = dinuc[0]
        sec_dinuc = dinuc[1]

        if first_dinuc not in trans:
            trans[first_dinuc] = [sec_dinuc]
        else:
            trans[first_dinuc] += sec_dinuc

    for key, value in trans.items():
        trans[key] = [(value.count('A') / len(value)), (value.count('C') / len(value)), (value.count('G') / len(value)), 1 - (value.count('A') / len(value)) - (value.count('C') / len(value)) - (value.count('G') / len(value))]

    print ('dinucleotide dictionary done...')

    #Generate simulations
    total_list = []
    count = 0

    #This while loop determines how many simulations will be completed
    while count < 5000:

        sim_n = []
        codon_list = []

        sim = []

        np.random.seed()

        #Calculate next base
        first_base = np.random.choice(['A', 'C', 'G', 'T'], p=[nt_freq['A'], nt_freq['C'], nt_freq['G'], (1 - nt_freq['A'] - nt_freq['C'] - nt_freq['G'])])
        sim.append(first_base)

        #For one gene simulation
        while (len(sim) < 3):

            #Generate next base
            prev_base = sim[-1]
            next_base = np.random.choice(['A', 'C', 'G', 'T'], p = [trans[prev_base][0], trans[prev_base][1], trans[prev_base][2], trans[prev_base][3]])
            sim.append(next_base)

        sim = "".join(sim)
        total_list.append([sim])
        end_time = time.time()
        total_time = end_time - start_time
        print (str(count) + '/10,000 simulations on this thread')
        count += 1

    taa, tga, tag = get_stop_trinuc_usage(total_list)

    return taa, tga, tag


def create_csv(headers, csv_total):

    filename = "euk_sims_expanded.csv"
    subdir = "eukaryotes/CSV"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
