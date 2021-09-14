### IMPORTS ###

import os
import numpy as np
import csv

### SOURCE ###

source = os.path.abspath("./eukaryotes/CSV/euk_sims_expanded.csv")

### MAIN ###

def main():

    raw = open(source).read()
    lines = raw.split('\n')[1:]

    counter = 0
    for line in lines:
        print (line)
        if line != '':
            split = line.split(',')
            taa = float(split[3])
            tga = float(split[6])
            tag = float(split[9])
            if taa > tga and taa > tag:
                counter += 1
    print (counter, len(lines))

### RUN ###

if __name__ == '__main__':
    main()
