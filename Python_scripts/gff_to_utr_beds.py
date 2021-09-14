### IMPORTS ###

import os

### MAIN ###

gff_folder = 'archaea/genome_assemblies_genome_gff/one-per-genus'
out_folder = 'archaea/BED'

def main():

    for root, dirs, filenames in os.walk(gff_folder):
        for f in filenames:

            path = os.path.join(gff_folder, f)
            gff = open(path).read()

            #Get accession
            accession = f.split('.')[0]

            #Split into sections & IGR filter
            lines = [i for i in gff.split('\n') if 'CDS\t' in i]
            filtered_ids = IGR_filter(lines)

            print ('Genes passing the IGR filter...', filtered_ids[:3])

            #For each gene in the filtered list, get exon annotations
            utr_info = get_utrs(lines, filtered_ids)

            #Create files
            bed_path = os.path.join(out_folder, accession + '.bed')
            create_bed(utr_info, bed_path)


### FUNCTIONS ###

def IGR_filter(list):

    ids = []

    #For each set of coordinates...
    for i in list:

        split_list = i.split('\t')

        #Gene info
        chromosome = split_list[0]
        start = split_list[3]
        stop = split_list[4]
        strand = split_list[6]
        gene_id = split_list[8].split(';')[0]

        #Split genes by strand - adjust coordinates to base 1 / reverse complement / to include start and stop codons
        if strand == '+':
            ids.append([gene_id, int(start)-1, int(stop)+3, strand])
        elif strand == '-':
            ids.append([gene_id, int(start)-4, int(stop), strand])

    #Set empty list of filtered ids
    gene_samples = []

    #Calculate 3' UTRs and filter gene ids
    for i,n in enumerate(ids):
        if ids[i][3] == '+':
            if i != len(ids)-1:
                IGR = ids[i+1][1] - ids[i][2]
                if 30 < IGR:
                    gene_samples.append(ids[i][0].replace('ID=gene:', ''))

        elif ids[i][3] == '-':
            if i != 0:
                IGR = ids[i][1] - ids[i-1][2]
                if 30 < IGR:
                    gene_samples.append(ids[i][0].replace('ID=gene:', ''))

    return gene_samples


def get_utrs(list, filtered_ids):

    utr_total = ''

    for entry in list:

        if '\tCDS\t' in entry:

            #Get gene id
            split_list = entry.split('\t')
            gene_id = split_list[8].split(';')[0]

            #Next, utr information
            chromosome = split_list[0]
            start = split_list[3]
            stop = split_list[4]
            strand = split_list[6]

            #Now get coordinates
            if gene_id in filtered_ids:

                if strand == '+':
                    start_co = str(int(stop)-3)
                    stop_co = str(int(start_co)+30)

                elif strand == '-':
                    start_co = str(int(start)-28)
                    stop_co = str(int(start_co)+30)
                else:
                    print ('no strand info!')

                line = "\t".join([chromosome, start_co, stop_co, gene_id, '.', strand, "\n"])
                utr_total += line

    return utr_total


def create_bed(file, bed_path):

    f = open(bed_path, 'a')
    f.write(file)
    f.close()

### RUN ###
if __name__ == '__main__':
    main()
