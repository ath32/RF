### IMPORTS ###
import os
import get_isoforms_from_gff as gio

### SOURCE ###

gff_file = os.path.abspath("./eukaryotes/GFF/Paramecium_tetraurelia.ASM16542v1.51.gff3")
exon_path = os.path.abspath("./eukaryotes/BED/Paramecium_tetraurelia.ASM16542v1.51_exon.bed")
cds_path = os.path.abspath("./eukaryotes/BED/Paramecium_tetraurelia.ASM16542v1.51_cds.bed")

### MAIN ###

def main():

    #Define raw data
    gff = open(gff_file).read()

    #Split into sections & IGR filter
    lines = gff.split('\n')
    filtered_ids = IGR_filter(lines)

    #Get isoforms
    transcripts_list = gio.get_isoforms_from_gff(gff_file)

    #For each gene in the filtered list, get exon annotations
    chunks = gff.split('###')
    cds_info = get_cds(chunks, filtered_ids, transcripts_list)
    exon_info = get_exons(chunks, filtered_ids, transcripts_list)

    #Create files
    create_bed(exon_info, exon_path)
    create_bed(cds_info, cds_path)

### FUNCTIONS ###

def IGR_filter(list):

    ids = []

    #For each set of coordinates...
    for i in list:
        split_list = i.split(';')

        #Split to access gene info
        if '\tgene\t' in split_list[0] and 'Mito' not in split_list[0]:
            chunk1 = split_list[0].strip()
            chunk2 = chunk1.split('\t')

            #Gene info
            chromosome = chunk2[0]
            start = chunk2[3]
            stop = chunk2[4]
            strand = chunk2[6]
            gene_id = chunk2[8]

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
                print (IGR)
                if 0 < IGR:
                    gene_samples.append(ids[i][0].replace('ID=gene:', ''))

        elif ids[i][3] == '-':
            if i != 0:
                IGR = ids[i][1] - ids[i-1][2]
                print (IGR)
                if 0 < IGR:
                    gene_samples.append(ids[i][0].replace('ID=gene:', ''))

    return gene_samples


def get_exons(list, filtered_ids, transcripts_list):

    exon_total = ''

    for entry in list[1:]:
        if '\tgene\t' in entry:

            #First get gene id
            entry_split = entry.split('ID=gene:')
            next = entry_split[1]
            next_split = next.split(';')
            gene_id = next_split[0]

            #Next, get the exon information
            lines = entry.split('\n')
            exon_lines = [i for i in lines if '\texon\t' in i]

            for exon in exon_lines:
                tab_split = exon.split('\t')
                chromosome = tab_split[0]
                dot = tab_split[5]
                strand = tab_split[6]
                other = tab_split[8]
                other_split = other.split(';')
                transcript = other_split[0].replace('Parent=transcript:', '')

                #Now get coordinates
                if gene_id in filtered_ids and transcript in transcripts_list:

                    if strand == '+':
                        start = str(int(tab_split[3].strip())-1)
                        stop = str(int(tab_split[4].strip()))

                    elif strand == '-':
                        start = str(int(tab_split[3].strip())-1)
                        stop = str(int(tab_split[4].strip()))
                    else:
                        print ('no strand info!')

                    line = "\t".join([chromosome, start, stop, gene_id, dot, strand, "\n"])
                    exon_total += line

    return exon_total


def get_cds(list, filtered_ids, transcripts_list):

    cds_total = ''

    for entry in list[1:]:
        if '\tgene\t' in entry:

            #First get gene id
            entry_split = entry.split('ID=gene:')
            next = entry_split[1]
            next_split = next.split(';')
            gene_id = next_split[0]

            #Next, get the exon information
            lines = entry.split('\n')
            cds_lines = [i for i in lines if '\tCDS\t' in i]

            for cds in cds_lines:
                tab_split = cds.split('\t')
                chromosome = tab_split[0]
                dot = tab_split[5]
                strand = tab_split[6]
                other = tab_split[8]
                other_split = other.split(';')
                transcript = other_split[1].replace('Parent=transcript:', '')

                #Now get coordinates
                if gene_id in filtered_ids and transcript in transcripts_list:

                    if strand == '+':
                        start = str(int(tab_split[3].strip())-1)
                        stop = str(int(tab_split[4].strip()))

                    elif strand == '-':
                        start = str(int(tab_split[3].strip())-1)
                        stop = str(int(tab_split[4].strip()))
                    else:
                        print ('no strand info!')

                    line = "\t".join([chromosome, start, stop, gene_id, dot, strand, "\n"])
                    cds_total += line

    return cds_total


def create_bed(file, bed_path):

    f = open(bed_path, 'a')
    f.write(file)
    f.close()

### RUN ###
if __name__ == '__main__':
    main()
