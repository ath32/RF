### IMPORTS ###
import os

### SOURCE ###

example = os.path.abspath("./3_Cterm/Homo_sapiens/Ensembl/Homo_sapiens.GRCh38.99.gff3")

### MAIN ###

def get_isoforms_from_gff(gff_path):

    #Create list to contain the longest transcripts
    transcripts = []
    #Define raw data
    gff = open(gff_path).read()
    #Split it into gene sections
    split = gff.split('###')
    for i in split:
        dict = {}
        lines = i.split('\n')
        exon_lines = [line for line in lines if '\texon\t' in line]
        if len(exon_lines) > 0:
            for exon_bit in exon_lines:
                exon_split = exon_bit.split('\t')
                start = exon_split[3]
                stop = exon_split[4]
                extra = exon_split[8]
                extra_split = extra.split(';')
                transcript = extra_split[0].replace('Parent=transcript:', '')
                if transcript in dict:
                    dict[transcript] += (int(stop) - int(start))
                else:
                    dict[transcript] = (int(stop) - int(start))
        if len(dict) > 0:
            longest = max(dict, key=dict.get)
            transcripts.append(longest)
    return transcripts
