### IMPORTS ###
import pandas as pd
import numpy as np
import os
import csv
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import scipy.stats as stats
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

### CHOOSE SOURCE FOLDER ###

source = 'bacteria/Intragenome/CSV_UTR'

### MAIN ###

def main():

    stops = ['TAA', 'TGA', 'TAG']

    csv_total = []

    #Parse folder
    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            if 'results' not in f:

                try:

                    #Accession
                    acc = f.split('_')[0]
                    print (acc)

                    #Get raw file and access sequences
                    path = os.path.join(source, f)
                    data = pd.read_csv(path, header = 0)
                    data = data.set_index('Gene_id')
                    data['intercept'] = 1

                    #Choose columns
                    y1 = data['TAA']
                    y2 = data['TGA']
                    y3 = data['TAG']
                    X = data[['GC3', 'intercept']]

                    #Model
                    logit_model1=sm.Logit(y1, X).fit()
                    coefficient1 = logit_model1.params[0]
                    pval1 = logit_model1.pvalues[0]

                    logit_model2=sm.Logit(y2, X).fit()
                    coefficient2 = logit_model2.params[0]
                    pval2 = logit_model2.pvalues[0]

                    logit_model3=sm.Logit(y3, X).fit()
                    coefficient3 = logit_model3.params[0]
                    pval3 = logit_model3.pvalues[0]

                    #Calculate standard dev
                    std= np.std(data['GC3'])
                    mean_gc3 = np.mean(data['GC3'])

                    csv_total.append([acc, coefficient1, pval1, coefficient2, pval2, coefficient3, pval3, std, mean_gc3])
                    print (acc, 'done!')

                except:
                    pass

    headers = ['Accession', 'TAA_coef', 'TAA_p', 'TGA_coef', 'TGA_p', 'TAG_coef', 'TAG_p', 'STD_GC3', 'Mean_GC3']
    create_csv(headers, csv_total, 'logreg_results_utr.csv', 'bacteria/Intragenome/CSV3')

### FUNCTONS ###

def create_csv(headers, csv_total, filename, subdir):

    filepath = os.path.join(subdir, filename)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
