"""
:File: ProcessTDT.py
:Author: Jack A. Kosmicki
:Last updated: 2014-09-08

read output file from TDT.py

Usage:
    ProcessTDT.py doTDTstat <inputFile> <outputFile_Name> [options]

Options:
    --genes=G_FILE  file with names of the genes you want to look at
    --norm          flag to compute a normalized variance and chi-squared stat
                    default is false
    --M=Val         Missingness threshold [default: 0.99]
    -h, --help      Print this message and exit.
    -v, --version   Print the version and exit.
"""

import sys
import numpy as np
import pandas as pd
from collections import OrderedDict
from scipy import stats
from docopt import docopt

def read_calcTDT(inputFile, outputFile_Name, missThresh):

    df = pd.read_csv(inputFile, sep='\t', header=0)
    # remove all rows that begin with CHROM
    mask = df.applymap(lambda x: x in ['CHROM', 'REF'])
    df = df[-mask.any(axis=1)]
    df['AN'] = df['AN'].apply(int)
    df['transmitted'] = df['transmitted'].apply(int)
    df['untransmitted'] = df['untransmitted'].apply(int)

    max = df['AN'].max(axis=0)

    df = df.ix[df.index[df['AN'] > max * missThresh]]  # check for missingness
    df = df.ix[df.index[df['transmitted'] + df['untransmitted'] != 0]]


    df['TDT_SCORE'] = ((df['transmitted'] - df['untransmitted'])**2) / (df['transmitted'] + df['untransmitted'])
    df['TDT_PVAL'] = 1- stats.chi2.cdf(df['TDT_SCORE'], 1)

    df.to_csv(outputFile_Name, sep='\t', index=False)

if __name__ == "__main__":
    args = docopt(__doc__, version='0.1')

    missThresh = float(args['--M'])   # default is 0.99

    if args['doTDTstat']:
        read_calcTDT(args['<inputFile>'], args['<outputFile_Name>'], missThresh)