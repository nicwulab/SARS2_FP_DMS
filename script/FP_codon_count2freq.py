#!/usr/bin/python
import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def wrapper(count_file):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'pos' not in colname and 'codon' not in colname:
            df = count_to_freq(df, colname)
    samples = [colname.replace('Rep1_','') for colname in colnames if 'Rep1' in colname]
    for sample in samples:
      df[sample+'_freq'] = (df['Rep1_'+sample+'_freq']+df['Rep2_'+sample+'_freq'])/2
    return (df)

def main():
    outfile = "result/FP_DMS_codon_freq.tsv"
    df = wrapper("result/FP_DMS_count_codon.tsv")
    print ('writing: %s' % outfile)
    df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
