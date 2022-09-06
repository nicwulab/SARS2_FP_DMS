#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def fitness_calculate(df, rep, freq_cutoff):
    df['fit_'+rep] = np.log10(df[rep+'_P1_freq']/df[rep+'_ipt_freq'])
    df_high_freq = df[df['avg_ipt_freq'] >= freq_cutoff]
    fit_summary = df_high_freq.groupby('mut_class')['fit_'+rep].mean()
    fit_summary = fit_summary.reset_index()
    fit_silent   = (float(fit_summary.loc[fit_summary['mut_class']=='silent']['fit_'+rep]))
    fit_nonsense = (float(fit_summary.loc[fit_summary['mut_class']=='nonsense']['fit_'+rep]))
    fit_WT       = (float(fit_summary.loc[fit_summary['mut_class']=='WT']['fit_'+rep]))
    df['fit_'+rep] = (df['fit_'+rep]-fit_nonsense)/(fit_WT-fit_nonsense)
    return (df)

def wrapper(count_file, freq_cutoff):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname:
            df = count_to_freq(df, colname)
    df['avg_ipt_freq'] = (df['Lib1_ipt_freq']+df['Lib2_ipt_freq'])/2
    df = fitness_calculate(df, 'Lib1', freq_cutoff)
    df = fitness_calculate(df, 'Lib2', freq_cutoff)
    df['fit'] = (df['fit_Lib1'] + df['fit_Lib2'])/2
    return (df)

def main():
    freq_cutoff = 0
    outfile_1 = "result/FP_DMS_fit.tsv"
    outfile_2 = "result/FP_DMS_fit_by_resi.tsv"
    df = wrapper("result/FP_DMS_count_aa.tsv", freq_cutoff)
    print ('writing: %s' % outfile_1)
    df.to_csv(outfile_1, sep="\t", index=False)
    df['resi'] = df['mut'].str[0:-1]
    all_resi   = df[df['mut_class'] != 'WT'].groupby('resi').size().reset_index()
    df_by_resi = df
    df_by_resi = df_by_resi[df_by_resi['avg_ipt_freq'] >= freq_cutoff]
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'WT']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'silent']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'nonsense']
    df_by_resi_mean  = df_by_resi.groupby('resi')['fit'].mean().reset_index(name='mean_fit')
    df_by_resi_count = df_by_resi.groupby('resi').size().reset_index(name='count')
    df_by_resi = pd.merge(df_by_resi_mean, df_by_resi_count, on='resi', how='outer')
    df_by_resi = pd.merge(all_resi, df_by_resi, on='resi', how='outer')
    df_by_resi = df_by_resi.sort_values(by='resi', key=lambda x:x.str[1::].astype(int))
    df_by_resi['pos'] = df_by_resi['resi'].str[1::].astype(int)
    df_by_resi = df_by_resi[['resi','pos','count','mean_fit']]
    print ('writing: %s' % outfile_2)
    df_by_resi.to_csv(outfile_2, sep="\t", index=False)

if __name__ == "__main__":
    main()
