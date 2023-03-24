import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def setup_dfs():
    FP_DMS_fit_df = pd.read_table('result/FP_DMS_fit.tsv')
    FP_DMS_fit_df = FP_DMS_fit_df[['mut', 'mut_class', 'fit_Calu3_noAb', 'fit_E6_noAb']]

    BA1_DMS_df = pd.read_csv('data/BA1_DMS_muteffects_observed.csv')
    BA1_DMS_df = BA1_DMS_df[['reference_site', 'wildtype', 'mutant', 'effect']]

    Tree_aa_fitness_df = pd.read_csv('data/Tree_aa_fitness.csv')
    Tree_aa_fitness_df = Tree_aa_fitness_df.rename(columns={'aa_site': 'reference_site'})
    Tree_aa_fitness_df = Tree_aa_fitness_df.rename(columns={'aa': 'mutant'})

    FP_DMS_fit_df = FP_DMS_fit_df[FP_DMS_fit_df['mut_class'] == 'missense']
    FP_DMS_fit_df[['wildtype', 'reference_site', 'mutant']] = FP_DMS_fit_df['mut'].str.extract(r'^([A-Z])(\d+)([A-Z])$')
    mask = FP_DMS_fit_df['wildtype'] != FP_DMS_fit_df['mutant']
    FP_DMS_fit_df = FP_DMS_fit_df[mask]

    BA1_merged_df = pd.merge(FP_DMS_fit_df, BA1_DMS_df, on=['reference_site', 'wildtype', 'mutant'], how='inner')

    FP_DMS_fit_df['reference_site'] = FP_DMS_fit_df['reference_site'].astype(str)
    Tree_aa_fitness_df['reference_site'] = Tree_aa_fitness_df['reference_site'].astype(str)

    tree_FP_DMS_fit_merged = pd.merge(FP_DMS_fit_df, Tree_aa_fitness_df, on=['reference_site', 'mutant'], how='inner')

    BA1_DMS_df['reference_site'] = BA1_DMS_df['reference_site'].astype(str)

    return BA1_merged_df, tree_FP_DMS_fit_merged


def plot(df, x, y, title, xlabel, ylabel, output_file):

    plt.rc('font', size=7) 
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    a4_dims = (1.2, 1.2)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.scatterplot(data=df, x=x, y=y, alpha=0.5, s=5, color='#4D4D4D', linewidth=0, marker='o').set(
        title=title)
    plt.xlabel(xlabel, fontfamily='Arial', fontsize=7, fontweight='bold')
    plt.ylabel(ylabel, fontfamily='Arial', fontsize=7, fontweight='bold')
    plt.title('')
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines['bottom'].set_color('#000000')
    ax.spines['left'].set_color('#000000')
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.clf()


def main():

    BA1_merged_df, tree_FP_DMS_fit_merged = setup_dfs()
    plot(BA1_merged_df, "fit_Calu3_noAb", "effect", 'BA1_DMS_Fit',
         'fitness in Calu3 (noAb)', 'BA1 functional score\n(Dadonaite et al.)','graph/cor_BA1_vs_Calu3.png')
    plot(BA1_merged_df, "fit_E6_noAb", "effect", 'BA1_DMS_Fit',
         'fitness in E6 (noAb)', 'BA1 functional score\n(Dadonaite et al.)', 'graph/cor_BA1_vs_E6.png')
    plot(tree_FP_DMS_fit_merged, "fit_Calu3_noAb", "fitness", 'tree_FP_DMS_fit',
         'fitness in Calu3 (noAb)', 'fitness from phylogeny\n(Bloom & Neher)', 'graph/cor_tree_vs_Calu3.png')
    plot(tree_FP_DMS_fit_merged, "fit_E6_noAb", "fitness", 'tree_FP_DMS_fit',
         'fitness in E6 (noAb)', 'fitness from phylogeny\n(Bloom & Neher)', 'graph/cor_tree_vs_E6.png')
    print (BA1_merged_df.corr(method='pearson'))
    print (tree_FP_DMS_fit_merged.corr(method='pearson'))

if __name__ == "__main__":
  main()
