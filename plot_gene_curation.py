import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

with open('test.tsv', 'r') as file:
    data = pd.read_csv('test.tsv', delim_whitespace = True)
    sns.set_style('whitegrid')
    plt.rcParams["font.family"] = "Times New Roman"

    genes = data['gene'].unique()

    plt.figure(figsize=(8,8))

    for gene in genes:

        cur_tc = data.loc[(data['status'] == 'curated') & (data['gene'] == gene)]['rel_tree_certainty'].iloc[0]
        raw_tc = data.loc[(data['status'] == 'raw') & (data['gene'] == gene)]['rel_tree_certainty'].iloc[0]

        cur_stdev = data.loc[(data['status'] == 'curated') & (data['gene'] == gene)]['stdev_raw_length'].iloc[0]
        raw_stdev = data.loc[(data['status'] == 'raw') & (data['gene'] == gene)]['stdev_raw_length'].iloc[0]

        plt.arrow(raw_stdev, raw_tc, cur_stdev-raw_stdev, cur_tc-raw_tc, color = 'darkgrey', width = 0.0015, length_includes_head = True)

        plt.scatter(raw_stdev, raw_tc, color = 'turquoise', s = 100)
        plt.annotate(gene + ' raw', (raw_stdev + 0.01, raw_tc + 0.01))
        plt.scatter(cur_stdev, cur_tc, color = 'steelblue', s = 100)
        plt.annotate(gene + ' curated', (cur_stdev + 0.01, cur_tc + 0.01))

    plt.ylabel('Tree certainty', fontsize = 16)
    plt.xlabel('Protein length σ', fontsize = 16)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.savefig('plot_tc_stdev.png', dpi = 1000)
