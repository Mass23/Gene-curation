import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

with open('test.tsv', 'r') as file:
    data = pd.read_csv('test.tsv', delim_whitespace = True)
    print(data)

    genes = data['gene'].unique()
    # Relative tree certainty improvement = TC_cur - TC_raw / TC_raw
    # Relative stdev length improvement = stdev_cur - stdev_raw / stdev_raw
    # Relative identity improvement = stdev_cur - stdev_raw / stdev_raw
    # Relative trim length improvement = length_cur - length_raw / length_raw

    data_improvement = {'gene': [], 'tree_certainty_rel_imp': [], 'identity_rel_imp': [], 'stdev_prot_length_rel_imp': [], 'trimmed_length_rel_imp': []}

    for gene in genes:
        data_improvement['gene'] = gene
        data_improvement['tree_certainty_rel_imp'].append(      (float(data[(data['gene'] == gene) & (data['status'] == 'curated')]['rel_tree_certainty'])    -   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['rel_tree_certainty']))     /   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['rel_tree_certainty']))
        data_improvement['identity_rel_imp'].append(            (float(data[(data['gene'] == gene) & (data['status'] == 'curated')]['aln_identity'])          -   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['aln_identity']))           /   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['aln_identity']))
        data_improvement['stdev_prot_length_rel_imp'].append(   (float(data[(data['gene'] == gene) & (data['status'] == 'curated')]['stdev_raw_length'])      -   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['stdev_raw_length']))       /   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['stdev_raw_length']))
        data_improvement['trimmed_length_rel_imp'].append(      (float(data[(data['gene'] == gene) & (data['status'] == 'curated')]['trimmed_length'])        -   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['trimmed_length']))         /   float(data[(data['gene'] == gene) & (data['status'] == 'raw')]['trimmed_length']))

    improvement_data = pd.DataFrame(data=data_improvement)

#################################################################
    plt.rcParams["font.family"] = "Times New Roman"
    fig, axs = plt.subplots(2, 2, sharex=True)
    fig.set_size_inches(8, 6)

#################################################################
    ax1 = sns.regplot(x='tree_certainty_rel_imp', y='identity_rel_imp', data=improvement_data, ax=axs[0][0], ci=None)

    ax1.set_ylim(-1.1,1.1)
    ax1.set_xlim((-0.1,0.7))
    ax1.set_xlabel('', fontsize=12)
    ax1.set_ylabel('Identity improvement', fontsize=12)

    ax1.text(-0.4,0.8,'A', weight = 'bold', size = 30)
    ax1.tick_params(labelsize=6)

#################################################################
    ax2 = sns.regplot(x='tree_certainty_rel_imp', y='stdev_prot_length_rel_imp', data=improvement_data, ax=axs[0][1], ci=None)

    ax2.set_ylim(-1.1,1.1)
    ax2.set_xlim((-0.1,0.7))
    ax2.set_xlabel('', fontsize=12)
    ax2.set_ylabel('Protein length improvement', fontsize=12)

    ax2.text(-0.4,0.8,'B', weight = 'bold', size = 30)
    ax2.tick_params(labelsize=6)

#################################################################
    ax3 = sns.regplot(x='tree_certainty_rel_imp', y='trimmed_length_rel_imp', data=improvement_data, ax=axs[1][0], ci=None)

    ax3.set_ylim(-1.1,1.1)
    ax3.set_xlim((-0.1,0.7))
    ax3.set_xlabel('Tree certainty improvement', fontsize=12)
    ax3.set_ylabel('Alignment length improvement', fontsize=12)

    ax3.text(-0.4,0.8,'C', weight = 'bold', size = 30)
    ax3.tick_params(labelsize=6)

#################################################################
    ax4 = sns.kdeplot(improvement_data['tree_certainty_rel_imp'], ax=axs[1][1], shade=True)

    ax4.set_xlim((-0.1,0.7))
    ax4.set_xlabel('Tree certainty improvement', fontsize=12)
    ax4.set_ylabel('Frequency', fontsize=12)

    ax4.text(-0.4,4.35,'D', weight = 'bold', size = 30)
    ax4.tick_params(labelsize=6)
    ax4.get_legend().remove()

    plt.subplots_adjust(top=0.9, bottom=0.15, left=0.2, right=0.95, hspace=0.05, wspace=0.5)
    plt.savefig('gene_curation.png', dpi = 1200)
