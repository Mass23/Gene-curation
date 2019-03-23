import glob
import subprocess
import statistics
from Bio import SeqIO

##### Functions of the pipeline
def Mafft(raw_fasta, gene_name, status):

    raw_length_list = []
    for i in SeqIO.parse(gene_name + '_' + status + '_trim.fasta', 'fasta'):
        raw_length_list.append(len(i.seq))
    stdev_length = statistics.stdev(raw_length_list)

    mafft_cmd = ['mafft', '--maxiterate', '1000', '--globalpair', raw_fasta, '>', gene_name + '_' + status + '_aln.fasta']

    try:
        subprocess.call(' '.join(mafft_cmd), shell = True)
        return(stdev_length)
    except:
        print('Mafft alignment:', gene_name, status, 'failed!')
        return('NA')


def Trimal(gene_name, status):
    trimal_cmd = ['trimal', '-automated1', '-in', gene_name + '_' + status + '_aln.fasta', '-out', gene_name + '_' + status + '_trim.fasta']

    try:
        subprocess.call(' '.join(mafft_cmd), shell = True)

        trim_length_list = []
        aln = []
        for i in SeqIO.parse(gene_name + '_' + status + '_trim.fasta', 'fasta'):
            trim_length_list.append(len(i.seq))
            aln.append(i.seq)
        trimmed_length = min(trim_length_list)

        match = 0
        for i in range(0, trimmed_length):
            current_amino_acid_list = [seq[i] for seq in aln]
            if len(set(current_amino_acid_list)) == 1:
                match += 1
        try:
            identity = match / trimmed_length
        except:
            identity = 'NA'

        return(trimmed_length, identity)

    except:
        print('Trimal trimming:', gene_name, status, 'failed!')
        return('NA')


def Raxml(gene_name, status):
    cmd_tree = ['raxml', '-p', '12345', '-m', 'PROTGAMMAWAG', '-#', '100', '-s', gene_name + '_' + status + '_trim.fasta', '-f', 'a', '-x', '12345',
           '-n', + gene_name + '_' + status, '-o', 'Drosophila_melanogaster']
    try:
        subprocess.call(' '.join(cmd_tree), shell = True)
    except:
        print('RAxML tree:', gene_name, status, 'failed!')

    cmd_tc = ['raxml', '-b', '12345', '-m', 'PROTGAMMAWAG', '-#', '100', '-f', 'i',
           '-n', gene_name + '_TC_' + status, '-z', 'RAxML_bootstrap.' + gene_name + '_' + status, '-t', 'RAxML_BestTree.' + gene_name + '_' + status, '-L', 'MR']
    try:
        subprocess.call(' '.join(cmd_tc), shell = True)
    except:
        print('RAxML tree certainty:', gene_name, status, 'failed!')

    with open('RAxML_info.' + gene_name + '_TC_' + status) as tc_file:
        for line in tc_file.readlines():
            if line.startswith('Relative tree certainty for this tree:'):
                rel_tree_certainty = line.split(' ')[-1]
    return(rel_tree_certainty)

def main():
    with open('gene_curation.csv', 'w') as out:
        out.write('gene\tstatus\tfasta_file\trel_tree_certainty\ttrimmed_length\taln_identity\tstdev_raw_length\n')

        genes_list = glob.glob('*')
        print(genes_list)

        for gene in genes_list:
            name = str(gene)

            raw_file = glob.glob(gene + '/*RAW*')
            curated_file = glob.glob(gene + '/*CURATED*')

            stdev_prot_raw = Mafft(raw_file, name, 'raw')
            stdev_prot_curated = Mafft(curated_file, name, 'curated')

            trimmed_length_raw = Trimal(name, 'raw')[0]
            trimmed_length_curated = Trimal(name, 'curated')[0]

            identity_raw = Trimal(name, 'raw')[1]
            identity_curated = Trimal(name, 'curated')[1]

            tree_certainty_raw = Raxml(name, 'raw')
            tree_certainty_curated = Raxml(name, 'curated')

            out.write('\t'.join[name, 'raw', gene, tree_certainty_raw, trimmed_length_raw, identity_raw, stdev_prot_raw] + '\n')
            out.write('\t'.join[name, 'curated', gene, tree_certainty_curated, trimmed_length_curated, identity_curated, stdev_prot_curated] + '\n')

if __name__== '__main__':
    main()
