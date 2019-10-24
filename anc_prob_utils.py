
import os
import pandas as pd
from collections import OrderedDict
from .constants import AA_CODON_3rdpos_BASES

def get_anc_base_comp_df(
    collapse_seq_flist, anc_prob_file, node_name1, node_name2, 
    aminoacid='', collapse_seq_dir=''):

    # Add gene_id column to ancestor probability table
    if not collapse_seq_dir:
        collapse_seq_dir = os.path.dirname(collapse_seq_flist)
    seq_len_d = get_seq_len_d(collapse_seq_flist, collapse_seq_dir)
    anc_prob_table = read_anc_prob_table(anc_prob_file)
    anc_prob_table = add_gene_id(anc_prob_table, seq_len_d)

    # Count base composition
    anc_prob_df = melt_anc_prob_table(anc_prob_table, node_name1, node_name2)
    
    if aminoacid:
        bases = AA_CODON_3rdpos_BASES[aminoacid]
    else:
        bases = list('TCAG')

    anc_base_comp_df1 = count_ancestor_base_comp(
        anc_prob_df, node_name1, bases)
    anc_base_comp_df2 = count_ancestor_base_comp(
        anc_prob_df, node_name2, bases)

    return anc_base_comp_df1, anc_base_comp_df2

def get_seq_len(aln_file):
    tmp_seq_list = []
    seq_len = None
    
    with open(aln_file, 'r') as f:
        for l in f:
            if l.startswith('>'):
                if tmp_seq_list:
                    seq = ''.join(tmp_seq_list)
                    
                    if seq_len:
                        assert seq_len == len(seq), \
                            f'{seq_len} != {len(seq)} in {aln_file}'
                    else:
                        seq_len = len(seq)
                        
                    tmp_seq_list = []
            else:
                tmp_seq_list.append(l.rstrip())
                
    if tmp_seq_list:
        seq = ''.join(tmp_seq_list)
        if seq_len:
            assert seq_len == len(seq)
        else:
            seq_len = len(seq)
            
    tmp_seq_list = []
            
    return seq_len

def get_seq_len_d(collapse_seq_flist, collapse_seq_dir=''):
    with open(collapse_seq_flist, 'r') as f:
        flist = [
            fname.rstrip() for fname in f if not fname.startswith('itemnum')
        ]
        
    gene_id_parser = lambda fname: int(fname.split('.')[0].split('_')[1])

    collapse_seq_dir = os.path.dirname(collapse_seq_flist)

    return OrderedDict(
        {
            gene_id_parser(fname): 
                get_seq_len(
                    os.path.join(collapse_seq_dir, fname)) for fname in flist
        }
    )

def read_anc_prob_table(anc_prob_file):
    pos_list = []
    freq1_list = []
    freq2_list = []
    config_list = []
    anc_probs_list = []

    with open(anc_prob_file, 'r') as f:
        for l in f:
            parts = l.rstrip().split(':\t\t')
            pos, freq1, freq2, config = parts[0].split('\t')
            
            pos_list.append(int(pos))
            freq1_list.append(int(freq1))
            freq2_list.append(int(freq2))
            config_list.append(config)
            anc_probs_list.append(parts[1])
            
    return pd.DataFrame(
        {
            'pos': pos_list,
            'freq1': freq1_list,
            'freq2': freq2_list,
            'config': config_list,
            'anc_probs': anc_probs_list
        }
    )
    
def lookup_geneinfo(pos, seq_len_d):
    total_pos = 0
    
    for geneinfo, seq_len in seq_len_d.items():
        total_pos += seq_len
        
        if pos <= total_pos:
            return geneinfo
        
    raise f'Not found: {pos} > {total_pos}'

def add_gene_id(anc_table, seq_len_d):
    anc_table['gene_id'] = anc_table.apply(
        lambda x: lookup_geneinfo(x['pos'], seq_len_d),
        axis=1
    )

    return anc_table

def list_to_DF(series):
    parts = series['anc_probs'].split('\t')
    anc_probs = pd.DataFrame(
        [parts[i:i+3] for i in range(0, len(parts), 3)],
        columns=['ms', 'm', 'prob']
    )
    anc_probs['gene_id'] = series['gene_id']
    anc_probs['pos'] = series['pos']
    anc_probs['freq1'] = series['freq1']
    anc_probs['freq2'] = series['freq2']
    anc_probs['config'] = series['config']
    
    return anc_probs

def melt_anc_prob_table(anc_table, node_name1, node_name2):
    anc_df = pd.concat(
        list(anc_table.apply(list_to_DF, axis=1)), 
        axis=0
    )
    anc_df['prob'] = anc_df['prob'].astype(float)
    anc_df = anc_df.loc[:, 
        ['pos', 'gene_id', 'freq1', 'freq2', 'config', 
         node_name1, node_name2, 'prob']]

    return anc_df

def count_ancestor_base_comp(anc_df, node_name, bases=['T', 'C', 'A', 'G']):
    anc_base_comp_df = anc_df\
        .loc[:, ['gene_id', node_name, 'prob']]\
        .groupby(['gene_id', node_name])\
        .sum()\
        .pivot_table(values='prob', index='gene_id', columns=node_name)

    if not bases:
        bases = ['T', 'C', 'A', 'G']
    # Filter bases columns
    anc_base_comp_df = anc_base_comp_df.loc[:,bases]

    anc_base_comp_df.fillna(0, inplace=True)

    # Add total site column
    anc_base_comp_df['total'] = anc_base_comp_df.loc[:, bases].sum(axis=1)

    # Add GC content column
    gc_col = set(bases) & set(['G', 'C'])
    anc_base_comp_df['GC_cont'] = \
        anc_base_comp_df.loc[:, bases].sum(axis=1) / anc_base_comp_df['total']

    return anc_base_comp_df
