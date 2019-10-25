
import re
import pickle
import pandas as pd
from collections.abc import Mapping
from .constants import *

class Database(Mapping):
    """ This class inherits Mapping class. __iter__, __getitem__ and __len__ 
    functions are overwritten. This is a base class of SFS class. """
    def __init__(self, df, description=''):
        self.df = df
        self.description = description

    @property
    def columns(self):
        return self.df.columns
        
    def filter(self, sort_by='', ascending=True, **kwargs):
        '''
        Search rows which have specifies items from a given dataframe.
        Please pass key words for searching to **kwargs.
        For example, if you want to get items that is greater than equal (>=)
        100 in column "A", please specify **kwargs as "A=gte100". 
        Please see below for details.
        If nothing passed to **kwargs, return input dataframe.
        Paramters
        ---------
            df: DataFrame (pandas)
                input dataframe
            **kwargs:
                key is for column, value is for filtering values (items)
                You can use indicators below for filtering way.
                "gt" for ">"
                "gte" for ">="
                "lt" for "<"
                "lte" for "<="
                "ne" for "!="
                "c/" for "contains"
                "" for "=="
                If you pass tuple to value, this function search and filter 
                items recursively.
        Dependencies
        ------------
            pandas
            re
        '''
        res_df = self.df
        def f(res_df, k, v):
            if isinstance(v, str):
                if v == '*':
                    pass
                elif re.search('^gt-*\d+', v):
                    v = float(re.search('^gt(-{0,1}\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] > v]
                elif re.search('^gte-*\d+', v):
                    v = float(re.search('^gte(-{0,1}\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] >= v]
                elif re.search('^lt-*\d+', v):
                    v = float(re.search('^lt(-{0,1}\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] < v]
                elif re.search('^lte-*\d+', v):
                    v = float(re.search('lte(-{0,1}\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] <= v]
                elif re.search('^ne-*\d+', v):
                    v = float(re.search('ne(-{0,1}\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] != v]
                elif re.search('^c\/', v):
                    v = re.search('^c\/(.+)\/$', v).group(1)
                    res_df = res_df[res_df[k].str.contains(v)]
                elif re.search('^nc\/', v):
                    v = re.search('^nc\/(.+)\/$', v).group(1)
                    res_df = res_df[~res_df[k].str.contains(v)]
                else:
                    res_df = res_df[res_df[k] == v]
            elif isinstance(v, list):
                res_df = res_df[res_df[k].isin(v)]
            else:
                res_df = res_df[res_df[k] == v]

            return res_df

        for k,v in kwargs.items():
            if isinstance(v, tuple): # "and"
                res_df_list = []
                for i in v:
                    tmp_res_df = f(res_df, k, i)
                    res_df = pd.merge(res_df, tmp_res_df, how='inner')
            elif isinstance(v, (str, list)): # list means "or" condition
                res_df = f(res_df, k, v)
            elif isinstance(v, int):
                res_df = res_df[res_df[k] == v]
        if sort_by:
            res_df.sort_values(by=sort_by, ascending=ascending, inplace=True)
            
        return res_df
    
    def groupby(self, **kwargs):
        return self.df.groupby(**kwargs)

    def head(self, *kwargs):
        return self.df.head(*kwargs)

    def tail(self, *kwargs):
        return self.df.tail(*kwargs)

    def __len__(self):
        return len(self.df.index)
    
    def __iter__(self):
        return self.df.iterrows()
    
    def __getitem__(self, key):
        if key == '*':
            return self.df
        else:
            return self.df.loc[key, :]
        # elif isinstance(key, (tuple, list)):
        #     for k in key:
        #         yield self.df.loc[k, :]
    
    def __repr__(self):
        return '<{name}: {desc} ({row} records x {col} columns)>'.format(
            name=type(self).__name__, desc=self.description, row=self.__len__(),
            col=len(self.columns)
        )

class EvoGenDatabase(Database):    
    def __init__(
        self, 
        annot:pd.DataFrame, # Gene info table.
        seq_d:dict={}, # Dictionary or DataFrame of DNA sequence
        anc_state_df:pd.DataFrame=None, # A list of base composition for anocestor node
        sfs_d:dict={}, # Gene SFS
        repl_Nan_with=-9,
        description='', **kwagrs): # Description on this dataset
        """ 
        Attributes
        ----------
        df: pd.DataFrame
            A row is a gene and you can list any kind of information in columns
        """
        # Load main DataFrame
        if seq_d:
            self.df = self._add_seq_to_annot(annot, seq_d)
            self.df.fillna(repl_Nan_with, inplace=True)

        else:
            self.df = annot
            
        if anc_state_df:
            try:
                assert kwargs['inplace'] == True, 'inplace argument has to be '\
                    'False since class does inplace automatically.'
            except KeyError:
                pass
            self.df = self.merge_df(anc_state_df)

        if sfs_d:
            self.load_sfs(sfs_d)
        self.description = description

    # ----- Property ----- #
    @property
    def reduced_df(self):
        # Returns only columns except SFS data
        # columns starting with "sfs" indicate SFS data columns
        cols = [
            'info_id', 'chr', 'tss', 'cds_start_tssd', 'cds_end_tssd', 'cds_len', 
            'tr_len', 'strand', 'alt_spl', 'exon_num', 'gene_name', 'fbgn', 'fbtr',
            'error_code'
        ]
        return self.df[cols]

    # ----- Reading methods ----- #
    @classmethod
    def from_files(
        cls, annot_path, seq_path='', seq_fmt='', sfs_path='', repl_Nan_with=-9, 
        description='', **kwargs):
        annot = pd.read_csv(annot_path, **kwargs)
        seq_d = cls.read_seq_file(seq_path, seq_fmt)
        sfs = parse_gene_sfs(sfs_path) if sfs_path else {}

        return cls(annot, seq_d, sfs, repl_Nan_with, description)

    @classmethod
    def read_seq_file(cls, seq_path:str, seq_fmt:str, **kwargs):
        if seq_fmt == 'fasta':
            return cls._read_fasta_to_dict(seq_path)
        if seq_fmt == 'pickle':
            return cls._unpickle_seq_d(seq_path)
        if seq_fmt == 'df':
            return pd.read_csv(seq_path, **kwargs).todict()
        raise IndexError(f'Unknown format argument {seq_fmt} found. '\
            'Currently "fasta", "pickle" or "df" are only supported.')

    @staticmethod
    def _read_fasta_to_dict(fasta_path):
        seqnames = []
        seqs = []
        tmp_seq = []

        with open(fasta_path, 'r') as f:
            for l in f:
                if l.startswith('>'):
                    if tmp_seq:
                        seqnames.append(seqname)
                        seqs.append(''.join(tmp_seq))
                    
                    seqname = l.rstrip()[1:]
                    tmp_seq = []
                else:
                    tmp_seq.append(l.rstrip())

        if tmp_seq:
            seqnames.append(seqname)
            seqs.append(''.join(tmp_seq))

        return dict(zip(seqnames, seqs))

    @staticmethod
    def _unpickle_seq_d(pickle_path):
        with open(pickle_path, 'rb') as f:
            seq_d = pickle.load(f)

        return seq_d

    @staticmethod
    def parse_gene_sfs(sfs_path):
        pass

    @staticmethod
    def _add_seq_to_annot(annot, seq_d, seq_key_col='seq_id'):
        if 'seq' in annot.columns:
            raise KeyError('Column named "seq" alread exists in annot table.')

        # Put nucleotide sequence in df
        annot['seq'] = annot.apply(lambda x: seq_d[x[seq_key_col]], axis=1)
        # Check if seq length and registered len matches
        assert annot.apply(lambda x: x['cds_len'] == len(x['seq']), axis=1).all()

        return annot

    @staticmethod
    def _merge_anc_state_to_annot(annot, anc_state_df, **kwargs):
        return pd.merge(annot, anc_state_df, **kwargs)

    def merge_df(self, df2, **kwargs):
        return pd.merge(self.df, df2, **kwargs)

    def load_sfs(self, sfs_d):
        self._sfs_d = sfs
        # Conncet to df table
        # TODO: Think how store SFS metadata such as codon and amino acid types?
        # Do we need all SFS data or pooled across genes or across mutations,
        # or all sfs data has to be separated?
        # Maybe store SFS object in a dictionary.

    # ----- Add ----- #
    def apply(self, func, **kwargs):
        return self.df.apply(func, axis=1, **kwargs)

    def add_col(self, func, col_name, **kwargs):
        self.df[col_name] = self.apply(func, **kwargs)

    def add_GC3_AAs(self, cod_type):
        col_name = f'GC3_{cod_type}'
        gc_content = lambda seq: -9 if len(seq) == 0 else \
            (seq.count('G') + seq.count('C')) / len(seq)
        cds_to_cod_list = lambda cds: [cds[n:n+3] for n in range(0, len(cds), 3)]
        match_AA = lambda codon: GENETIC_CODE_i[codon] in CODON_TYPES[cod_type]
        
        func = lambda x: -9 if 'N' in x['seq'] or x['cds_len'] % 3 != 0 else \
            gc_content(
                ''.join([codon[2] 
                for codon in list(filter(match_AA, cds_to_cod_list(x['seq'])))])
            )

        self.add_col(func, col_name)

    def add_GC3_iAA(self, aa):
        col_name = f'GC3_{aa}'
        gc_content = lambda seq: -9 if len(seq) == 0 else \
            (seq.count('G') + seq.count('C')) / len(seq)
        cds_to_cod_list = lambda cds: [cds[n:n+3] for n in range(0, len(cds), 3)]
        match_AA = lambda codon: GENETIC_CODE_i[codon] == aa
        
        func = lambda x: -9 if 'N' in x['seq'] or x['cds_len'] % 3 != 0 else \
            gc_content(
                ''.join([codon[2] 
                for codon in list(filter(match_AA, cds_to_cod_list(x['seq'])))])
            )

        self.add_col(func, col_name)

    # ----- Iteration method ----- #
    def gen_seq(self, sort_by='', ascending=True, **kwargs):
        """ Generate registered nucleotide sequences. kwargs accepts options 
        for filtering CDS. 
        
        Returns
        -------
        i: int
            index of DataFrame
        seq_id: int
            key of sequence dictionary
        sequence: str
            nucleotide sequence registered in a dictionary
        """
        res_df = self.filter(sort_by, ascending, **kwargs)
        seq_id_list = res_df['seq_id'].tolist()
        id_list = list(res_df.index)
        
        for i, seq_id in zip(id_list, seq_id_list):
            yield i, seq_id, self._d[seq_id]

    # ----- Output method ----- #
    def to_fasta(fasta_path, seq_name_encoder=None, itemnum=False, 
                 sort_by='', ascending=True, **kwargs):
        if not seq_name_encoder:
            seq_name_encoder = lambda x: str(x['info_id'])

        fasta_lines = []
        for _, item in enumerate(self.gen_seq(sort_by, ascending, **kwargs)):
            i, seq_id, seq = item
            seq_name = seq_name_encoder(self[i])
            fasta_lines.append('>{}\n{}'.format(seq_name, seq))

        with open(fasta_path, 'w') as f:
            if itemnum:
                print('itemnum: {}'.format(len(fasta_lines)), file=f)
            print('\n'.join(fasta_lines), file=f)


