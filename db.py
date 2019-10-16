
import re
import pandas as pd
from collections.abc import Mapping

class Database(Mapping):
    """ This class inherits Mapping class. __iter__, __getitem__ and __len__ 
    functions are overwritten. This is a base class of SFS class. """
    def __init__(self, df, description=''):
        self.df = df
        self.description = description
        
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
                elif re.search('^gt\d+', v):
                    v = float(re.search('^gt(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] > v]
                elif re.search('^gte\d+', v):
                    v = float(re.search('^gte(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] >= v]
                elif re.search('^lt\d+', v):
                    v = float(re.search('^lt(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] < v]
                elif re.search('^lte\d+', v):
                    v = float(re.search('lte(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] <= v]
                elif re.search('^ne\d+', v):
                    v = float(re.search('ne(\d+\.*\d*)$', v).group(1))
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
        return '<{name}: {desc} ({size} records)>'.format(
            name=type(self).__name__, desc=self.description, size=self.__len__()
        )

class EvoGenDatabase(Database):
    def __init__(
        self, 
        df:pd.DataFrame, # Gene info table.
        seq_d:pd.DataFrame, # Dictionary or DataFrame of DNA sequence 
        sfs_d:dict={}, # Gene SFS
        description=''): # Description on this dataset
        """ 
        Attributes
        ----------
        df: pd.DataFrame
            A row is a gene and you can list any kind of information in columns
        """
        self.df = self._merge_seq_and_metadata(df, seq_d)
        if sfs:
            self.load_sfs(sfs)
        self.description = description

    @property
    def reduced_df(self):
        # Returns only columns except SFS data
        # columns starting with "sfs" indicate SFS data columns
        cols = [col for col in self.df.columns if not col.startswith('sfs')]
        return self._df[cols]

    @classmethod
    def from_FASTA(
        cls, df_path, seq_path, seq_fmt, sfs_path='', description='', **kwargs):
        df = pd.read_csv(df_path, **kwargs)
        seq_d = read_seq_file(seq_path, seq_fmt)
        sfs = parse_gene_sfs(sfs_path) if sfs_path else {}

        return cls(df, seq_d, sfs, description)

    @staticmethod
    def read_seq_file(seq_path:str, seq_fmt:str, **kwargs):
        if seq_fmt == 'fasta':
            return _read_fasta_to_dict(seq_path)
        if seq_fmt == 'pickle':
            return _unpickle_seq_d(seq_path)
        if seq_fmt == 'df':
            return pd.read_csv(seq_path, **kwargs)
        raise IndexError(f'Unknown format argument {seq_fmt} found. '\
            'Currently "fasta", "pickle" or "df" are only supported.')

    @staticmethod
    def _read_fasta_to_dict(fasta_path):
        pass

    @staticmethod
    def _unpickle_seq_d(pickle_path):
        with open(pickle_path, 'rb') as f:
            seq_d = pickle.load(f)

        return seq_d

    @staticmethod
    def parse_gene_sfs(sfs_path):
        pass

    @staticmethod
    def _merge_seq_and_metadata(df, seq_d):
        pass

    def load_sfs(self, sfs):
        self._sfs_d = sfs
        # Conncet to df table
        # TODO: Think how store SFS metadata such as codon and amino acid types?
        # Do we need all SFS data or pooled across genes or across mutations,
        # or all sfs data has to be separated?
        # Maybe store SFS object in a dictionary.

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

