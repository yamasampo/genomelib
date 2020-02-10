""" Classes for SFS data """

import numpy as np
import pandas as pd
from .db import Database
from .constants import FREQUENCY_POOL_CLASS, MUTATION_POOL_CLASS, MUTATIONS
from barplot.barplot import plot
from myBasic.num import bootstrap

class GeneSFS(object):
    def __init__(self, gene_ids, mut_sfs_d, description=''):
        self.gene_ids, self.mut_sfs_d = np.array(gene_ids), mut_sfs_d
        self.description = description

    @classmethod
    def from_geneSFS_path(cls, geneSFS_path, print_gene_num=False):
        """ Returns a gene id list and a dictionary containing a matrix  of SFS 
        for each gene as a list.
        
        Parameter
        ---------
        file_path: str
            A path to a gene by gene SFS file (TMg format)
        
        Returns
        -------
        gene_id_list: list
            A list of gene indices that are retrieved from discription lines 
            (starts with ">")
        mut_sfs_d: dict
            A dictionary whose keys are mutations and values are a matrix containing 
            a list of gene by gene SFS. The order of SFS appearing in the matrix
            is the same as gene_id_list.
        
        """
        tmp_sfs_d = {}
        gene_sfs = []
        gene_id_list = []
        gene_id = None

        with open(geneSFS_path, 'r') as f:
            for line in f:
                if line.startswith('/*') or line.startswith('itemnum'):
                    continue
                    
                if line.startswith('>'):
                    if gene_id:
                        gene_id_list.append(gene_id)
                    if '_' in line.rstrip()[1:]:
                        gene_id = line.rstrip()[1:]
                    else:
                        gene_id = int(line.rstrip()[1:])
                    
                else:
                    items = line.rstrip().split('\t')
        #             if items[1:] == ['0' for _ in range(len(items[1:]))]:
        #                 continue
                    
                    mutation = items[0]
                    gene_sfs = [float(d) for d in items[1:]]
                    # if mutation not in {'TC', 'CT', 'AG', 'GA'}:
                    #     continue
                    
                    try:
                        tmp_sfs_d[mutation].append(gene_sfs)
                    except KeyError:
                        tmp_sfs_d[mutation] = [gene_sfs]
                        
            if gene_id:
                gene_id_list.append(gene_id)
                
        # convert list to numpy array instance
        mut_sfs_d = dict.fromkeys(tmp_sfs_d.keys())
        for mut, sfs_matrix in tmp_sfs_d.items():
            mut_sfs_d[mut] = np.array(sfs_matrix)
            if print_gene_num:
                print('{}\tfound in {} genes'.format(mut, len(sfs_matrix)))

        return GeneSFS(gene_id_list, mut_sfs_d)

    @staticmethod
    def pool_geneSFS(geneSFS_list):
        pooled_gsfs = GeneSFS([], [])

        for gsfs in geneSFS_list:
            # gsfs = GeneSFS.from_geneSFS_path(geneSFS_path)

            if len(pooled_gsfs.gene_ids) > 0:
                pooled_gene_id = list(set(pooled_gsfs.gene_ids) & 
                                      set(gsfs.gene_ids))

                pooled_gsfs = pooled_gsfs.extract_gene_subset(pooled_gene_id)
                subset_gsfs = gsfs.extract_gene_subset(pooled_gene_id)

                assert pooled_gsfs.gene_ids.all() == subset_gsfs.gene_ids.all()

                for mutation, sfs_matrix in subset_gsfs.mut_sfs_d.items():
                    try:
                        pooled_gsfs.mut_sfs_d[mutation] += np.array(sfs_matrix)
                    except KeyError:
                        pooled_gsfs.mut_sfs_d[mutation] = np.array(sfs_matrix)

            else:
                pooled_gsfs = GeneSFS(gsfs.gene_ids, gsfs.mut_sfs_d)
                
        return pooled_gsfs

    def extract_gene_subset(self, gene_list, gene_match_func=None):
        out_d = dict.fromkeys(self.mut_sfs_d.keys())
        if not gene_match_func:
            gene_match_func = lambda gene_name: gene_name
            
        out_gene_idx = [i 
            for i, gene in enumerate(self.gene_ids) 
                if gene_match_func(gene) in gene_list]
        
        for mutation, sfs_matrix in self.mut_sfs_d.items():
            out_d[mutation] = sfs_matrix[out_gene_idx]
            
        gene_id = np.array(self.gene_ids)
        return GeneSFS(gene_id[out_gene_idx], out_d)

    def extract_mutation_subset(self, mutations, new_description=''):
        new_mut_d = {}
        for mut, sfs_matrix in self.mut_sfs_d.items():
            if mut in set(mutations):
                new_mut_d[mut] = sfs_matrix
        
        return GeneSFS(self.gene_ids, new_mut_d, new_description)

    def to_SFS(self, species, cod_type, rep_num):
        """ Bootstrap genes from a given list of genes for rep_num times. This
        returns FDboot object. This function assumes that the orders of appearance of 
        genes in gene_id_list and mut_sfs_d are identical.

        Parameters
        ----------
        gene_id_list: list-like
        mut_sfs_d: dict
        rep_num: int
        species: str
        description: str
        gene_id_outdir: str (file path)
        df_outpath: str (file path)

        """
        # get a matrix of bootstrap indices for gene_id_list
        # shape is rep_num x len(gene_id_list)
        boot_indices = bootstrap(data=np.array(self.gene_ids), 
                                rep_num=rep_num, output='index')
        # add original data to replication indices
        orig_id = np.array([range(len(self.gene_ids))])
        boot_indices2 = np.concatenate((orig_id, boot_indices))

        # initialize lists
        count_list = []
        freq_list = []
        mut_list = []
        rep_list = []
        bin_list = []
        # loop for each replication
        for rep, indices in enumerate(boot_indices2):
            for mutation, sfs_matrix in self.mut_sfs_d.items():
                boot_sfs_dat = sfs_matrix[indices] # get SFS matrix for sets of genes
                boot_pool_dat = np.sum(boot_sfs_dat, axis=0) # pool across genes
                
                for f in range(1, len(boot_pool_dat)+1):
                    count_list.append(boot_pool_dat[f-1])
                    freq_list.append(f)
                    mut_list.append(mutation)
                    rep_list.append(rep)
                    bin_list.append(1)
                    
        df = pd.DataFrame({
            'rep': rep_list,
            'mutation': mut_list,
            'frequency': freq_list,
            'count': count_list
        })
        
        return boot_indices2, SFS(df, species, rep_num, cod_type, self.description)

    def to_file(self, path):
        lines = []

        for i, gene_id in enumerate(self.gene_ids):
            lines.append(f'>{gene_id}')

            for mutation in self.mut_sfs_d.keys():
                sfs = self.mut_sfs_d[mutation][i]
                lines.append(
                    '\t'.join([mutation] + list(map(str, sfs)))
                )

        with open(path, 'w') as f:
            print('\n'.join(lines), file=f)

    def __len__(self):
        return len(self.gene_ids)

    def __repr__(self):
        return '<{name}: {desc} ({row} genes)>'.format(
            name=type(self).__name__, desc=self.description, row=self.__len__()
        )


class SFS(Database):
    """ Constructs Site Frequency Spectrum (SFS) object. """
    def __init__(self, df, species, rep_num, cod_type, description=''):        
        # assert species
        if species == 'Dm':
            freq = 14
        elif species == 'Ds':
            freq = 21
        else:
            raise Exception('Unknown species {}. Please input Dm or Ds.'\
                .format(species))
        self.species = species

        # Filter out unnecessary mutations
        if cod_type == '2f':
            self.mutations = ['AG', 'CT', 'GA', 'SW', 'TC', 'WS']
        else:
            self.mutations = MUTATIONS
        df = df[df['mutation'].isin(self.mutations)]

        # Set main attributes
        self.df = df
        self.description = description

        # assert frequency classes
        assert self.df['frequency'].max() == freq, \
            '{} frequency classes found instead of {}'.format(
                self.df['frequency'].max(), freq)
        
        # assert replications
        assert self.df['rep'].max() == rep_num, \
            '{} replications found instead of {}'.format(
                self.df['rep'].max(), rep_num)

        self.rep_num = rep_num
        # pool by mutation class
        self.reasign_df_pooled_by_mutation()
        
    @property
    def pooled_df_by_frequency(self):
        columns = ['rep', 'mutation', 'frequency', 'count']
        concat_list = []
        
        for out_freq, group_freq in FREQUENCY_POOL_CLASS[self.species].items():
            tmp_df = self.filter(frequency=group_freq)\
                          .groupby(by=['rep', 'mutation']).sum().reset_index()
            tmp_df['frequency'] = out_freq

            concat_list.append(tmp_df.loc[:, columns])

        return pd.concat(concat_list)

    @property
    def base_table(self):
        return self.get_orig_table(scale='count', sfs_type='fdd', pool_freq=False)
    
    def reasign_df_pooled_by_mutation(self):
        assert set(self.df['mutation']) != set(self.mutations), 'You have enough mutation classes'
        columns = ['rep', 'mutation', 'frequency', 'count']
        concat_list = [self.df.loc[:, columns]]
        
        for pool_mut, mut_list in MUTATION_POOL_CLASS.items():
            if pool_mut not in set(self.mutations):
                continue

            tmp_df = self.filter(mutation=mut_list)\
                          .groupby(by=['rep', 'frequency']).sum().reset_index()
            tmp_df['mutation'] = pool_mut

            concat_list.append(tmp_df.loc[:, columns])

        self.df = pd.concat(concat_list)

    def get_table(self, filter_kw, groupby_kw, func=None, add_column=None, 
                  pivot_table_kw={
                      'values': ['count'], 
                      'index': ['mutation'], 
                      'columns': ['frequency']}):
        """ Returns SFS table. This is a base function to get SFS table. You can 
        pool the SFS table (output of this function) by frequency and/or mutation 
        using another function "pool".
        """
        if not func:
            func = lambda x: x

        # Get subset of dataframe
        subset_df = self.filter(**filter_kw)
        
        # Get temporary table
        tmp_df = subset_df.groupby(**groupby_kw)\
                        .agg(func)

        # Add column just before pivot_table
        if add_column:
            # Check if only 1 data is trying to be added
            assert len(add_column) == 1, 'Only 1 additional data is allowed.'
            
            # Parse argument
            column_name = list(add_column.keys())[0]
            data = list(add_column.values())[0]

            # Add data
            tmp_df[column_name] = data
        
        return tmp_df.pivot_table(**pivot_table_kw)['count']

    def get_orig_table(self, scale, sfs_type, pool_freq):
        """ Returns SFS table for the original data (rep=0) for every mutations 
        including WS, SW, WW, SS and neu (WWSS).

        Parameter
        ---------
        scale: "count" or "prop" (default: "count")
            specifies the scaling scheme of output SFS table. "prop" is passed, 
            frequency values will become propotion of sites. 
            "prop" is not aveilable for now (190404).

        Note
        ----
        This function is 1.5 times faster than a function in the previous 
        classes.FrequencyDistributionBootstrapncy class.

        FrequencyDistributionBootstrapncy ->
            %timeit sfs.get_fdd(bin_num=1, rep_num=0)
            163 ms ± 725 µs per loop (mean ± std. dev. of 7 runs, 1 loop each) 
        SFS ->
            %timeit sfs.get_orig_table()
            90.1 ms ± 2.35 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
        """
        # Get individual mutation classes
        orig_table = self.get_table(
            func=lambda x: x, filter_kw={'rep': 0}, 
            groupby_kw={'by': ['mutation', 'frequency']})
        
        # Pool mutation classes
        if sorted(orig_table.index.tolist()) != sorted(self.mutations):
            concat_list = [orig_table]
            pooled_table = self.pool(by='mutation', table=orig_table)
            concat_list.append(pooled_table)
            orig_table = pd.concat(concat_list)
        
        # Pool frequencies
        if pool_freq:
            orig_table = self.pool(by='frequency', table=orig_table)
            
        return self.format_table(orig_table, scale, sfs_type)

    def get_mean_table(self, scale, sfs_type, pool_freq, filter_kw={'rep': 'ne0'}):
        # get base SFS table
        mean_table = self.get_table(
            filter_kw, 
            groupby_kw={'by': ['mutation', 'frequency']}, 
            func=np.mean)
        
        # Pool mutation classes
        if sorted(mean_table.index.tolist()) != sorted(self.mutations):
            concat_list = [mean_table]
            pooled_table = self.pool(by='mutation', table=mean_table)
            concat_list.append(pooled_table)
            mean_table = pd.concat(concat_list)
        
        # Pool frequencies
        if pool_freq:
            mean_table = self.pool(by='frequency', table=mean_table)
        
        return self.format_table(mean_table, scale, sfs_type)

    def get_CI_tables(self, ci, scale, sfs_type, pool_freq):
        """ Returns low CI table and high CI tables. 
        
        Parameter
        ---------
        ci: int or float
            percentage for Confidence Intervals. The value has to be within 
            0 ~ 100.
        """
        return self._CI_table(ci, 'lower', scale, sfs_type, pool_freq), \
               self._CI_table(ci, 'upper', scale, sfs_type, pool_freq)

    def _CI_table(self, ci, pos, scale, sfs_type, pool_freq):
        # Get percentile
        if pos == 'lower':
            per = (100-ci) / 100 / 2
        elif pos == 'upper':
            per = 1 - ((100-ci) / 100 / 2)
        else:
            raise Exception('pos has to be "lower" or "upper".')

        # Define function that returns a values at the position
        f =lambda rep_dat: rep_dat.quantile(per)

        # Get CIs for every individual mutation
        if set(self.df['mutation']) != set(self.mutations):
            self.reasign_df_pooled_by_mutation()
        
        # Set parameters
        if pool_freq:
            if self.species == 'Dm':
                div_freq = 8
            else:
                div_freq =  9
        else:
            if self.species == 'Dm':
                div_freq = 14
            else:
                div_freq =  21
        
        if sfs_type == 'fd':
            filter_kw = {'rep': 'ne0', 'frequency': 'ne{}'.format(div_freq)}
        else:
            filter_kw = {'rep': 'ne0'}
            
        groupby_kw = {'by': ['mutation', 'frequency']}
        pivot_table_kw={
            'values': ['count'], 'index': ['mutation'], 'columns': ['frequency']}
        
        # Get subset data
        if pool_freq:
            subset_df = Database(self.pooled_df_by_frequency).filter(**filter_kw)
        else:
            subset_df = self.filter(**filter_kw)
            
        if scale == 'prop':
            # compute proportion of site count for each mutation and rep
            subset_df = subset_df.set_index(['mutation', 'rep']) / \
                subset_df.groupby(['mutation', 'rep']).sum()
            subset_df['count'] = subset_df['count'].fillna(0)

            # format frequency values
            subset_df['frequency'] = subset_df['frequency'] * sum(range(1, div_freq))
            subset_df['frequency'] = subset_df['frequency'].astype(int)            
            
        ci_table = subset_df.groupby(by=['mutation', 'frequency'])\
                            .agg(f)\
                            .pivot_table(**pivot_table_kw)['count']
        
        return ci_table

    def get_yerr_tables(self, ci, scale, sfs_type, pool_freq, data_table):
        """ Returns two tables for distance between data and lower CI ,and 
        between data and upper CI """
        
        ciL_table, ciH_table = self.get_CI_tables(ci, scale, sfs_type, pool_freq)
        
        assert ciL_table.shape == data_table.shape, \
            'ciL_table shape {} is different from data_table shape {}.'\
            .format(ciL_table.shape, data_table.shape)
        
        assert ciH_table.shape == data_table.shape, \
            'ciH_table shape {} is different from data_table shape {}.'\
            .format(ciH_table.shape, data_table.shape)

        return abs(data_table - ciL_table), abs(data_table - ciH_table)

    def pool(self, by, table):
        """ Pool SFS according to a given mode 
        Parameter
        ---------
        by: str
            specify based on what SFS will be pooled.
            This has to be "frequency" or "mutation".
        """
        if by == 'frequency':
            new_sfs = self._pool_by_freq(table, self.species)
        
        elif by == 'mutation':
            new_sfs = self._pool_by_mutation(table)
            
        else:
            raise Exception('by has to be "frequency" or "mutation".')
        
        return new_sfs
    
    @staticmethod
    def _pool_by_freq(table, species):
        """ Pool frequency of segregating sites. """
        concat_list = []

        for out_freq, in_freqs in FREQUENCY_POOL_CLASS[species].items():
            try:
                tmp_s = pd.Series(
                    table[in_freqs].apply('sum', axis=1), 
                    name=out_freq)

            except KeyError:
                continue

            concat_list.append(
                tmp_s
            )

        return pd.concat(concat_list, axis=1)
    
    @staticmethod
    def _pool_by_mutation(table):
        # Pool mutation classes
        concat_list = []
        mutations = table.index
        
        for pool_mut, mut_list in MUTATION_POOL_CLASS.items():
            retain_mutations = list(set(mut_list) & set(mutations))

            try:
                tmp_s = pd.Series(
                    table.loc[retain_mutations].apply('sum', axis=0), 
                    name=pool_mut)

            except KeyError:
                continue

            concat_list.append(
                tmp_s
            )
        
        return pd.concat(concat_list, axis=1).T
    
    @staticmethod
    def format_table(table, scale, sfs_type):
        # TODO: Add assertion for SFS table format

        # Return table as specified format
        if scale == 'count':
            if sfs_type == 'fd':
                return table.iloc[:, :-1]
            elif sfs_type == 'fdd':
                return table
            else:
                raise Exception('sfs_type has to be "fd" or "fdd".')

        elif scale == 'prop':
            if sfs_type == 'fd':
                return table.iloc[:, :-1]\
                             .apply(lambda x: x / x.sum(), axis=1)
            elif sfs_type == 'fdd':
                return table.apply(lambda x: x / x.sum(), axis=1)
            else:
                raise Exception('sfs_type has to be "fd" or "fdd".')

        else:
            raise Exception('scale has to be "count" or "prop".')

    def plot(self, ax, mutation1, mutation2, scale, sfs_type, 
             pool_freq, data_type='orig', ci=95, test='', style='default', **plot_kw):

        assert mutation1 in self.mutations, '{} not in SFS dataset {}'\
            .format(mutation1, self.description)
        assert mutation2 in self.mutations, '{} not in SFS dataset {}'\
            .format(mutation2, self.description)

        if data_type == 'orig':
            data_table = self.get_orig_table(scale, sfs_type, pool_freq)
        elif data_type == 'mean':
            data_table = self.get_mean_table(
                scale, sfs_type, pool_freq, filter_kw={'rep': 'ne0'})
        else:
            raise Exception('Input your data_table or data_type as either '\
                            'of "orig" or "mean"')
        
        yerrL, yerrH = self.get_yerr_tables(ci, scale, sfs_type, pool_freq, data_table)
        
        data_mtx = [
            data_table.loc[mutation1].values,
            data_table.loc[mutation2].values
        ]
        data1_yerrs = [
            yerrL.loc[mutation1].values,
            yerrH.loc[mutation1].values
        ]
        data2_yerrs = [
            yerrL.loc[mutation2].values,
            yerrH.loc[mutation2].values
        ]
        
        if style == 'default':
            plot_data = plot(ax, data_mtx, yerrs=[data1_yerrs, data2_yerrs], 
                 with_yerrs=True, xticklabels=data_table.columns, 
                 data_series=[mutation1, mutation2],
                 color_list=['gray', 'yellow'], return_data=True, ci=ci,
                 edgecolor='black', linewidth=0.4, alpha=0.9,
                 error_kw={'elinewidth': 1}, space=0.1, **plot_kw)
            
            ax.patch.set_facecolor('#f6f5ea')
            ax.patch.set_alpha(0.5)
        
        if scale == 'count':
            ax.set_ylabel('Number of sites')
        elif scale == 'prop':
            ax.set_ylabel('Proportion of sites')
            
        ax.set_xlabel('Mutation frequency')
        
        return plot_data
    