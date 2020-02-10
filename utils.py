import os
import numpy as np

def to_filelist_file(dir_path, file_name='0.filelist'):
    '''Return a file list in a given directory'''
    flist= [a for a in os.listdir(dir_path) if 
        not a.startswith('.') and not a.startswith('0')]

    with open(os.path.join(dir_path, file_name), 'w') as f:
        print('itemnum: '+str(len(flist)), file=f)
        print('\n'.join(flist), file=f)

    return flist

def to_genelist_file(genelist, path):
    with open(path, 'w') as f:
        print('itemnum: {}'.format(len(genelist)), file=f)
        print('\n'.join(list(map(str, genelist))), file=f)

def genelists_to_genelist_files(out_dir, base_genelist, rep_indices):
    base_genelist  = np.array(base_genelist)

    for rep, ix in enumerate(rep_indices):
        genelist = base_genelist[ix]
        path = os.path.join(out_dir, f'gene_id_rep_{rep}.txt')

        to_genelist_file(genelist, path)

    flist = to_filelist_file(out_dir)
    return len(flist)
