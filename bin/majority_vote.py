#!/usr/bin/env python
import scanpy as sc
import sys
import pandas as pd
import numpy as np

methods = sys.argv[1].split(',')
adata_files = sys.argv[2].split(',')
threshold = float(sys.argv[3])

assert threshold >= 0, 'threshold must be non-negative'
assert len(methods) == len(adata_files), 'Something went wrong, the number of methods do not align with number of anndata files'

result = sc.read_h5ad(adata_files[0])
result.obs['predicted_doublet']
adata_list = []

# read in the different adata_files from different doublet methods
for m, f in zip(methods, adata_files):
    adata = sc.read_h5ad(f)
    assert np.all(result.obs.index == result.obs.index) and np.all(result.var.index == result.var.index), f'anndata from doublet {m} do not match'
    result.obs[m+'_predicted_doublet'] = adata.obs['predicted_doublet'] 
    if 'doublet_score' in result.obs:
        result.obs[m+'_doublet_score'] = adata.obs['doublet_score'] 

# columns of doublet predictions
filter_cols = result.obs.columns[result.obs.columns.str.contains(pat = '_predicted_doublet')]


# can also change the rest to filter by average doublet score, though not all methods may have this 
# filter_cols = result.obs.columns[result.obs.columns.str.contains(pat = '_doublet_score')]

#simple majority voting based on threshold, if threshold is >= 1, then filtering is based on number of methods that marks douplet
def simple_vote(meta, cols, threshold) -> np.array:
    
    votes = meta[cols].sum(axis = 1) 
    vote_ratio = votes / len(cols)
    new_predict = np.zeros(meta.shape[0], dtype = 'bool')
    if threshold > 0.99:
        new_predict [votes > threshold] = True
    else:
        new_predict [vote_ratio > threshold] = True
    return new_predict

result.obs['predicted_doublet'] = simple_vote(result.obs, filter_cols, threshold)


result.write_h5ad("doublet_major.h5ad")
