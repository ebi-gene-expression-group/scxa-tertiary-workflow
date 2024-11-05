
import scanpy as sc
import anndata
from numpy import all
import logging
import os

adata = sc.read('input.h5')


gene_name = 'index'
qc_vars = list()



gene_names = getattr(adata.var, gene_name)
# Define the directory containing your source files
source_dir = '.'  # Adjust to the appropriate path

ad_s = sc.read('r_source.h5')
if not all(adata.obs.index.isin(ad_s.obs.index)):
  logging.error("Specified object for .raw must contain all .obs from main object.")
  sys.exit(1)
else:
  adata.raw = ad_s[adata.obs.index]
del ad_s

ad_s = sc.read('x_source_0.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  if "filtered" == '':
    logging.error("%sth destination layer for %sth X source not specified" % ("0", "0"))
    sys.exit(1)
  adata.layers["filtered"] = ad_s.X
else:
  logging.error("X source 0 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('x_source_1.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  if "normalised" == '':
    logging.error("%sth destination layer for %sth X source not specified" % ("1", "1"))
    sys.exit(1)
  adata.layers["normalised"] = ad_s.X
else:
  logging.error("X source 1 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s


ad_s = sc.read('obs_source_0.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_0"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 0 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_1.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_1"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 1 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_2.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_2"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 2 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_3.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_3"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 3 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_4.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_4"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 4 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_5.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_5"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 5 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_6.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_6"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 6 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_7.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_7"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 7 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('obs_source_8.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obs.keys() if "louvain" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.obs:
        suffix = "_8"

    adata.obs[[k_to_copy+suffix]] = ad_s.obs[[k_to_copy]]
    if k_to_copy in ad_s.uns.keys():
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Observation source 8 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s




embedding_sources = [file for file in os.listdir(source_dir) if file.startswith('embedding_source_') and file.endswith('.h5')]

for idx, embedding_file in enumerate(sorted(embedding_sources)):
    ad_s = sc.read(os.path.join(source_dir, embedding_file))
    if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
          # Copy tsne embeddings
          keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
          for k_to_copy in keys_to_copy:
              suffix = ''
              if k_to_copy in adata.obsm:
                  suffix = f"_{idx}"
              adata.obsm[k_to_copy + suffix] = ad_s.obsm[k_to_copy]
          # Copy umap embeddings
          keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
          for k_to_copy in keys_to_copy:
              suffix = ''
              if k_to_copy in adata.obsm
                  suffix = f"_{idx}"
              adata.obsm[k_to_copy + suffix] = ad_s.obsm[k_to_copy]
    else:
        logging.error(f"Embedding source {idx} AnnData file is not compatible to be merged to main AnnData file, different cell names.")
        sys.exit(1)
    del ad_s


uns_sources = [file for file in os.listdir(source_dir) if file.startswith('uns_source_') and file.endswith('.h5')]
for idx, uns_file in enumerate(sorted(uns_sources)):
  ad_s = sc.read(os.path.join(source_dir, uns_file))
  if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
    keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
    for k_to_copy in keys_to_copy:
      suffix=''
      if k_to_copy in adata.uns:
          suffix=f"_{idx}"
      adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
  else:
    logging.error(f"Uns source {idx} AnnData file is not compatible to be merged to main AnnData file, different cell names.")
    sys.exit(1)
  del ad_s


if len(qc_vars) > 0:
    pct_top = [50]
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=pct_top, inplace=True)

if 'n_genes' not in adata.obs.columns:
    sc.pp.filter_cells(adata, min_genes=0)
if 'n_counts' not in adata.obs.columns:
    sc.pp.filter_cells(adata, min_counts=0)
if 'n_cells' not in adata.var.columns:
    sc.pp.filter_genes(adata, min_cells=0)
if 'n_counts' not in adata.var.columns:
    sc.pp.filter_genes(adata, min_counts=0)

adata.write('output.h5', compression='gzip')
    
