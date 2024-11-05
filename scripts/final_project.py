
import scanpy as sc
import anndata
from numpy import all
import logging

adata = sc.read('input.h5')


gene_name = 'index'
qc_vars = list()



gene_names = getattr(adata.var, gene_name)


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


ad_s = sc.read('embedding_source_0.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_0"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_0"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 0 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_1.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_1"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_1"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 1 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_2.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_2"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_2"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 2 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_3.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_3"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_3"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 3 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_4.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_4"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_4"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 4 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_5.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_5"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_5"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 5 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_6.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_6"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_6"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 6 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_7.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_7"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_7"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 7 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_8.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_8"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_8"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 8 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_9.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_9"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_9"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 9 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_10.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_10"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_10"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 10 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_11.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_11"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_11"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 11 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_12.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_12"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_12"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 12 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_13.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_13"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_13"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 13 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_14.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_14"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_14"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 14 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_15.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_15"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_15"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 15 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('embedding_source_16.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.obsm.keys() if "tsne" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_16"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
  keys_to_copy = (k for k in ad_s.obsm.keys() if "umap" in k)
  for k_to_copy in keys_to_copy:
    suffix = ''
    if k_to_copy in adata.obsm:
        suffix = "_16"
    adata.obsm[k_to_copy+suffix] = ad_s.obsm[k_to_copy]
else:
  logging.error("Embedding source 16 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s

ad_s = sc.read('uns_source_0.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_0"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 0 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_1.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_1"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 1 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_2.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_2"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 2 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_3.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_3"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 3 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_4.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_4"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 4 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_5.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_5"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 5 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_6.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_6"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 6 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_7.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_7"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 7 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
  sys.exit(1)
del ad_s
ad_s = sc.read('uns_source_8.h5')
if adata.n_obs == ad_s.n_obs and all(adata.obs_names == ad_s.obs_names):
  keys_to_copy = (k for k in ad_s.uns.keys() if "marker" in k)
  for k_to_copy in keys_to_copy:
    suffix=''
    if k_to_copy in adata.uns:
        suffix="_8"
    adata.uns[k_to_copy+suffix] = ad_s.uns[k_to_copy]
else:
  logging.error("Uns source 8 AnnData file is not compatible to be merged to main AnnData file, different cell names.")
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
    
