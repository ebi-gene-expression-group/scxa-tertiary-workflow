import scanpy as sc
from numpy import all
import anndata
import logging
import sys

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
