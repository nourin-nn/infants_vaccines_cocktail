# matplotlib inline
# use Python 3.8.3 64-bit (base:conda)
import sys
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as scn

samples = ['JB19050', 'JB19051', 'JB19052', 'JB19053', 'JB19054', 'JB19055', 
           'JB19058', 'JB19059', 'JB19060', 'JB19061', 'JB19062', 'JB19063', 
           'JB19066', 'JB19067', 'JB19068', 'JB19069', 'JB19070', 'JB19071']
ind_data_path = os.path.join(project_path, 'data/Individual_samples')
analysis_path = os.path.join(project_path, 'results/analysis')
for smpl in samples:
    counts_matrix = scipy.io.mmread(os.path.join(ind_data_path, smpl, 'filtered_feature_bc_matrix/matrix.mtx.gz')).T.tocsc()
    print('Counts matrix shape: {} rows (cells), {} columns (genes)'.format(counts_matrix.shape[0], counts_matrix.shape[1]))    
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    print('scrub_doublets output: {} doublet_scores (cells), {} predicted_doublets (cells)'.format(len(doublet_scores), len(predicted_doublets)))
    adata = scn.read_10x_mtx(os.path.join(ind_data_path, smpl, 'filtered_feature_bc_matrix'), var_names='gene_symbols') 
    res_summ_obs = pd.DataFrame(list(zip(adata.obs_names, doublet_scores, predicted_doublets)), 
                                columns=['barcode','score_obs', 'doublet'])
    res_summ_sim = pd.DataFrame(list(zip(scrub.doublet_scores_sim_)), 
                                columns=['score_sim'])
    res_summ_obs.to_csv(os.path.join(analysis_path, 'obs_doublets_' + smpl + '.tsv'), index=False, sep="\t")
    res_summ_sim.to_csv(os.path.join(analysis_path, 'sim_doublets_' + smpl + '.tsv'), index=False, sep="\t")
