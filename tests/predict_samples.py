##########################################################
## OncoMerge:  predict_samples.py                       ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

##########################################
## Load Python packages for classifiers ##
##########################################

# General
import pandas as pd
import scanpy as sc
import ccAFv2


#####################
## Test prediction ##
#####################

# Load up test dataset
U5hNSC = sc.read_h5ad('../data/U5_normalized_ensembl_all_genes.h5ad')

# Run ccAFv2 to predict cell labels
U5hNSC_labels = ccAFv2.predict_labels(U5hNSC, species='human', gene_id='ensembl')

# Save into scanpy object
U5hNSC.obs['ccAFv2'] = pd.Categorical(U5hNSC_labels[0], categories=['Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown'], ordered=True)

# Run UMAP of U5 hNSCs
sc.pp.highly_variable_genes(U5hNSC, n_top_genes=3000)
sc.tl.pca(U5hNSC)
sc.pp.neighbors(U5hNSC)
sc.tl.umap(U5hNSC)
sc.tl.tsne(U5hNSC)

# Reordering
#U5hNSC.obs['ccAFv2'] = U5hNSC.obs['ccAFv2'].cat.reorder_categories(['Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1', 'Unknown'], ordered=True)

cmap1 = {"Neural G0": "#d9a428", "G1": "#f37f73", "Late G1": "#1fb1a9",  "S": "#8571b2", "S/G2": "#db7092", "G2/M": "#3db270" ,"M/Early G1": "#6d90ca",  "Unknown": "#d3d3d3"}

# Plot UMAP of U5 hNSCs
sc.pl.tsne(U5hNSC, color=['ccAFv2'], palette=cmap1, save='ccAFv2_tSNE.pdf')
sc.pl.umap(U5hNSC, color=['ccAFv2'], palette=cmap1, save='ccAFv2_UMAP.pdf')



#PCW8 = sc.read_h5ad('../data/W8-1_normalized_ensembl.h5ad')
#PCW8_labels = ccAFv2.predict_labels(PCW8, species='human', gene_id='ensembl')

