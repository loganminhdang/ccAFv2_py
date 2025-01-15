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

# docker run -it -v '/home/cplaisier/ccAFv2_py:/files' cplaisier/ccafv2_py

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
PCW8 = sc.read_h5ad('../data/W8-1_normalized_ensembl.h5ad')

# Run ccAFv2 to predict cell labels
PCW8_labels = ccAFv2.predict_labels(PCW8, threshold=0.9, species='human', gene_id='ensembl')

# Save into scanpy object
PCW8.obs['ccAFv2'] = PCW8_labels[0]

# Run UMAP of U5 hNSCs
sc.pp.highly_variable_genes(PCW8, n_top_genes=2000)
sc.tl.pca(PCW8)
sc.pp.neighbors(PCW8)
sc.tl.umap(PCW8)

# Dictionary to provdie colors to plotting
cmap1 = {"G0/G1": "#e34234", "Neural G0": "#d9a428", "G1": "#f37f73", "Late G1": "#1fb1a9",  "S": "#8571b2", "S/G2": "#db7092", "G2/M": "#3db270" ,"M/Early G1": "#6d90ca",  "Unknown": "#d3d3d3"}

# Plot UMAP of U5 hNSCs
sc.pl.umap(PCW8, color=['ccAFv2'], palette=cmap1, save='ccAFv2_UMAP_PCW8.pdf')

