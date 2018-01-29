# Campen_Brine_2018
# 3_estimate_alpha_diversities

"""
This file deals with the sequencing processing of samples for characterising
Blood Falls outflow brine in comparison to non-outflow samples collected 
on/near Blood Falls, and samples from the deeper water column of West lobe
Lake Bonney. 

This file describes the process of estimating alpha diversity metrics for the samples.
It also descibes correlating these metrics with sequence counts and geochemical parameters.

To run this analysis you will require:

* mothur shared file
* mothur constaxonomy file
* sample metadata file

"""

#========== General Setup ==========#

# import necessary libraries
import math
import os
import shutil

from collections import OrderedDict
from copy import deepcopy
from decimal import Decimal
from itertools import combinations, product

import numpy as np
import pandas as pd
import seq_experiment as se

from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.spatial.distance import squareform

from seq_experiment.io.mothur import MothurIO
from seq_experiment.utils import mothur

# setup directories
# NOTE: assumes code is being run from the code directory
base_dir = os.path.join('..')
base_dir = os.path.abspath(base_dir)
data_dir = os.path.join(base_dir, 'data')

processed_dir = os.path.join(data_dir, 'processed')
results_dir = os.path.join(base_dir, 'results')

#========== Read in Data ==========#

# create SeqExp object from mothur files
sxp = se.SeqExp(
    features=MothurIO.read_shared_file(os.path.join(processed_dir, 'TIO.shared')),
    classifications=MothurIO.read_constaxonomy_file(os.path.join(processed_dir, 'TIO.cons.taxonomy')),
)

# fix classification names and remove consensus agreement numbers
sxp.classifications.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
sxp.classifications = mothur.tidy_mothur_classifications(sxp.classifications)

# read in and format metadata before mergning with SeqExp
metadata = pd.read_csv(os.path.join(data_dir, 'sample_metadata.csv'), index_col=0)
metadata.index.name = None
sxp = sxp.merge(component='metadata', right=metadata)
# set date column to type datetime
sxp.metadata['date'] =  pd.to_datetime(sxp.metadata['date'])

# remove OTUs classified as unknown from the SeqExp
sxp = sxp.drop(by='features', items=sxp[sxp.classifications['Kingdom'] == 'unknown'].feature_names)

# rename OTUs of Class Chloroplast to Phylum Chloroplast so it is more apparant what they are, 
# instead of grouping them under cyannobacteria
new_classifications = sxp.classifications
for row in sxp.classifications.iterrows():
    feature = row[0]
    data = row[1]
    
    # fix Chloroplasts
    if data['Class'] == 'Chloroplast':
        new_classifications.set_value(col='Kingdom', index=feature, value='Chloroplast')
        new_classifications.set_value(col='Phylum', index=feature, value='Chloroplast_ph')
        new_classifications.set_value(col='Class', index=feature, value='Chloroplast_cl')
        
    if data['Phylum'] == 'Woesearchaeota_(DHVEG-6)':
        for level in sxp.classifications.columns[2:]:
            new_classifications.set_value(col=level, index=feature, value='Woesearchaeota_(DHVEG-6)_%s' % level.lower()[:2])
        
sxp.classifications = new_classifications
    
def replace_tax(column, match, replace):
    
    new_col = column
    for idx, value in enumerate(column):
        if match in value:
            new_col[idx] = column[idx].replace(match, replace)
            
    return new_col

sxp.classifications.apply(lambda col: replace_tax(col, 'Woesearchaeota_(DHVEG-6)', 'DHVE6'))

# generate prokaryote only subset sequence data
sxp_prok = deepcopy(sxp).drop(by='features', items=sxp[sxp.classifications['Kingdom'].isin(['Eukaryota', 'Chloroplast'])].feature_names)

#========== Estimate Alpha Diversity Metrics ==========#

def shannon_index(x):
    """Estimate Shannon diversity."""
    
    # eliminate zeros to remove OTUs not present in this sample
    x = x[x != 0]   
    
    # convert to proportions
    x = x / (sum(x))   
    
    # multiply by ln
    x = -x * np.log(x)    
    
    # sum for shannon
    H = x.sum()
    
    return H
    
def hill_diversity(x, q=1):
    """
    Calculates Hill's diversity numbers for a given quotient.

    Hill's numbers are also called `Effective Species Number` and are the number of features
    required in an evenly abundant population to acheive a given observed sample diversity.

    """
    
    # eliminate zeros to remove features not present in this sample
    x = x[x != 0]
    
    # get proportions
    x = np.divide(x, np.sum(x))
    
    if q == 1:
        # cannot calculate when q is 1 as it causes a divide by zero error
        # instead we approximate from the value obtained as q approaches 1
        # which is just the exponent of the shannon diversity
        return np.exp(shannon_index(x))
    else:
        x = np.power(x, q)
        sum_ = np.sum(x)
        div = np.power(sum_, 1/(1-q))   
        return div
    
def simpson_index(x):    
    # is just the reciprical of hill_diversity number where q=2
    return 1 / hill_diversity(x, q=2)

def calculate_alphas(sxp_):
    """Estimates alpha diversities for samples in a SeqExp object."""

    # init resutls container for all metrics across all samples
    alpha_lists = []

    # iterate over each sample in the SeqExp
    for samp in sxp_.sample_names:

        # init results container for the alpha diversity metrics for this sample
        alpha_list = []
        
        # calcualte total sequenes
        alpha_list.append(sxp_.features[samp].sum())

        # calculate number of features, Shannon, and Simpson indexes
        alpha_list.append(hill_diversity(sxp_.features[samp], 0))
        alpha_list.append(shannon_index(sxp_.features[samp]))
        alpha_list.append(simpson_index(sxp_.features[samp]))

        # add in hills numbers
        for i in [1, 2]:
            div = hill_diversity(sxp_.features[samp], i)
            alpha_list.append(div)

        alpha_lists.append(alpha_list)

    # format results into pandas.DataFrame and return
    results_df = pd.DataFrame(alpha_lists, index=sxp.sample_names, columns=['Number Seqs', 'Number OTUs', 'Shannon', 'Simpson', 'Hills q=1', 'Hills q=2'])
    return results_df

def format_alpha_values(df_):
    """Conditionally formats number of decimal places for alpha diversity data."""
    
    # init results dataframe
    new_df = df_.copy(deep=True)
    
    # tidy the formatting of the numbers in the diversity indicies
    for index, row in df_.iterrows():
        for col in row.index:
            
            # get decimals
            num = df_.get_value(index, col)
            parts = str(num).split('.')
            
            # format decimal places and update in results dataframe if it is a float
            if len(parts) > 1:
                if len(parts[1]) > 2:
                    num = np.around(num, 2)
                    new_df.set_value(index, col, num)

    return new_df

# calcualte alpha metrics and format into a single df
alphas_all_df = calculate_alphas(sxp)
alphas_all_df = format_alpha_values(alphas_all_df)
alphas_all_df['Number OTUs'] = alphas_all_df['Number OTUs'].astype(int)
alphas_all_df.columns = pd.MultiIndex.from_tuples([('All', column) for column in alphas_all_df.columns])

alphas_prok_df = calculate_alphas(sxp_prok)
alphas_prok_df = format_alpha_values(alphas_prok_df)
alphas_prok_df['Number OTUs'] = alphas_prok_df['Number OTUs'].astype(int)
alphas_prok_df.columns = pd.MultiIndex.from_tuples([('Prokaryotes', column) for column in alphas_prok_df.columns])
        
alphas_df = pd.concat([alphas_all_df, alphas_prok_df], axis=1)
alphas_df.index.name = 'Sample'

# save to the results folder
with open(os.path.join(results_dir, 'alphas.csv'), 'w') as out_handle:
    
    # print header line
    out_handle.write('\"Alpha diversities of samples collected at Blood Falls and the West lobe of Lake Bonney, negative controls, and mock community:\"\n')
    alphas_df.to_csv(out_handle)

alphas_df

#========== END ==========#
