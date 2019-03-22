# Campen_Brine_2018
# 2_correlate_geochemistry

"""
This file deals with the sequencing processing of samples for characterising
Blood Falls outflow brine in comparison to non-outflow samples collected 
on/near Blood Falls, and samples from the deeper water column of West lobe
Lake Bonney. 

This file describes the process of correlating the geochemical parameters of the brine.

To run this analysis you will require:

* sample metadata file

"""

#========== General Setup ==========#

# import necessary libraries
import os
import shutil

from collections import OrderedDict
from itertools import combinations

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.spatial.distance import squareform

# setup directories
# NOTE: assumes code is being run from the code directory
base_dir = os.path.join('..')
base_dir = os.path.abspath(base_dir)
data_dir = os.path.join(base_dir, 'data')

processed_dir = os.path.join(data_dir, 'processed')
results_dir = os.path.join(base_dir, 'results')

#========== Calculate correlations ==========#

def correlate_columns(df, combos=None, round_=3):
    """Calulcates spearmans rho and significance for all column combinations in a dataframe."""
    
    # get all pairwise combinations to compare id user has not specified
    if combos is None:
        combos = list(combinations(df.columns, 2))

    # init results dataframes       
    index = cols = df.columns
    data = pd.DataFrame(np.zeros(shape=(len(index), len(cols))), columns=cols, index=index)  # dummy data
    r_df = data.copy(deep=True)
    p_df = data.copy(deep=True)
        
    # correlate columns, iterating over tuples of column combinations
    for x_y in combos:

        # combine columns of interest into new pd Dataframe for easy manipulation
        # drop any samples that don't have both of the parameters being compared
        compare_df = pd.concat([df[x_y[0]], df[x_y[1]]], axis=1).dropna()

        # calucalte statistics and save to results containers
        r, p = stats.spearmanr(compare_df.iloc[:,0], compare_df.iloc[:,1])
        r_df.set_value(x_y[1], x_y[0], r)
        p_df.set_value(x_y[1], x_y[0], p)
    
    return r_df, p_df

def correct_p_vals(df, correction_method='fdr_bh'):
    """Corrects p-values."""
    
    corrected_p = multipletests(df.values.flatten(), method=correction_method)[1]
    corrected_p_df = np.array(corrected_p).reshape((len(df.index), len(df.columns)))
    corrected_p_df = pd.DataFrame(corrected_p_df, index=df.index, columns=df.columns)
    
    return corrected_p_df

# read in geochemistry from sample metadata for BF samples only (not negatives or mock)
metadata = pd.read_csv(os.path.join(data_dir, 'sample_metadata.csv'), index_col=0)
geochem = metadata.iloc[:, 4:]

# correlate geochemical parameters
from decimal import Decimal

# get spearman correlations
geo_corr_r, geo_corr_p = correlate_columns(geochem)
geo_corr_p = correct_p_vals(geo_corr_p)

# combine results into single df where correlation is lower triangle and p val in upper triangle
geo_corr_res = pd.DataFrame((geo_corr_r + geo_corr_p.transpose()), index=geochem.columns, columns=geochem.columns)

# format values to ensure they are to 3 decimal places
places = Decimal(10) ** -3 # same as Decimal('0.001')
geo_corr_res = geo_corr_res.apply(lambda x: [Decimal(i).quantize(places) for i in x])

# change values on diagonal to empy string
for i in range(len(geochem.columns)):
    geo_corr_res.iloc[i,i] = ''

# save output
with open(os.path.join(results_dir, 'geo_corr.csv'), 'w') as out_handle:
    
    # write header text
    out_handle.write(
        '\"Spearmans\'s correlations (rho, lower triangle) between geoochemical paramertes collected for Blood Falls ' +
        'and West lobe Lake Bonney samples, with p-values corrected for false discovery rate by ' +
        'the Benjamini-Hochberg method (upper triangle):\"\n'
    )

    # write out data
    geo_corr_res.to_csv(out_handle)

    # write out additional text
    out_handle.write(
			'\n\"DIC: Dissolved Inorganic Carbon, DOC: Dissolved Organic Carbon, DIN: Dissolved ' +
			'Inorganic Nitrogen, DO: Dissolved Oxygen\"\n'
    	)

#========== END ==========#
