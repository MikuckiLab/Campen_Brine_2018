# Campen_Brine_2018
# 4_EMB_core_microbiome

"""
This file deals with the sequencing processing of samples for characterising
Blood Falls outflow brine in comparison to non-outflow samples collected 
on/near Blood Falls, and samples from the deeper water column of West lobe
Lake Bonney. 

This file describes the process of determining the EMB core microbiome.

To run this analysis you will require:

* mothur shared file
* mothur constaxonomy file
* mothur database file
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

#========== Determine EMB Core Microbiome ==========#

# get EMB shared OTUs
emb_sxp = sxp.loc[:,sxp.metadata['type'] == 'End Member Brine']
emb_sxp_relabund = emb_sxp.relabund(100)
emb_shared_relabund = emb_sxp_relabund.loc[emb_sxp_relabund.features.min(axis=1) > 0].relabund(100)

# calculate descriptive statistics
otu_means = pd.DataFrame(emb_shared_relabund.features.mean(axis=1), columns=['mean'])
otu_stds = pd.DataFrame(emb_shared_relabund.features.std(axis=1), columns=['std'])

# combine sample abundances with stats and format columns
emb_abunds = pd.concat([
    np.around(emb_shared_relabund.features, 2),
    np.around(otu_means, 2),
    np.around(otu_stds, 2)
], axis=1)
emb_abunds.columns = pd.MultiIndex.from_tuples([('OTU Abundances (%)', x) for x in emb_abunds.columns])

# get classifications and format columns
emb_classifications = emb_shared_relabund.classifications.copy(deep=True)
emb_classifications.columns = pd.MultiIndex.from_tuples([('Silva v128 Classifications', x) for x in emb_classifications.columns])

# combine abundances and clasifications and sort in descending order by mean abundance
core_otus = pd.concat([
    emb_abunds,
    emb_classifications
], axis=1)
core_otus = core_otus.sort_values(('OTU Abundances (%)', 'mean'), ascending=False)

# save to file and display
with open(os.path.join(results_dir, 'emb_core_mircobiome.csv'), 'w') as out_handle:
    
    # write header line
    out_handle.write('\"End Member Brine core microbiome OTUs (99% level) relative abundances, mean abundances ' +
                     'and standard deviation, and Silva v128 classifications:\"\n')
    
    # write out data
    core_otus.to_csv(out_handle)

#---------- Add EMB core micobiome membership into corrected database file ----------#

# read in database and replace taxonomies with the new ones, also dropping unknowns
tio_db = pd.read_csv(os.path.join(processed_dir, 'TIO.database.csv'), index_col=0)
abunds = tio_db.iloc[:, 0:20].copy(deep=True)
seqs = tio_db[['repSeq']].copy(deep=True)
classifications = sxp.classifications.copy(deep=True)
core_membership = pd.DataFrame(tio_db.index.isin(core_otus.index), index=tio_db.index)

# format columns
abunds.columns = pd.MultiIndex.from_tuples([('OTU Abundances (99%)', x) for x in abunds.columns])
seqs.columns = pd.MultiIndex.from_tuples([('', 'RepSeq')])
classifications.columns = pd.MultiIndex.from_tuples([('Silva v128 classifications', x) for x in classifications.columns])
core_membership.columns = pd.MultiIndex.from_tuples([('', 'EMBCore')])

# combine back together and add rep sequences
tio_db = pd.concat([core_membership, abunds, classifications, seqs], axis=1).dropna()
tio_db.index.name = 'OTU #'


# save to file
results_file = os.path.join(results_dir, 'corrected_database.csv')
with open(results_file, 'w') as out_handle:
    
    # print header text then write data to file
    out_handle.write('\"Corrected mothur database containing OTU (99%) abundances, EMB core microbiome membership, Silva v128 ' +
                     'classifications, and representitive fasta sequences:\"\n')
    tio_db.to_csv(out_handle)    

#========== END ==========#
