# Campen_Brine_2018
# 5_identify_biomarkers

"""
This file deals with the sequencing processing of samples for characterising
Blood Falls outflow brine in comparison to non-outflow samples collected 
on/near Blood Falls, and samples from the deeper water column of West lobe
Lake Bonney. 

This file describes the process of identifying taxonomic biomarkers of the EMB and NOF samples
using Linear discriminant analysis effect size (LefSe). We will do this using mothur's
implementation of the LefSe algorithm. We will also correlate the abundances of the taxa
identified as biomakers with the geochemical parameters.

To run this analysis you will require:

* Mothur executable
* sample metadata file
* mothur shared file
* mothur constaxonomy file
* sample metadata file

"""

#========== General Setup ==========#

# import necessary libraries
import itertools
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

from mothur_py import Mothur
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

#========== Mothur Setup ==========#

tmp_dir =  os.path.join(data_dir, 'tmp')
if not os.path.isdir(tmp_dir):
    os.mkdir(tmp_dir)

mothur_out_dir = os.path.join(tmp_dir, '4_identify_biomarkers')
if not os.path.isdir(mothur_out_dir):
    os.mkdir(mothur_out_dir)

# create Mothur instance and set configuration
m = Mothur()
m.verbosity = 1
m.suppress_logfile = True
m.mothur_seed = 4321
m.current_dirs['output'] = m.current_dirs['tempdefault'] = mothur_out_dir

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

# drop unwated samples and then drop features that are absent from all samples
sxp = sxp.drop(by='samples', items=sxp.loc[:, sxp.metadata['type'].isin(['Mock', 'Negative'])].sample_names)
sxp = sxp.loc[sxp.features.max(axis=1) > 0]

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

# create subset of full data without Eukaryotic sequences or Chloroplast sequences
sxp_prok = deepcopy(sxp)
sxp_prok = sxp_prok.loc[sxp_prok.classifications['Kingdom'] != 'Eukaryota']
sxp_prok = sxp_prok.loc[sxp_prok.classifications['Kingdom'] != 'Chloroplast']

# get EMB and NOF from of full sequence data
emb_nof_sxp = sxp.loc[:, sxp.metadata['type'].isin(['End Member Brine', 'Non Outflow'])]
emb_nof_sxp = emb_nof_sxp.loc[emb_nof_sxp.features.max(axis=1) > 0]

# get EMB and NOF from of prokaryote only sequence data
emb_nof_sxp_prok = sxp_prok.loc[:, sxp_prok.metadata['type'].isin(['End Member Brine', 'Non Outflow'])]
emb_nof_sxp_prok = emb_nof_sxp_prok.loc[emb_nof_sxp_prok.features.max(axis=1) > 0]

#========== Lefse ==========#

def to_mothur_shared(sxp_, label, out_file):
        """Exports SeqExp objects features to a Mothur shared file."""

        shared = sxp_.features
        shared = shared.transpose()
        shared = shared.reset_index()
        
        shared['label'] = label
        shared['numOtus'] = len(sxp_.features)
        new_column_order = ['label', 'Group', 'numOtus', *sxp_.features.index]
        shared = shared[new_column_order]

        shared.to_csv(out_file, sep='\t', header=True, index=False)

        return

def generate_shared_files(sxp, levels, out_dir, file_base, otu_level='0.01'):
    """
    Generates shared files for a SeqExp object at a range of taxonomic levels.
    
    Groups SeqExp objects at the taxonomic levels specified, and saves as shared files using
    the `file_base` prefix to the `out_dir`.
    
    Returns a tuple of the path to the shared file, and the assosiated SeqExp object.
    
    """

    # init results containers
    grouped_sxps = OrderedDict()
    shared_files = OrderedDict()
    
    # create shared files at the different taxonomic levels except Kingdom, but including OTU level
    for level in levels:
        
        sxp_ = deepcopy(sxp)
        
        # get SeqExp grouped at the given taxonomic level and save for later reference
        # also set label for shared file
        if level == 'OTU':
            sxp_grouped = sxp_
            label = otu_level
        else:
            sxp_grouped = sxp_.groupby_classification(level)
            label = level

        sxp_grouped.features.index.name = None
        sxp_grouped.features.columns.name = 'Group'

        # write shared file
        out_file = os.path.join(out_dir, '{0}_{1}.shared'.format(file_base, label))
        to_mothur_shared(sxp_grouped, label, out_file)
        
        # save to results conainers
        grouped_sxps[label] = sxp_grouped
        shared_files[label] = out_file
                
    # return grouped SeqExp objects
    return grouped_sxps, shared_files

def run_lefse(shared_files, mothur_obj, out_dir, file_base):
    """Runs mothur's LefSe implementation for the SeqExp at the specified taxonomic levels"""
    
    # init results conatiners
    mothur_objs = OrderedDict()
        
    # run lefse
    for label, shared_file in shared_files.items():
               
        # make copy of mothur object and feed in dirs/files and save to results containers
        m_ = deepcopy(m)
        m_.current_dirs['output'] = mothur_out_dir

        # run mothur and save resulting object to results container along with the SeqExp
        m_.lefse(shared=shared_file, design=design_file)
        mothur_objs[label] = m_
        
    return mothur_objs

def read_lefse_summaries(results):
    """Reads in LefSe results."""

    # init results container
    formatted_summaries = OrderedDict()
    
    # init results containers
    # lefse_summaries = OrderedDict()

    # iterate levels that LefSe was performed on
    for label in results.keys():

        # read in LefSe summary, format, and save to result container
        summary = pd.read_table(os.path.join(mothur_out_dir, '{0}_{1}.{1}.lefse_summary'.format(FILE_BASE_ALL, label)))
        summary = summary.set_index('OTU')
        summary.index.name = None
        formatted_summaries[label] = summary
        
    return formatted_summaries

def get_top_lefse_hits(summaries, lda_cutoff, p_cutoff):

    #init results container
    top_hits = OrderedDict()

    # iterate levels that LefSe was performed on
    for label in summaries.keys():

        # filter LefSe results by our criteria and save to results container
        summary = summaries[label]
        filtered_summary = summary[summary['LDA'] > LDA_CUTOFF]
        filtered_summary = summary[summary['pValue'] < PVAL_CUTTOF]

        top_hits[label] = filtered_summary
        
    return top_hits

def split_lefse_by_type(summaries, types):
        
    #init results container
    new_summaries = OrderedDict()
    for type_ in types:
        new_summaries[type_] = OrderedDict()
    
    for label in summaries.keys():
        
        # get summary
        summary = summaries[label]
        
        # split results by sample type
        for type_ in types: 
            new_summary = summary[summary['Class'] == type_]
            new_summaries[type_][label] = new_summary
    
    return new_summaries

def filter_lefse_abundance(summaries, grouped_sxps, types, abund_cutoff, class_mapping):
    
    #init results container
    new_summaries = OrderedDict()
    for type_ in types:
        new_summaries[type_] = OrderedDict()
        
    # iterate each level for each type
    for type_ in types:
        for label in summaries[type_].keys():
            
            # get filtered summary
            summary = summaries[type_][label]

            # get taxa whose mean abundance is greater than the abundance cutoff
            sxp_ = grouped_sxps[label].relabund(100)
            sxp_ = sxp_.loc[:, sxp_.metadata['type'] == class_mapping[type_]]
            sxp_ = sxp_.loc[summary.index]
            mean = sxp_.features.mean(axis=1)
            abundant_taxa = mean[mean > abund_cutoff]

            # filter summary by abundant taxa, then add in abundance data
            new_summary = summary.loc[abundant_taxa.index]
            new_summary['MeanRelativeAbundance'] = abundant_taxa.values

            # save results
            new_summaries[type_][label] = new_summary
    
    return new_summaries

def format_lefse_results(summaries, types):
    
    # init results container
    level_groups = OrderedDict()

    # group data at each level
    for type_ in types:
        for label in summaries[type_].keys():
            
            if level_groups.get(label) is None:
                level_groups[label] = OrderedDict()

            # get data for each sample type
            summary = summaries[type_][label].copy(deep=True)
            summary.columns = pd.MultiIndex.from_tuples([(type_, x) for x in summary.columns])
            
            level_groups[label][type_] = summary
            
    # combine data for each sample type and format
    merged_groups = OrderedDict()
    for label in level_groups.keys():
        
        merged_group = pd.concat([level_groups[label][type_] for type_ in level_groups[label].keys()], axis=1)
        
        # compact down what sample the taxa is a biomarker for into a single columns
        biomarker_type = pd.DataFrame(merged_group[('EMB', 'Class')].fillna('NOF'))  # shortcut as we know the column names
        biomarker_type.columns = pd.MultiIndex.from_tuples([('', 'BiomarkerSampleType')])
        
        merged_group = pd.concat([biomarker_type, merged_group], axis=1)
        merged_group = merged_group.drop([(x, 'Class') for x in types], axis=1)
        merged_group.index.name = label
        
        merged_groups[label] = merged_group
        
    return merged_groups
    

def run_lefse_pipeline(sxp, mothur_obj, levels, out_dir, file_base):

    # generate SeqExp's grouped at the specified taxonomic levels and the respective shared files
    grouped_sxps, shared_files = generate_shared_files(
        sxp=sxp, 
        levels=levels,
        out_dir=out_dir,
        file_base=file_base
    )
    # run lefse
    lefse_results = run_lefse(
        shared_files=shared_files,
        mothur_obj=mothur_obj,
        out_dir=out_dir,
        file_base=file_base
    )
    lefse_summaries = read_lefse_summaries(results=lefse_results)
    lefse_filtered_hits = get_top_lefse_hits(lefse_summaries, lda_cutoff=LDA_CUTOFF, p_cutoff=PVAL_CUTTOF)
    lefse_filtered_type_hits = split_lefse_by_type(lefse_filtered_hits, types=['EMB', 'NOF'])
    lefse_filtered_type_abund_hits = filter_lefse_abundance(lefse_filtered_type_hits, grouped_sxps, types=['EMB', 'NOF'], 
                                                            abund_cutoff=ABUND_CUTOFF, class_mapping=TYPE_MAPPING)
    formatted_results = format_lefse_results(lefse_filtered_type_abund_hits, types=['EMB', 'NOF'])
    
    return formatted_results

# create design file for mothur and format names
design = emb_nof_sxp.metadata['type']
design = design.apply(lambda x: 'EMB' if x == 'End Member Brine' else 'NOF')

design_file = os.path.join(mothur_out_dir, 'emb_v_nof.design')
design.to_csv(design_file, sep='\t', header=False)

# init constants
LDA_CUTOFF = 2
PVAL_CUTTOF = 0.05
ABUND_CUTOFF = 0.5
TYPE_MAPPING = {
    'EMB': 'End Member Brine',
    'NOF': 'Non Outflow'
}
# setup constants
LEVELS = emb_nof_sxp.classifications.columns[1:].tolist()
LEVELS.append('OTU')
OTU_LEVEL = '0.01'
FILE_BASE_ALL = 'emb_nof'
FILE_BASE_PROK = 'emb_nof_prok'

# run LefSe using full sequence data set
lefse_res_all = run_lefse_pipeline(
    sxp=emb_nof_sxp,
    mothur_obj=m,
    levels=LEVELS,
    out_dir=mothur_out_dir,
    file_base=FILE_BASE_ALL
)

# run LefSe using prokaryote only sequence data set
lefse_res_prok = run_lefse_pipeline(
    sxp=emb_nof_sxp_prok,
    mothur_obj=m,
    levels=LEVELS,
    out_dir=mothur_out_dir,
    file_base=FILE_BASE_ALL
)

# combine data from the full data set and the prokaryote only data set
# init results container
combined_lefse_results = OrderedDict()

for label in lefse_res_all.keys():

    all_data = lefse_res_all[label].copy(deep=True)
    all_data.columns = pd.MultiIndex.from_tuples([('All Sequences', x[0], x[1]) for x in all_data.columns])
    prok_dat = lefse_res_prok[label].copy(deep=True)
    prok_dat.columns = pd.MultiIndex.from_tuples([('Prokaryote Sequences', x[0], x[1]) for x in prok_dat.columns])

    combined_lefse_result = pd.concat([all_data, prok_dat], axis=1).fillna('NA')
    combined_lefse_result.index.name = label
    combined_lefse_results[label] = combined_lefse_result

#========== Correlate LefSe Biomarker Abundance and Geochemical Parameters ==========#

def correlate_columns(df, combos=None, round_=3):
    """Calulcates perasons r and significance for all column combinations in a dataframe."""
    
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

    # tidy values
    r_df, p_df = r_df.round(round_), p_df.round(round_)
    
    return r_df, p_df

def correct_p_vals(df, correction_method='fdr_bh'):
    """Corrects p values."""
    
    # we need to correct the p values to account for the false discovery rate assosiated with multiple testing
    corrected_p = multipletests(df.values.flatten(), method=correction_method)[1]
    corrected_p_df = np.array(corrected_p).reshape((len(df.index), len(df.columns)))
    corrected_p_df = pd.DataFrame(corrected_p_df, index=df.index, columns=df.columns)
    
    return corrected_p_df

STARS_DICT = OrderedDict([('*', 0.05), ('**', 0.01), ('***', 0.001)])

def add_sig_stars(r_df, p_df, stars_dict=STARS_DICT):
    """Conditonally adds starts to correlation values to indicate differing levels of significance."""

    new_r_df = r_df.copy(deep=True).astype(str)

    for index, row in r_df.iterrows():
        for col in row.index:
            
            # get values using Decimal to maintain the correct number of decimal places
            p = p_df.get_value(index, col)
            r = r_df.get_value(index, col)
            
            places = Decimal(10) ** -3 # same as Decimal('0.01')
            r_str = str(Decimal(r).quantize(places))

            # conditonally add stars to correlation 
            if p < stars_dict['*']:
                r_str += '*'
            if p < stars_dict['**']:
                r_str += '*'
            if p < stars_dict['***']:
                r_str += '*'

            # set new df val
            new_r_df.set_value(index, col, r_str)
        
    return new_r_df

# setup geochem
geochem = emb_nof_sxp.metadata.iloc[:, 4:]

# perform correlations
# init results container
lefse_corrs = OrderedDict()
lefse_corrs['all'] = OrderedDict()
lefse_corrs['prok'] = OrderedDict()

# iterate over labels
for label in combined_lefse_results.keys():
    
    # get lefse data
    lefse_hits_all = lefse_res_all[label]
    lefse_hits_prok = lefse_res_prok[label]
    
    # get abundance data, and orient correctly
    if label == OTU_LEVEL:
        sxp_all_ = deepcopy(sxp)
        sxp_prok_ = deepcopy(sxp_prok)
    else:
        sxp_all_ = deepcopy(sxp).groupby_classification(label)
        sxp_prok_ = deepcopy(sxp_prok).groupby_classification(label)
        
    abund_all = sxp_all_.loc[lefse_hits_all.index].features.transpose()
    abund_prok = sxp_prok_.loc[lefse_hits_prok.index].features.transpose()

    for data_set, df in [('all', abund_all), ('prok', abund_prok)]:    
        # correlate geochemistry with taxonomic abundances and format results
        r_df, p_df = correlate_columns(pd.concat([geochem, df], axis=1), combos=itertools.product(geochem, df))
        r_df = r_df.loc[df.columns, geochem.columns]
        p_df = correct_p_vals(p_df)
        corrs = add_sig_stars(r_df, p_df)

        lefse_corrs[data_set][label] = corrs

# add in extra information to the correlation tables including type the biomarker is for and dataset it is for
# init results container
combined_corrs = OrderedDict()

# iterate over labels
for label in combined_lefse_results.keys():
    
    # get correlation data
    corr_all = lefse_corrs['all'][label]
    corr_prok = lefse_corrs['prok'][label]
    
    # get lefse results
    lefse_all = lefse_res_all[label]
    lefse_prok = lefse_res_prok[label]
    
    # add in what data set they were from and reorder columns
    corr_all['BiomarkerSampleType'] = lefse_all[('', 'BiomarkerSampleType')]
    corr_all = corr_all[['BiomarkerSampleType'] + corr_all.columns[:-1].tolist()]
    corr_prok['BiomarkerSampleType'] = lefse_prok[('', 'BiomarkerSampleType')]
    corr_prok = corr_prok[['BiomarkerSampleType'] + corr_prok.columns[:-1].tolist()]
        
    # get overlapping and non-overlapping features
    shared_features = corr_all.index.intersection(corr_prok.index).tolist()
    all_only = corr_all.index.difference(shared_features).tolist()
    prok_only = corr_prok.index.difference(shared_features).tolist()
    
    # subset correlation data by new feature lists
    corr_shared = corr_all.loc[shared_features]
    corr_all_only = corr_all.loc[all_only]
    corr_prok_only = corr_prok.loc[prok_only]
    
    # add in utility columns
    corr_shared['group'] = '0_shared'
    corr_all_only['group'] = '1_all'
    corr_prok_only['group'] = '2_prok'
    
    # format index's to multilevel
    corr_shared.index = pd.MultiIndex.from_tuples([('Shared', x) for x in corr_shared.index])
    corr_all_only.index = pd.MultiIndex.from_tuples([('All only', x) for x in corr_all_only.index])
    corr_prok_only.index = pd.MultiIndex.from_tuples([('Prokaryotes only', x) for x in corr_prok_only.index])
    
    # combine and save to result container
    combined = pd.concat([corr_shared, corr_all_only, corr_prok_only])
    combined['feature'] = combined.index.tolist()
    combined = combined.sort_values(['group', 'BiomarkerSampleType', 'feature'])
    combined = combined.drop(['group', 'feature'], axis=1)
    combined.index.names = ['DataSet', 'Feature']
    
    combined_corrs[label] = combined

# save intermediate data to processed dir
for label in combined_corrs.keys():
    df = combined_corrs[label]
    df.astype(str).to_csv(os.path.join(processed_dir, 'lefse_%s_corr.csv' % label))

# merge LefSe scores with correlation results
# init results container
lefse_res_corrs = OrderedDict()

# iterate over labels
for label in combined_lefse_results.keys():
    
    # get correlations and format
    correlations = combined_corrs[label]
    correlations = correlations.reset_index()
    correlations = correlations.set_index('Feature')
    correlations.columns = pd.MultiIndex.from_tuples([('Geochemistry Correlations (Spearman\'s rho)', '', x) for x in correlations.columns])
    
    biomaker_sample_types = correlations[[('Geochemistry Correlations (Spearman\'s rho)', '', 'BiomarkerSampleType')]]
    correlations = correlations.drop(biomaker_sample_types.columns, axis=1)
    biomaker_sample_types.columns = pd.MultiIndex.from_tuples([('', '', x[-1]) for x in biomaker_sample_types.columns])   
    
    # get lefse results and correlations
    lefse_results = combined_lefse_results[label]
    lefse_results_all = lefse_results['All Sequences'].drop(('', 'BiomarkerSampleType'), axis=1)
    lefse_results_prok = lefse_results['Prokaryote Sequences'].drop(('', 'BiomarkerSampleType'), axis=1)
    
    # get features in correct order from correlations and reorder LefSe results
    features = correlations.index
    lefse_results_all = lefse_results_all.loc[features]
    lefse_results_all.columns = pd.MultiIndex.from_tuples(
        [('LefSe Results (All sequences)', x[0], x[1]) for x in lefse_results_all.columns]
    )
    lefse_results_prok = lefse_results_prok.loc[features]
    lefse_results_prok.columns = pd.MultiIndex.from_tuples(
        [('LefSe Results (Prokaryote sequences)', x[0], x[1]) for x in lefse_results_prok.columns]
    )
    
    # merge lefse results again
    merged_results = pd.concat([lefse_results_all, lefse_results_prok], axis=1)
    
    # merge tables
    merged = pd.concat([biomaker_sample_types, merged_results, correlations], axis=1)
    merged = merged.reset_index()

    data_set = merged[[('Geochemistry Correlations (Spearman\'s rho)', '', 'DataSet')]]
    data_set.columns = pd.MultiIndex.from_tuples([('', '', x[-1]) for x in data_set.columns])

    merged = merged.drop(('Geochemistry Correlations (Spearman\'s rho)', '', 'DataSet'), axis=1)
    merged = pd.concat([data_set, merged], axis=1)
    merged.index = merged[('Feature')]
    merged = merged.drop(('Feature', ''), axis=1)
    merged.index.name = label

      # save to results container
    lefse_res_corrs[label] = merged

# save final output
results_file = os.path.join(results_dir, 'lefse_results.csv')

with open(results_file, 'w') as out_handle:
    
    # write out header text
    out_handle.write('\"Lefse results and correlations with geochemical parameters for Biomarkers of End Member Brine and Non Outflow samples in the full sequence data set and prokaryote only data set at each taxonomic level.\"\n')
    
    # save first result data frame in full, then skip the headers for the rest
    lefse_res_corrs['Phylum'].to_csv(out_handle)
    
    for label in lefse_res_corrs.keys():
        if label != 'Phylum':
            
            out_handle.write('%s\n' % label)
            df = lefse_res_corrs[label]
            df.astype(str).to_csv(out_handle, header=None)
            
    out_handle.write('\n')

#========== END ==========#