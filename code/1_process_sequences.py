# Campen_Brine_2018
# 1_process_sequences

"""
This file deals with the sequencing processing of samples for characterising
Blood Falls outflow brine in comparison to non-outflow samples collected 
on/near Blood Falls, and samples from the deeper water column of West lobe
Lake Bonney. 

This file describes the optimised sequence processing pipeline.

To run this analysis you will require:

* Mothur executable
* Uchime executable (included with mothur)
* Raw fastq sequence files

Additionally you will need the following, and to specify their location in 
the mothur_vars dictionary below.

* Mothur compatible SILVA reference alignment
* Mock community referece fasta file

This notebook was run with Mothur version 1.39.5.

"""

#========== General Setup ==========#

# import necessary libraries
import os
import shutil

import numpy as np
import pandas as pd
import seq_experiment as se

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from mothur_py import Mothur


# configure directories
base_dir = os.path.join('..')
base_dir = os.path.abspath(base_dir)
data_dir = os.path.join(base_dir, 'data')

processed_dir = os.path.join(data_dir, 'processed')
results_dir = os.path.join(base_dir, 'results')

tmp_dir =  os.path.join(data_dir, 'tmp')
if not os.path.isdir(tmp_dir):
	os.mkdir(tmp_dir)

mothur_out_dir = os.path.join(tmp_dir, '1_process_sequences')
if not os.path.isdir(mothur_out_dir):
	os.mkdir(mothur_out_dir)

#========== Mothur Setup ==========#

# flag whether to assess error rate or not
# NOTE:: if True, this version of mothur will crash. This does not prevent
# calculation of errors, it simply requires the user to monitor the 
# analysis and ensure it does not stall
ASSESS_ERRORS = False

# create Mothur instance and set configuration
m = Mothur()
m.verbosity = 1
m.suppress_logfile = True
m.mothur_seed = 4321
m.current_dirs['output'] = m.current_dirs['tempdefault'] = mothur_out_dir

# setup other variables for mothur execution
mothur_vars = {
    'max_processors': 1,  # the more processors the more memory used
    # for make.contigs
    'pdif': 3,
    'oligos_file': os.path.join(data_dir, '515806.oligos'),
    # for pcr.seqs()
    'max_homopolymer': 12,
    'max_ambiguous': 0,
    # for screen.seqs()
    'screen_start': 8,
    'screen_end': 9582,
    # for precluster
    'preclust_difs': 2,
    # for chimera.uchime
    'chimera_derep': True,
    # for seq.error()
    'mock_group': 'Mock',
    'ref_mock': os.path.join(data_dir, 'HMP_MOCK.v35.fasta'),
    'mock_aligned': False,
    # for classify.seqs()
    'classify_cutoff': 80,
    # for dist.seqs(); also used in downstream commands
    'otu_level': 0.01,
    # for cluster.seqs()
    'clust_method': 'opti',
}

# specify folder containing the SILVA reference files and add to mothur_vars
silva_dir = os.path.join(data_dir, 'silva', 'silva_v128')
mothur_vars['ref_align'] = os.path.join(silva_dir, 'silva.nr_v128.V4.align')
mothur_vars['ref_tax'] = os.path.join(silva_dir, 'silva.nr_v128.tax')

#---------- Create .files file ----------#

# setup input dir
raw_seq_dir = os.path.join(data_dir, 'raw_seqs')

# get seq_ids and matching sample names from metadata file
metadata = pd.read_csv(
	os.path.join(data_dir, 'sample_metadata.csv'), 
	index_col=0
)
seq_files_df = pd.DataFrame(metadata.loc[:, 'seq_id'])

# add in seq files, using seq_id entries in table to generate correct file names
for read in ['R1', 'R2']:
    for sample in seq_files_df.index:
        seq_id = seq_files_df.get_value(sample, 'seq_id')
        seq_files_df.set_value(sample, read, os.path.join(raw_seq_dir, '%s-%s.fastq' % (seq_id, read)))
                
# tidy up df and write out to file
seq_files_df = seq_files_df.drop('seq_id', axis=1).reset_index()
files_file = os.path.join(mothur_out_dir, 'TIO.files')
seq_files_df.to_csv(files_file, sep='\t' , header=False , index=False)

#========== Sequence Processing ==========#

# make contigs from forward and reverse reads
m.make.contigs(
	file=files_file, 
	processors=mothur_vars['max_processors'], 
	oligos=mothur_vars['oligos_file'], 
    pdiffs=mothur_vars['pdif']
)

# remove sequences with too many homopolymers or ambiguities
m.screen.seqs(
	fasta='current', 
	group='current', 
	maxambig=mothur_vars['max_ambiguous'], 
	maxhomop=mothur_vars['max_homopolymer']
)

# get unqiue sequences to reduce computational load
m.unique.seqs()

# count number of sequences per sample
m.count.seqs(name='current', group='current')

# align sequence to reference alignment
m.align.seqs(fasta='current', reference=mothur_vars['ref_align'])

# remove aligned sequences not fully overlapping region of interest 
m.screen.seqs(
	fasta='current', 
	count='current', 
	summary='current', 
	start=mothur_vars['screen_start'],
	end=mothur_vars['screen_end']
)

# tidy up alignment
m.filter.seqs(fasta='current', vertical=True, trump='.')

# get unqiue sequences to reduce computational load
m.unique.seqs(fasta='current', count='current')

# denoise sequence data by preclustering
m.pre.cluster(
	fasta='current', 
	count='current', 
	diffs=mothur_vars['preclust_difs']
)

# tidy up some file names
# rename fasta file and update current files
new_fasta_file = os.path.join(m.current_dirs['output'], 'TIO.an.fasta')    
m.rename.file(input=m.current_files['fasta'], new=new_fasta_file)
m.current_files['fasta'] = new_fasta_file

# rename count file and update current files
new_count_file = os.path.join(m.current_dirs['output'], 'TIO.an.count_table ')
m.rename.file(input=m.current_files['count'], new=new_count_file)
m.current_files['count'] = new_count_file

# check for chimeric sequences
m.chimera.uchime(
	fasta='current', 
	count='current', 
	dereplicate=mothur_vars['chimera_derep'], 
	processors=1
)

# remove chimeric sequences
m.remove.seqs(fasta='current', accnos='current')

# classify sequences by reference taxonomy
m.classify.seqs(
	fasta='current', 
	count='current', 
	reference=mothur_vars['ref_align'], 
	taxonomy=mothur_vars['ref_tax'], 
    cutoff=mothur_vars['classify_cutoff'], 
    processors=mothur_vars['max_processors']
)

#---------- Assessing erros ----------#

## NOTE:: Mothur will crash when running this command, we just power through

if ASSESS_ERRORS:

	from copy import deepcopy

	# create copy of mothur object for mock community only
	m_mock = deepcopy(m)

	# set in/out dirs for easier downstream parsing of results
	m_mock.current_dirs['output'] = os.path.join('output', 'mock')
	m_mock.current_dirs['input'] = '.'

	# run mothur and save resulting object
	m_mock.get.groups(
		count='current', 
		fasta='current', 
		groups=mothur_vars['mock_group']
	)

	try:
	    m_mock.seq.error(
	    	fasta='current', 
	    	count='current', 
	    	reference=mothur_vars['ref_mock'], 
	    	aligned=mothur_vars['mock_aligned']
	    )
	except RuntimeError:
	    # mothur can error on seq.error command so need to 
	    # forcefully ignore it and keep going. If there is 
	    # a prompt requiring you to close the non-responding
	    # program, do so.
	    pass

#---------- Continue sequence processing ----------#

# calculate sequence distances
m.dist.seqs(fasta='current', cutoff=mothur_vars['otu_level'])

# cluster sequences into OTUs
m.cluster(
	column='current', 
	count='current', 
	cutoff=mothur_vars['otu_level'], 
	method=mothur_vars['clust_method']
)

# create shared file
m.make.shared(
	list='current', 
	count='current', 
	label=mothur_vars['otu_level']
)

# classify OTUs based on the consensus classification of the sequences
m.classify.otu(
	list='current', 
	count='current', 
	taxonomy='current', 
	label=mothur_vars['otu_level']
)

# get representitive fasta sequence for each OTU
m.get.oturep(
	list='current', 
	count='current', 
	fasta='current', 
	method='abundance', 
	label=mothur_vars['otu_level']
)

# format mothur results into the database file format
m.create.database(
	list='current', 
	label=mothur_vars['otu_level'], 
	repfasta=m.current_files['fasta'], 
	constaxonomy='current', 
	count='current'
)

#========== Tidy Output ==========#

shutil.copyfile(
	src=m.current_files['fasta'], 
	dst=os.path.join(processed_dir, 'TIO.rep.fasta')
)
shutil.copyfile(
	src=m.current_files['shared'], 
	dst=os.path.join(processed_dir, 'TIO.shared')
)
shutil.copyfile(
	src=m.current_files['constaxonomy'], 
	dst=os.path.join(processed_dir, 'TIO.cons.taxonomy')
)

# tidy up database file then export to csv
# read in database file
database = pd.read_table(m.output_files['database'][0])

# tidy classifications formatting
classifications = database['OTUConTaxonomy']
classifications = classifications.str.split(';', expand=True).drop(6, axis=1)

new_classifications = classifications.copy(deep=True)
for column in classifications:
    new_classifications[column] = new_classifications[column].str.rsplit('(', 1, expand=True)[0]       

classifications = new_classifications
classifications.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']

# remake database dataframe with tidier classifications formatting
database = pd.concat([database['OTUNumber'], database[seq_files_df['sample']], database['repSeq'], classifications], axis=1)
database = database.set_index('OTUNumber')

# write out database file
database.to_csv(os.path.join(processed_dir, 'TIO.database.csv'))

#========== END ==========#
