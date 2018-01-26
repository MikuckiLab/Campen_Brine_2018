### Repository layout

legend:
* `|-` denotes file
* `|- <directory_name>/` denotes directory whose contents are listed
* `|+ <directory_name>/` denotes directory whose contents are not listed
* `| |-` denotes file/directory is within another directory


	this_is_outflow
	|- README	   				# Top level readme file (this file)
	|- .gitignore  				# configuration for git of which files to ignore
	|
	|- data/					# raw and intermediate processed files
	| |- 515806.oligos 				# oligos file for running mothur's make.contigs command
	| |- HMP_MOCK.v35.fasta 		# sequences assosiated with the Mock Community sample, used for seq.error
	| |- sample_metadata.csv 		# table of sample metadata
	| |+ raw_seqs/					# folder containing raw sequence files
	| |+ processed/					# folder containing intermediate processed files
	| |+ silva/						# folder containing downloaded Silva reference files (if downloaded)
	| |+ tmp/						# temporary files created during analysis that can be deleted safely
	|
	|+ results/					# final results tables and figures output from analysis
	|
	|- code/					# folder contianing scripts for sequence processing and data analysis
	| |- 0_download_silva.py  				# download Silva v128 reference files (**SEE BELOW**)
	| |- 1_process_sequences.py 			# run mothur sequence processing pipeline (**SEE BELOW**)
	| |- 2_correlate_geochemistry.py 		# correlate geochemical parameters
	| |- 3_estimate_alpha_diversities.py 	# estimate alpha diversities and correlate with geochemistry
	| |- 4_identify_biomarkers.py 			# identify biomarkers using LefSe and correlate with geochemistry
	


**CAUTION:** Due to a bug in mothur v1.39.5 the outputs will vary slightly between each run of any mothur commands. As a result the sequence analysis will produce slightly different results each time (in terms of taxonomic classification). As a result, re-running these particular files as denoted with **SEE BELOW** in the above file structure will produce different final results tables/figures. It is therefore advised to not rerun these files unless absolutely necessary. The outputs of these files are already in the `data/processed` directory if needed. All other files have fully reproducibly outputs and can be rerun as much as desired.
