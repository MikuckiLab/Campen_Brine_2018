# High resolution characterization of an Antarctic groundwater microbiome: replicate sampling and delineation of subglacial outflow at Blood Falls 

R. Campen, W.B. Lyons, S. Tulaczyk, B. Dachwald, J. Kowalski, E.C. Pettit, K.A. Welch, J.A. Mikucki

---

### Abstract

Antarctic subglacial environments host microbial ecosystems and are proving to be geochemically and biologically diverse. The Taylor Glacier, Antarctica, periodically expels iron-rich brine through an englacial conduit sourced from a deep subglacial aquifer. The discharged groundwater creates a dramatic red feature on the glacier surface known as Blood Falls. Here we describe the microbial diversity of subglacial fluids collected from this feature over a decade-long time frame. This represents the first replicate sampling of an Antarctic subglacial system, and constitutes the first collection from within an englacial conduit using a novel, minimally invasive thermoelectric melting drill. Our comparative analysis allows for the delineation of a subglacial brine core microbiome that is stable. Previously undetected but abundant groups including the candidate bacterial phylum Atribacteria and the archaeal phylum Pacearchaeota. Many of the core taxa correlated strongly with brine geochemistry. Comparison with samples collected at the glacier surface when discharge was not active and from within its proglacial lake; this enabled the delineation of subsurface samples and provides insight into brine hydrology. The ability to evaluate sample purity will enable future investigations of surface collection of subglacial brine. Our work has implications for deconvoluting subsurface contributions to the composition of brines released from thick ice covers, such as may be encountered on Europa or Enceladus.

---

### Requirements

It is recommended to use the [Anaconda Python3](https://www.anaconda.com/download/) environment to ensure the majority of dependencies for running the code in this repository are met. To install the remaining dependencies run `pip install mothur-py seq-experiment`.

You will also require [mothur](https://github.com/mothur/mothur) (See below for version compatibility issues), and the Silva reference alignment which can be downloaded using the script `0_download_silva.py`. Before doing so ensure you have read to, and agreed to the Silva [Terms of Use/License Agreement](https://www.arb-silva.de/silva-license-information).

You will also need to download the raw sequence files which the NCBI Short Read Archive (SRA) under accession number `<insert_accession_numbers>`, and place them in the `data/raw_seq/` folder.

---

### Repository layout

```
this_is_outflow
|- README.md   			# Top level readme file (this file)
|- .gitignore  				# configuration for git of which files to ignore
|
|- data/				# raw and intermediate processed files
| |- 515806.oligos 			# oligos file for running mothur's make.contigs command
| |- HMP_MOCK.v35.fasta 		# sequences assosiated with the Mock Community sample used for seq.error
| |- sample_metadata.csv 		# table of sample metadata
| |+ raw_seqs/				# folder containing raw sequence files
| |+ processed/				# folder containing intermediate processed files
| |+ silva/				# folder containing downloaded Silva reference files (if downloaded)
| |+ tmp/				# temporary files created during analysis that can be deleted safely
|
|+ results/				# final results tables output from analysis
|
|- code/				# folder contianing scripts for sequence processing and data analysis
| |- 0_download_silva.py  		# download Silva v128 reference files (SEE BELOW)
| |- 1_process_sequences.py 		# run mothur sequence processing pipeline (SEE BELOW)
| |- 2_correlate_geochemistry.py 	# correlate geochemical parameters
| |- 3_estimate_alpha_diversities.py 	# estimate alpha diversities and correlate with geochemistry
| |- 4_EMB_core_microbiome.py 		# determine the EMB core microbiome OTU composition
| |- 5_identify_biomarkers.py 		# identify biomarkers using LefSe and correlate with geochemistry
```

legend:
* `|-` denotes file
* `|- <directory_name>/` denotes directory whose contents are listed
* `|+ <directory_name>/` denotes directory whose contents are not listed
* `| |-` denotes file/directory is within another directory


**CAUTION:** Due to a bug in mothur v1.39.5 (and earlier) the outputs will vary slightly between each run of any mothur commands. This should be resolved in the next release version, or can be avoided by using a build of the `Threads_373` branch of mothur from commit [`1e8aa08`](https://github.com/mothur/mothur/commit/1e8aa085dc33d2d874b9819fc869d6b000eb2ab7) onwards where this bug was fixed. Re-running the files as denoted with `SEE BELOW` in the above file structure with a version of mothur with this bug will produce different final results tables/figures. It is therefore advised to not rerun these files unless absolutely necessary, or without a compatible version of mothur. The outputs of these files are already in the `data/processed/` directory if needed. All other files have fully reproducibly outputs and can be rerun as much as desired.
