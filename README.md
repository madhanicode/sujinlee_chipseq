# sujinlee_chipseq

Description:
This repository contains a loose collection of scripts used to analyze and generate graphs for ChIP-Seq data.
These scripts are customized for specific experiments and may not be generally applicable.
Scripts should not be used directly, and should serve primarily as reference.

Requirements:
Required packages are listed in dani_pipeline.yml
Python scripts use Python 3.9.9

Usage:
Starting from .bedgraph files output from Jordan Burke's ChIP-Seq pipeline
1. run get_background.sh
2. run WCE_bedgr_to_cen_tel.sh
3. Run either dani_analysis.py or dani_analysis_additional.py to generate desired graphs

Notes:
get_background.sh
This script calls bedtools to map all reads not in centromeric or telomeric regions, then uses those values to calculate an average background read depth per chromosome.
Change refdir to directory containing .bed file "NOT_CEN_TEL_labeled.bed"
get_background.sh bedgraph_dir output_dir

Authorship:
Scripts were originally written by Daniela Perry, then proofread and modified by Manning Huang. 
