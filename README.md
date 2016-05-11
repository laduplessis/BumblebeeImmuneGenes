## Summary

Scripts used for the paper "A depauperate immune repertoire precedes evolution of sociality in bees" (Barribeau, Sadd, du Plessis and many more), Genome Biology (2015). 

Use these scripts at your own risk! THey may need some tweaks or fixes for different gene sets or organisms. No documentation is available for these scripts and I don't plan on writing. They are effectively unsupported, but if you have questions please get in touch with me and I may still be able to help.  


## Example workflows 
Example workflows for the 4-taxa tree presented in the paper
	- workflow4final.sh
	- workflow4EW.sh



## Python scripts
Pre-processing, extracting the coding sequences and translations, setting up control files, post-processing and extracting useful information from output files.

- RestrictAlignments.py
- ExtractDNA.py
- GetNewAlignments.py
- RunPhyml.py: 
	Run PhyML to get branch length estimates (topology not optimized, uses JTT model)

- SetupPAML.py
	Setup control files for different PAML models (for all alignments and trees in the specified directory). Requires template control files. Outputs bash scripts for running the analyses in PAML (both locally or on the Brutus cluster)

	- Sites: M0, M1, M2, M3, M7, M8, M8A
	- Branch-Site
	- Clade models C and D

	Trees for branch-site and clade models are hardcoded for 4- and 5-taxa trees in the paper (but can be easily changed).


- ExtractParameter.py
- ExtractBest.py
- ExtractBEB.py
- FormatBEB.py
- GetInterPro.py
- MergeBEBInterPro.py


## R scripts
Likelihood-ratio tests, figures, tables, Venn diagrams etc.
