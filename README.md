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
- RunPhyml.py
	Run PhyML to get branch length estimates (topology not optimized, uses JTT model)
- SetupPAML.py
- ExtractParameter.py
- ExtractBest.py
- ExtractBEB.py
- FormatBEB.py
- GetInterPro.py
- MergeBEBInterPro.py


## R scripts
Likelihood-ratio tests, figures, tables, Venn diagrams etc.
