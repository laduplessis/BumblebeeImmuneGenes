## Summary

Scripts used for the paper "A depauperate immune repertoire precedes evolution of sociality in bees" (Barribeau, Sadd, du Plessis and many more), Genome Biology (2015). 

Use these scripts at your own risk! THey may need some tweaks or fixes for different gene sets or organisms. No documentation is available for these scripts and I don't plan on writing. They are effectively unsupported, but if you have questions please get in touch with me and I may still be able to help.  


## Example workflows 
Example workflows for the 4-taxa tree presented in the paper
	- workflow4final.sh
	- workflow4EW.sh



## Python scripts
Pre-processing, extracting the coding sequences and translations, setting up control files, post-processing and extracting useful information from output files.

- RestrictAlignments.py:
	Extract all alignments that contain sequences from all the species in the organism list. Only sequences from species in the list are extracted. Restricted alignments are saved in a new directory in fasta format.

- ExtractDNA.py:
	Extract cDNA sequences from OrthoDB alignments
	 
	- For now the filenames containing the cDNA are hardcoded
	- Translations are tested
	- Only matched if protein sequence matches the alignment sequence and the translation
	- Sometimes cDNA sequences need to be shifted since the cDNA contains untranslated regions - this is done automatically and the untranslated ends are removed
	- In some cases there is a detected gap between 2 exons, labelled by an X in the protein sequence.  The corresponding codon in the cDNA sequence usually translates to a stop codon.  This is also detected and the residue in question is removed from the sequence.
	
	Output one file with all the protein sequences (gapless, no stopcodon)
	Output one file with all the cDNA sequences (gapless, no stopcodon)
	Output one file with statistics of all the processed alignments (Alignments.csv)
	Output one file containing all sequences with translation errors (Errors.csv)
	Output one file with the complete identifiers and sources of all the processed sequences

- GetNewAlignments.py:
	Produce alignments that can be used with PAML

    - Extract sequences for each group from the files created with ExtractDNA.py
    - Skip groups that have missing sequences for some species
    - Check translation again with OrthoDB alignments
    - If there is more than one sequence of a species:

        * Perform a new alignment on the translated sequences (AA) with the chosen MSA program for 
          every possible combination of sequences so that the MSA has only one sequence per species
        * Use Gblocks (with strict settings) to restrict the alignments to the best matching areas only
        * Select the alignment with the longest well-aligned parts (as determined by Gblocks)
        * If there are more alignments with disjoint sets of sequences select alignments by longest well aligned parts until 
          the maximum number is selected

        Else:
            * Perform a new alignment on the translated sequences with the chosen MSA program
    - Trim the alignment:
            * Using Gblocks
    - Skip groups where the trimmed alignment is shorter than the minimum length
    - Save trimmed cDNA alignments in Phylip format so it can be used with PAML

	Need to have MSA programs and Gblocks installled.  Check that programs are in your path.  
	(For prographmsa the path is hardcoded for the moment)
	(Some problems with Muscle and Prank)

	Output:

    - Alignments.csv:  Exactly like Alignments.csv from ExtractDNA.py, but only for groups selected here
    - Skipped.csv:     Groups which have been skipped for whatever reason
    - Colsremoved.csv: List of columns removed from each of the alignments

    For each group that wasn't skipped:

    - \<group\>.aa.fa:                          Protein sequences without gaps
    - \<group\>.dna.fa:                         Corresponding cDNA sequences
    - \<group\>.\<msaprogram\>.fas:               New alignment (AA, ids remapped)
    - \<group\>.\<msaprogram\>.dna.fas            New alignment (cDNA, ids remapped)
    - \<group\>.\<msaprogram\>.dna.fas-gb         Trimmed alignment from Gblocks
    - \<group\>.\<msaprogram\>.dna.fas-gb.txt     Gblocks statistics
    - \<group\>.\<msaprogram\>.dna.phylip         Trimmed cDNA alignment in Phylip, ready for Paml
    - \<group\>.\<msaprogram\>.phylip             Trimmed AA alignment in phylip
    - \<group\>.map                             Mapping of ids to the Phylip files
    - \<group\>.ids                             Identifiers of the original sequences used
    - \<group\>.combinations.tar                Temporary files from finding the best combination of sequences to use




- RunPhyml.py: 
	Run PhyML to get branch length estimates (topology not optimized, uses JTT model)

- SetupPAML.py:
	Setup control files for different PAML models (for all alignments and trees in the specified directory). Requires template control files. Outputs bash scripts for running the analyses in PAML (both locally or on the Brutus cluster). Where applicable make multiple control files for different random initializations of omega.

	- Sites: M0, M1, M2, M3, M7, M8, M8A
	- Branch-Site
	- Clade models C and D

	Trees for branch-site and clade models are hardcoded for 4- and 5-taxa trees in the paper (but can be easily changed).


- ExtractParameter.py:
	Extract parameters from the PAML output files into easily readable csv files.

- ExtractBest.py:
	Looks at output files from ExtractParameters.py to extract the replicate with the maximum-likelihood for each model.

- ExtractBEB.py:
	Extract BEB posterior values from PAML branch-site and M8 models.

- FormatBEB.py:
	Format BEB into BEB values for each of the sequences in the alignment. (So that the posterior estimate for the rate of evolution at each site in each full sequence (outside of the alignment) is known).

- GetInterPro.py:
	Download domains returned by InterPro Scan for AA sequences. A lot of functions based on http://www.ebi.ac.uk/Tools/webservices/download_lients/python/suds/iprscan5_suds.py

- MergeBEBInterPro.py:
	Merge the output from the BEB and Interpro scan into one file that can be plotted to see where sites under selection fall within domains.



## Data
Files needed as input for some scripts

- PAML templates
- Files containing sets of taxa for various analyses 
- Files containing topologies for various analyses  

## R scripts
Likelihood-ratio tests, figures, tables, Venn diagrams etc.
These need to be checked.
