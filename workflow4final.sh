#!/bin/bash

#####################################################################################################################################################
# 1) Restrict alignment to a set of species (Extract all alignments that have one or more of any of the species)
python RestrictAlignments.py -i ../DivergenceData/May-Immunity/ -l ../DivergenceData/Taxa4.txt -o ../DivergenceData/4Species/Restricted/



#####################################################################################################################################################
# 2) Extract DNA and AA sequences from restricted alignments (DNA and AA into two big files) 
# 		THIS STEP TAKES SOME TIME! ~30 minutes!
python ExtractDNA.py -i ../DivergenceData/4Species/Restricted/ -g ../DivergenceData/Genes/ -o ../DivergenceData/4Species/



#####################################################################################################################################################
# 3) Extract new trimmed alignments
#		Can choose to only extract alignments with certain numbers of each species
#		Alignmnts trimmed here (either with Gblocks or by simply removing ambiguous columns)
#		Rename sequences and save alignments in Phylip format
#		If the maximum number of sequences per species in a group is more than one this will take a LONG TIME (~2 HOURS for ProGraphMSA, ~5 HOURS for Probcons)

# ProGraphMSA Strict
python GetNewAlignments.py -i ../DivergenceData/4Species/ -a ../DivergenceData/4Species/Restricted/ -o ../DivergenceData/4Species/MSA_ProGraphMSA_Strict/ -s 1,1,0,0,1,1 -m prographmsa -g strict -l 50

# ProGraphMSA Relaxed
python GetNewAlignments.py -i ../DivergenceData/4Species/ -a ../DivergenceData/4Species/Restricted/ -o ../DivergenceData/4Species/MSA_ProGraphMSA_Relaxed/ -s 1,1,0,0,1,1 -m prographmsa -g relaxed -l 50

# ProGraphMSA No trimming
python GetNewAlignments.py -i ../DivergenceData/4Species/ -a ../DivergenceData/4Species/Restricted/ -o ../DivergenceData/4Species/MSA_ProGraphMSA_None/ -s 1,1,0,0,1,1 -m prographmsa -p 0 -l 50

# Probcons Strict
python GetNewAlignments.py -i ../DivergenceData/4Species/ -a ../DivergenceData/4Species/Restricted/ -o ../DivergenceData/4Species/MSA_Probcons_Strict/ -s 1,1,0,0,1,1 -m probcons -g strict -l 50


#####################################################################################################################################################
# 4) Get trees (based on AA sequences)

# ProGraphMSA Strict
python RunPhyml.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_Strict/ -m prographmsa -p phyml -t ../DivergenceData/InputTree4.txt 

# ProGraphMSA Relaxed
python RunPhyml.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_Relaxed/ -m prographmsa -p phyml -t ../DivergenceData/InputTree4.txt 

# ProGraphMSA No trimming
python RunPhyml.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_None/ -m prographmsa -p phyml -t ../DivergenceData/InputTree4.txt 

# Probcons Strict
python RunPhyml.py -i ../DivergenceData/4Species/MSA_Probcons_Strict/ -m probcons -p phyml -t ../DivergenceData/InputTree4.txt 



#####################################################################################################################################################
# 5) Make PAML control files

# ProGraphMSA Strict
python SetupPAML.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_Strict/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_Strict/ -c ../DivergenceData/PAML_Templates/ -m prographmsa -s 4 -R 3 -r 25

# ProGraphMSA Relaxed
python SetupPAML.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_Relaxed/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed/ -c ../DivergenceData/PAML_Templates/ -m prographmsa -s 4 -R 3 -r 25

# ProGraphMSA No trimming 
python SetupPAML.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_None/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_None/ -c ../DivergenceData/PAML_Templates/ -m prographmsa -s 4 -R 3 -r 25

# Probcons Strict
python SetupPAML.py -i ../DivergenceData/4Species/MSA_Probcons_Strict/ -o ../DivergenceData/4Species/PAML_Probcons_Strict/ -c ../DivergenceData/PAML_Templates/ -m probcons -s 4 -R 3 -r 25


#####################################################################################################################################################
# 6) Run PAML

#### Done on Brutus



#####################################################################################################################################################
# 7) Extract the parameters found by PAML for each of the models

# ProGraphMSA Strict
python ExtractParameters.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Strict/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/

# ProGraphMSA Relaxed
python ExtractParameters.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed_Results/

# ProGraphMSA None
python ExtractParameters.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_None/ -o ../DivergenceData/4Species/PAML_ProGraphMSA_None_Results/


#####################################################################################################################################################
# 8) Extract the best run

# ProGraphMSA Strict
python ExtractBest.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results

# ProGraphMSA Relaxed
python ExtractBest.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed_Results

# ProGraphMSA None
python ExtractBest.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_None_Results


#####################################################################################################################################################
# 9) Extract BEB results for M8 and Branch-site models for the best runs

# ProGraphMSA Strict
python ExtractBEB.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Strict -b ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/ -m prographmsa

# ProGraphMSA Relaxed
python ExtractBEB.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed -b ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed_Results/ -m prographmsa

# ProGraphMSA None
python ExtractBEB.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_None -b ../DivergenceData/4Species/PAML_ProGraphMSA_None_Results/ -m prographmsa


#####################################################################################################################################################
# 10) Copy alignments summary file to results file (all results files in one directory for easy analysis)

# ProGraphMSA Strict
cp ../DivergenceData/4Species/MSA_ProGraphMSA_Strict/Alignments.csv ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/

# ProGraphMSA Relaxed
cp ../DivergenceData/4Species/MSA_ProGraphMSA_Relaxed/Alignments.csv ../DivergenceData/4Species/PAML_ProGraphMSA_Relaxed_Results/

# ProGraphMSA None
cp ../DivergenceData/4Species/MSA_ProGraphMSA_None/Alignments.csv ../DivergenceData/4Species/PAML_ProGraphMSA_None_Results/








python FormatBEB.py -i ../DivergenceData/4Species/MSA_ProGraphMSA_Strict/ -b ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/BEB/ -m prographmsa -o ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/InterPro/

python MergeBEBInterPro.py -i ../DivergenceData/4Species/PAML_ProGraphMSA_Strict_Results/InterPro/ -d smart,pfam,phobius


