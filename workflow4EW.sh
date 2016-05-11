#!/bin/bash

#####################################################################################################################################################
# 1) Restrict alignment to a set of species (Extract all alignments that have one or more of any of the species)
python RestrictAlignments.py -i ../Data/EIGER_WENGEN/ -l ../Data/Taxa4.txt -o ../Data/4SpeciesEW/Restricted/

# 2) Restrict missing alignments (10 groups) manually and copy to Data/4Species/Restricted/
#####################################################################################################################################################

#####################################################################################################################################################
# 3) Extract DNA and AA sequences from restricted alignments (DNA and AA into two big files) 
# 		THIS STEP TAKES SOME TIME! ~30 minutes!
python ExtractDNA.py -i ../Data/4SpeciesEW/Restricted/ -g ../Data/Genes/ -o ../Data/4SpeciesEW/



#####################################################################################################################################################
# 4) Extract new trimmed alignments
#		Can choose to only extract alignments with certain numbers of each species
#		Alignmnts trimmed here (either with Gblocks or by simply removing ambiguous columns)
#		Rename sequences and save alignments in Phylip format
#		If the maximum number of sequences per species in a group is more than one this will take a LONG TIME (~2 HOURS for ProGraphMSA, ~5 HOURS for Probcons)

# ProGraphMSA Strict
python GetNewAlignments.py -i ../Data/4SpeciesEW/ -a ../Data/4SpeciesEW/Restricted/ -o ../Data/4SpeciesEW/MSA_ProGraphMSA_Strict/ -s 1,1,0,0,1,1 -m prographmsa -g strict -l 10

# Probcons Strict
#python GetNewAlignments.py -i ../Data/4Species/ -a ../Data/4Species/Restricted/ -o ../Data/4Species/MSA_Probcons_Strict/ -s 1,1,0,0,1,1 -m probcons -g strict -l 10 -c ../Data/4Species/MSA_ProGraphMSA_Strict/



#####################################################################################################################################################
# 5) Get trees (based on AA sequences)

# ProGraphMSA Strict
python RunPhyml.py -i ../Data/4SpeciesEW/MSA_ProGraphMSA_Strict/ -m prographmsa -p phyml -t ../Data/InputTree4.txt 

# Probcons Strict
#python RunPhyml.py -i ../Data/4Species/MSA_Probcons_Strict/ -m probcons -p phyml -t ../Data/InputTree4.txt 




#####################################################################################################################################################
# 6) Make PAML control files

# ProGraphMSA Strict
python SetupPAML.py -i ../Data/4SpeciesEW/MSA_ProGraphMSA_Strict/ -o ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict/ -c ../Data/PAML_Templates/ -m prographmsa -s 4 -R 3 -r 25

# Probcons Strict
#python SetupPAML.py -i ../Data/4Species/MSA_Probcons_Strict/ -o ../Data/4Species/PAML_Probcons_Strict/ -c ../Data/PAML_Templates/ -m probcons -s 4 -R 3 -r 25


#####################################################################################################################################################
# 7) Run PAML

#### Done on Brutus



#####################################################################################################################################################
# 8) Extract the parameters found by PAML for each of the models

# ProGraphMSA Strict
python ExtractParameters.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict/ -o ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_Results/

# Probcons Strict
#python ExtractParameters.py -i ../Data/4Species/PAML_Probcons_Strict/ -o ../Data/4Species/PAML_Probcons_Strict_Results/



#####################################################################################################################################################
# 9) Extract the best run

# ProGraphMSA Strict
python ExtractBest.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_Results

# Probcons Strict
#python ExtractBest.py -i ../Data/4Species/PAML_Probcons_Strict_Results



#####################################################################################################################################################
# 10) Extract BEB results for M8 and Branch-site models for the best runs

# ProGraphMSA Strict
python ExtractBEB.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict -b ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_Results/ -m prographmsa

# Probcons Strict
#python ExtractBEB.py -i ../Data/4Species/PAML_Probcons_Strict -b ../Data/4Species/PAML_Probcons_Strict_Results/ -m probcons



#####################################################################################################################################################
# 11) Copy alignments summary file to results file (all results files in one directory for easy analysis)

# ProGraphMSA Strict
cp ../Data/4SpeciesEW/MSA_ProGraphMSA_Strict/Alignments.csv ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_Results/
cat 4Species/MSA_ProGraphMSA_Strict/Alignments.csv 4SpeciesEW/MSA_ProGraphMSA_Strict/Alignments.csv > 4SpeciesEW/PAML_ProGraphMSA_Strict_Results/Alignments.csv 


# Probcons Strict
#cp ../Data/4Species/MSA_Probcons_Strict/Alignments.csv ../Data/4Species/PAML_Probcons_Strict_Results/


python ExtractParameters.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_StrictEW/ -o ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/
python ExtractBest.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW
python ExtractBEB.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict -b ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/ -m prographmsa



python FormatBEB.py -i ../Data/4SpeciesEW/MSA_ProGraphMSA_Strict/ -b ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/BEB/ -m prographmsa -o ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/InterPro/
python GetInterPro.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/InterPro/ -o ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_ResultsEW/InterPro/Raw/ -e louis.duplessis@env.ethz.ch -g
python MergeBEBInterPro.py -i ../Data/4SpeciesEW/PAML_ProGraphMSA_Strict_Results/InterPro/ -d Smart,pfam,phobius

