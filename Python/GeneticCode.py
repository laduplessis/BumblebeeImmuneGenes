# -*- coding: utf-8 -*-
import sys




AACodes   = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',  \
             'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',  \
             'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',  \
             'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}

AAList = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


StartCodonsBacteria = ['TTG','GTG','ATT','ATC','ATA','ATG','GTG']
StartCodonsMammals  = ['ATG']

StopCodonsBacteria  = ['TAA','TAG','TGA']
StopCodonsMammals   = ['TAA','TAG','TGA']

GeneticCode	    =  {'TTT':'F', \
			'TTC':'F', \
			'TTA':'L', \
			'TTG':'L', \
			'TCT':'S', \
			'TCC':'S', \
			'TCA':'S', \
			'TCG':'S', \
			'TAT':'Y', \
			'TAC':'Y', \
			'TAA':'*', \
			'TAG':'*', \
			'TGT':'C', \
			'TGC':'C', \
			'TGA':'*', \
			'TGG':'W', \
			'CTT':'L', \
			'CTC':'L', \
			'CTA':'L', \
			'CTG':'L', \
			'CCT':'P', \
			'CCC':'P', \
			'CCA':'P', \
			'CCG':'P', \
			'CAT':'H', \
			'CAC':'H', \
			'CAA':'Q', \
			'CAG':'Q', \
			'CGT':'R', \
			'CGC':'R', \
			'CGA':'R', \
			'CGG':'R', \
			'ATT':'I', \
			'ATC':'I', \
			'ATA':'I', \
			'ATG':'M', \
			'ACT':'T', \
			'ACC':'T', \
			'ACA':'T', \
			'ACG':'T', \
			'AAT':'N', \
			'AAC':'N', \
			'AAA':'K', \
			'AAG':'K', \
			'AGT':'S', \
			'AGC':'S', \
			'AGA':'R', \
			'AGG':'R', \
			'GTT':'V', \
			'GTC':'V', \
			'GTA':'V', \
			'GTG':'V', \
			'GCT':'A', \
			'GCC':'A', \
			'GCA':'A', \
			'GCG':'A', \
			'GAT':'D', \
			'GAC':'D', \
			'GAA':'E', \
			'GAG':'E', \
			'GGT':'G', \
			'GGC':'G', \
			'GGA':'G', \
			'GGG':'G'}
			
			
AAInfo  = {'L' : {'AA'    : 'Leu', \
                  'Codons': ['TTA','TTG','CTT','CTC','CTA','CTG'], \
                  'Hydro' : 3.8}, \
	  'R' : {'AA'    : 'Arg', \
                  'Codons': ['CGT','CGC','CGA','CGG','AGA','AGG'], \
                  'Hydro' : -4.5}, \
	  'S' : {'AA'    : 'Ser', \
                  'Codons': ['TCT','TCC','TCA','TCG','AGT','AGC'], \
                  'Hydro' : -0.8}, \
          'A' : {'AA'    : 'Ala', \
                  'Codons': ['GCT','GCC','GCA','GCG'], \
                  'Hydro' : 1.8}, \
	  'G' : {'AA'    : 'Gly', \
                  'Codons': ['GGT','GGC','GGA','GGG'], \
                  'Hydro' : -0.4}, \
	  'P' : {'AA'    : 'Pro', \
                  'Codons': ['CCT','CCC','CCA','CCG'], \
                  'Hydro' : 1.6}, \
	  'T' : {'AA'    : 'Thr', \
                  'Codons': ['ACT','ACC','ACA','ACG'], \
                  'Hydro' : -0.7}, \
	  'V' : {'AA'    : 'Val', \
                  'Codons': ['GTT','GTC','GTA','GTG'], \
                  'Hydro' : 4.2}, \
	  'I' : {'AA'    : 'Ile', \
                  'Codons': ['ATT','ATC','ATA'], \
                  'Hydro' : 4.5}, \
	  'C' : {'AA'    : 'Cys', \
                  'Codons': ['TGT','TGC'], \
                  'Hydro' : 2.5}, \
	  'D' : {'AA'    : 'Asp', \
                  'Codons': ['GAT','GAC'], \
                  'Hydro' : -3.5}, \
	  'E' : {'AA'    : 'Glu', \
                  'Codons': ['GAA','GAG'], \
                  'Hydro' : -3.5}, \
	  'F' : {'AA'    : 'Phe', \
                  'Codons': ['TTT','TTC'], \
                  'Hydro' : 2.8}, \
	  'H' : {'AA'    : 'His', \
                  'Codons': ['CAT','CAC'], \
                  'Hydro' : -3.2}, \
	  'K' : {'AA'    : 'Lys', \
                  'Codons': ['AAA','AAG'], \
                  'Hydro' : -4.5}, \
	  'N' : {'AA'    : 'Asn', \
                  'Codons': ['AAT','AAC'], \
                  'Hydro' : -3.5}, \
	  'Q' : {'AA'    : 'Gln', \
                  'Codons': ['CAA','CAG'], \
                  'Hydro' : -3.5}, \
	  'Y' : {'AA'    : 'Tyr', \
                  'Codons': ['TAT','TAC'], \
                  'Hydro' : -1.3}, \
	  'M' : {'AA'    : 'Met', \
                  'Codons': ['ATG'], \
                  'Hydro' : 1.9}, \
	  'W' : {'AA'    : 'Trp', \
                  'Codons': ['TGG'], \
                  'Hydro' : -0.9},
	  '*' : {'AA'    : '***', \
                  'Codons': ['TAA','TAG','TGA'], \
                  'Hydro' : 0}}



def idxToCodon(n):
	map = {0:'T',1:'C',2:'A',3:'G'}
		
	if n < 0: raise ValueError, "must be a positive integer"
	
	if n == 0: return 'TTT'
	cod = ''
	while n > 0:
		cod = map[n % 4] + cod
		n = n >> 2
	#
	if (len(cod) < 3):
		cod = 'T'*(3-len(cod))+cod
	return cod
		
	
def codonToIdx(cod):
	nuc = 'ACGT'
	map = {'T':'0','C':'1','A':'2','G':'3'}
	
	cod = cod.upper().replace('U','T')
	
	if (cod[0] in nuc and cod[1] in nuc and cod[2] in nuc):
		return int(map[cod[0]]+map[cod[1]]+map[cod[2]],4)
	else:
		return -1
#	


def getGCMap():
	GCMap = [[],[],[]]
	for i in range(0,64):
		for j in range(0,3):
			if (idxToCodon(i)[j] == 'G' or idxToCodon(i)[j] == 'C'):
				GCMap[j].append(1)
			else:
				GCMap[j].append(0)
	return GCMap


	

def translateCodon(codon):
	"""Translates a codon according the standard/bacterial genetic code
	
	(Bacterial and Standard Genetic code have the same translations)
	Only checks if codon length is not equal to 3.
	Works on RNA and DNA sequences.
	Returns X for all ambiguous codons, - for a codon that is all gaps.
	"""
	
	nuc = 'ACGT'
	codon = codon.upper().replace('U','T')
	
	if len(codon) != 3:
		sys.stderr.write('Error! Codon of wrong length!\n')
		raise Exception
		
	if (codon[0] in nuc and codon[1] in nuc and codon[2] in nuc):
		return GeneticCode[codon]
	elif (codon == '---'):
		return '-'
	else:
		return 'X'
##

def parseCodons(sequence):
	"""Iterable that returns the codons in a sequence
	
	No error checking!
	"""
	for i in range(0,len(sequence),3):
		yield sequence[i:i+3]
##

def parseCodonstoAA(sequence):
	"""Iterable that returns the codons in a sequence, translated to AA
	
	No error checking!
	"""
	for i in range(0,len(sequence),3):
		yield translateCodon(sequence[i:i+3])
##	


def translateSeq(seq):
	"""Translate a sequence to AA, no error checking done
	
	"""
	
	return "".join(list(parseCodonstoAA(seq))) 
##

def translateCDS(sequence,table,includestop = False):
	"""Translates a coding sequence
	
	Will check if the sequence starts with a valid start codon and ends with a stop codon
	Will return an error if the sequence is not at least 2 codons in length
	If includestop is set the returned sequence will end with *
	"""
	
	if (len(sequence) % 3 != 0 and len(sequence) < 6):
		sys.stderr.write("Error! Sequence length not divisible by 3 or too short (need at least 2 codons)!\n")
		raise
	#
	
	start = sequence[:3]
	stop  = sequence[-3:]
	
	translation = translateSeq(sequence[3:-3])
	if (translation.find('*') < 0):
		if (includestop):
			translation += '*'
	
		if   (table == 'Standard'):
			if  (start in StartCodonsMammals and stop in StopCodonsMammals):
				return 'M'+translation 	
						
		elif (table == 'Bacterial'):
			if  (start in StartCodonsBacteria and stop in StopCodonsBacteria):
				return 'M'+translation
		else:
			sys.stderr.write('Error! Unknown genetic code!\n')
			raise
	#
	
	sys.stderr.write("Error! Not a valid coding sequence!\n")
	raise
##






