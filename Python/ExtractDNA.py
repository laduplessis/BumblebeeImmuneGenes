import sys, os, time
from fnmatch import fnmatch
from optparse import OptionParser
from Bio import SeqIO


# Extract cDNA sequences from OrthoDB alignments
# 
# - For now the filenames containing the cDNA are hardcoded
# - Translations are tested
# - Only matched if protein sequence matches the alignment sequence and the translation
# - Sometimes cDNA sequences need to be shifted since the cDNA contains untranslated regions - this is done
#   automatically and the untranslated ends are removed
# - In some cases there is a detected gap between 2 exons, labelled by an X in the protein sequence.  The corresponding
#   codon in the cDNA sequence usually translates to a stop codon.  This is also detected and the residue in question
#   is removed from the sequence.
#
# Output one file with all the protein sequences (gapless, no stopcodon)
# Output one file with all the cDNA sequences (gapless, no stopcodon)
# Output one file with statistics of all the processed alignments (Alignments.csv)
# Output one file containing all sequences with translation errors (Errors.csv)
# Output one file with the complete identifiers and sources of all the processed sequences
# 
################################################################################################################################
# Parameters
################################################################################################################################

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)


parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing alignments in Fasta format [default = %default]")

parser.add_option("-g","--genepath",
                  dest = "genepath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing genes [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "path to store output files in [required]")

parser.add_option("-p","--pattern",
                  dest = "pattern",
                  default = "*",
                  metavar = "filename",
                  help = "Pattern for alignment filenames to search for (in quotes) [default = %default]")

(options,args) = parser.parse_args()

inputpath  = os.path.abspath(options.inputpath)
genepath   = os.path.abspath(options.genepath)
outputpath = os.path.abspath(options.outputpath)
pattern    = options.pattern


# Filenames... hardcoded for now
amell_rna = "amel_OGSv3.2_cds.fa"
amell_pep = "amel_OGSv3.2_pep.fa"
aflor_rna = "au1.af.codingseq"
aflor_pep = "au1.af.aa"
nvitr_rna = "Nvitr_cds.fa"
nvitr_pep = "Nvitr_protein.fa"
mrotu_rna = "megachile_rotundata.all.maker.transcripts.fasta"
mrotu_pep = "megachile_rotundata.all.maker.proteins.fasta"
bterr_rna = "Bterr_rna.gbk"
bterr_pep = "Bterr_protein.gbk"
bimpa_rna = "Bimpa_rna.gbk"
bimpa_pep = "Bimpa_protein.gbk"

missing_pep = "Missing_protein.fa"
missing_rna = "Missing_rna.fa"

StopCodons = ['TAA','TAG','TGA']

################################################################################################################################  

def pad(s,lpad,rpad,padchar = '#'):
  """Pads string s on left with lpad chars, on right with rpad chars all of padchar
  
  """
  return padchar*lpad+s+rpad*padchar
##
  
def chopStr(s,l):
  """Splits string into segments of length l
  
  """
    
  if (l >= len(s)):
    return [s]

  ss = []
  
  i = 0
  for i in range(0,len(s)/l):
    ss.append(s[i*l:min(len(s),i*l+l)])
    
  if (i*l+l < len(s)):
    ss.append(s[i*l+l:])
    
  return ss
##


def printAlign(seq1,seqs,offset,seglen = -1,printcomp = True):
  """Print the alignment by offsetting all sequences in seqs by offset characters
  
  Comparison done to seqs[0]
  
  Use seglen to specify how big the segments should be, negative to not split into segments
  Set printcomp to False to not print the comparison string in the middle
  """
  
  if (type(seqs) != type(list())):
    seqs = [seqs]
  
  seq1 = pad(seq1,max(-offset,0),0,'_')
  for i in range(0,len(seqs)):
      seqs[i] = pad(seqs[i],offset,0,'_')
      lendiff = len(seqs[i])-len(seq1)
      if (lendiff < 0):
          seqs[i] = pad(seqs[i],0,-lendiff,'_')
      else:
          seq1 = pad(seq1,0,lendiff,'_')

  
  comp = ['.']*len(seq1)

  for i in range(0,len(seq1)):
    for j in range(0,len(seqs)):
        if (seq1[i] != ' ' and seq1[i] == seqs[j][i]):
            if (comp[i] != ':'):
                comp[i] = '|'
        else:
            if (comp[i] == '|' and seqs[j][i] != '_'):
                comp[i] = ':'
  comp = ''.join(comp)
      
  if (seglen > 0):
    ssseq1 = chopStr(seq1,seglen)
    sscomp = chopStr(comp,seglen)
    ssseqs = []

    for seq in seqs:
        ssseqs.append(chopStr(seq,seglen))
    
    align = ''
    for i in range(0,len(ssseq1)):
        align += ssseq1[i]+'\n'
        if (printcomp):
          align += sscomp[i]+'\n'
        for j in range(0,len(seqs)):
          align += ssseqs[j][i]+'\n'
        align += '\n'
  else:
    align = seq1+'\n'
    if (printcomp):
        align += comp+'\n'
    for j in range(0,len(seqs)):
        align += seqs[j]+'\n'
    align += '\n'
    
  return align  
##



def ParseCodons(sequence):
    """Iterable that returns the codons in a sequence
    
    No error checking!
    """
    for i in range(0,len(sequence),3):
        yield sequence[i:i+3]
##

def GetORFs(rnaseq):

      orfs = []
      orfs.append(rnaseq)
      orfs.append(rnaseq[1:])
      orfs.append(rnaseq[2:])
      orfs.append(rnaseq.reverse_complement())
      orfs.append(rnaseq.reverse_complement()[1:])
      orfs.append(rnaseq.reverse_complement()[2:])

      #orfs = map(str, orfs)
      return orfs
#

def ScoreORFs(orfs, pepseq):

      scores  = []
      offsets = []
      for rnaseq in orfs:
            seq = str(rnaseq.upper().translate())
            maxscore  = 0
            maxoffset = 0
            for offset in range(-len(pepseq), len(pepseq)):
                  score = 0
                  if (offset >= 0):
                      for i in range(0, min(len(pepseq)-offset, len(seq))):
                          if (pepseq[i+offset] == seq[i]):
                              score += 1
                  else:
                      for i in range(0, min(len(pepseq), len(seq)+offset)):
                          if (pepseq[i] == seq[i-offset]):
                              score += 1

                  #print "%d\t%d" % (offset, score)

                  if (score > maxscore):
                        maxscore  = score
                        maxoffset = offset
                  #
            #
            scores.append(maxscore)
            offsets.append(maxoffset)
            #print printAlign(pepseq, seq, maxoffset, 100)
            #print
      #

      maxidx = sorted(range(6), key=lambda k: scores[k])[-1]

      return (maxidx, scores[maxidx], offsets[maxidx])
#


def CheckTranslation(sequences):
    """Attempt to fix the translation by checking ORFs of all transcripts

    """
    
    aaseq  = list(str(sequences[0]['pepsequence']))
    codons = ['']*len(aaseq)
    printseqs = []

    for rnaseq in sequences:

        orfs = GetORFs(rnaseq['rnasequence'])
        (orf, score, offset) = ScoreORFs(orfs, ''.join(aaseq))

        #sys.stdout.write("\t%d\t%d\t%d\n" % (orf, score, offset))

        translation = list(str(orfs[orf].upper().translate()))
        transcodons = list(ParseCodons(str(orfs[orf]).upper()))

        printseqs.append(pad(''.join(translation),offset,0,'_'))


        for i in range(0, min(len(aaseq)-offset, len(translation))):
            if (aaseq[offset+i] == translation[i]):
                if (codons[offset+i] == '' or codons[offset+i] == transcodons[i]):
                    codons[offset+i] = transcodons[i]
                else:
                    sys.stdout.write('\t\tResidue %d: Multiple codon matches' % (offset+i))
                    sys.stdout.write('\t%s\t%s\t%s\n' % (aaseq[offset+i], codons[offset+i], transcodons[i]))
                #
            # This is a possible gap caused by splicing together exons
            elif (aaseq[offset+i] == 'X' and translation[i] == '*'):
                codons[offset+i] = '*'
        #
    #


    print 
    if (offset >= 0):
        print printAlign(''.join(aaseq), printseqs, 0, 100)
    else:
        print printAlign(''.join(aaseq), printseqs, offset, 100)

    problems = 0
    for i in range(0,len(codons)):
        if (codons[i] == ''):
            sys.stdout.write('\t\tProblem with codon %d\t%s\t%s\n' % (i, aaseq[i], codons[i]))
            problems += 1

    if (problems == 0):
        sequences[0]['rnasequence'] = ''.join(codons).replace('*','')
        sequences[0]['pepsequence'] = ''.join(aaseq).replace('X','')
        return sequences[0]
    else:
        return None
    #
#






def SearchFasta(geneid, geneseq, pepfname, rnafname):
    """Extract genes from matching fasta cds and fasta peptide files 

        (Nasonia and Apis, Megachile)
    """

    # Search in RNA file
    rnafile = open(rnafname,'r')
    found = []
    for record in SeqIO.parse(rnafile,'fasta'):
        if (geneid in record.description):
            seqid = record.id

            # Remove stop codon if present
            if (str(record.seq[-3:]).upper() in StopCodons):
                record.seq = record.seq[:-3]

            found.append({'rnadescription' : record.description, 
                          'rnasequence'    : record.seq,
                          'rnafile'        : rnafname[rnafname.rfind('/')+1:],
                          'pepdescription' : '', 
                          'pepsequence'    : '',
                          'pepfile'        : pepfname[pepfname.rfind('/')+1:]})
    #
    rnafile.close()

    # Not found!
    if (found == []):
        sys.stdout.write("RNA SEQUENCE NOT FOUND!\n")
        return None
    else:
        sys.stdout.write("Found %d RNA sequence(s)..." % len(found))
    #

    # Search found sequences in protein file
    pepfile = open(pepfname, 'r')
    for record in SeqIO.parse(pepfile,'fasta'):          
        if (geneid in record.description):
            sys.stdout.write("ok...")

            # Remove stop codon if present
            if (str(record.seq[-1]) == '*'):
                record.seq = record.seq[:-1]

            # Check if it matches the alignment sequence, and check the translation
            alignseq = str(geneseq).replace('-','')
            if (alignseq == str(record.seq)):
                sys.stdout.write('protein matched...')
  
                for rnaseq in found:
                    rnaseq['pepdescription'] = record.description
                    rnaseq['pepsequence']    = record.seq
                    if (str(record.seq) == str(rnaseq['rnasequence'].upper().translate())):
                        sys.stdout.write("translation ok!")
                        pepfile.close()
                        return rnaseq
                sys.stdout.write("Problem with translation...attempting to fix\n")
                return CheckTranslation(found)
                #
            #
        #
    #
    pepfile.close()

    return None
##


def SearchGenbank(geneid, geneseq, gbfile):
    """Extract cDNA from Genbank file

        (Bombus)
    """

    rnafile = open(gbfile, 'r')

    for record in SeqIO.parse(rnafile,'genbank'):
        for feature in record.features:             
            for featurequal in feature.qualifiers:
                for featurecontent in feature.qualifiers[featurequal]:
                    if (geneid in featurecontent):
                        sys.stdout.write("ok...")

                        if (feature.type == 'CDS'):
                            sys.stdout.write("Found RNA sequence(s)...")
            
                            pepsequence = feature.qualifiers['translation'][0]
                            rnasequence = feature.extract(record).seq
                            
                            # Remove stop codon if present
                            if (str(rnasequence[-3:]) in StopCodons):
                                rnasequence = rnasequence[:-3]
                            if (str(pepsequence[-1]) == '*'):
                                pepsequence = pepsequence[:-1]

                            alignseq = str(geneseq).replace('-','')
                            if (alignseq == str(pepsequence)):
                                sys.stdout.write('protein matched...')
          
                                record = {'rnasequence'    : rnasequence,
                                            'gbfilename'     : gbfile[gbfile.rfind('/')+1:],
                                            'rnadescription' : record.id+' '+record.description,
                                            'pepdescription' : feature.qualifiers['protein_id'][0]+' '+feature.qualifiers['product'][0],
                                            'pepsequence'    : pepsequence}                              

                                if (str(pepsequence) == str(rnasequence.translate())):
                                    sys.stdout.write("translation ok!")
                                    rnafile.close()
                                    return record
                                else:
                                    return CheckTranslation([record])
            #
    #
    rnafile.close()
                    
    return None
##


# def SearchAmell(geneid, geneseq):
#     """Deprecated - Match AA sequence (used when gene ids are mismatched)

#     """

#     alignseq = str(geneseq).replace('-','')

#     rnafile = open(genepath+'/'+amell_rna)
#     for record in SeqIO.parse(rnafile,'genbank'):
#         for feature in record.features:
#             if (feature.type == 'CDS'):
#                 pepsequence = feature.qualifiers['translation'][0]

#                 # Remove stop codon if present
#                 if (str(pepsequence[-1]) == '*'):
#                     pepsequence = pepsequence[:-1]

#                 if (alignseq == str(pepsequence)):
#                     sys.stdout.write('protein matched...')
#                     rnasequence = feature.extract(record).seq

#                     # Remove stop codon if present
#                     if (str(rnasequence[-3:]) in StopCodons):
#                         rnasequence = rnasequence[:-3]
 

#                     if (str(pepsequence) == str(rnasequence.translate())):
#                         sys.stdout.write("translation ok!")
#                         rnafile.close()
#                         return rnasequence
#         #
#     #

#     rnafile.close()
#     return None
# ##




################################################################################################################################  
start = time.time()

if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

sys.stdout.write("Processing alignments...\n")
errfile = open(outputpath+'/Errors.csv','w')
pepfile = open(outputpath+'/AAsequences.fa','w')
dnafile = open(outputpath+'/DNAsequences.fa','w')
idsfile = open(outputpath+'/Stats.txt','w')
statfile= open(outputpath+'/Alignments.csv','w')
statfile.write('Alignment\tLength\tAMELL\tAFLOR\tNVITR\tMROTU\tBTERR\tBIMPA\n')
i = 1
one_each    = 0
missing_seq = 0
for filename in os.listdir(inputpath):
    if (fnmatch(filename,pattern)):
        sys.stdout.write("%3d. " % i+filename+":\n")
        alnfile = open(inputpath+'/'+filename,'r')

        amell = 0
        aflor = 0
        nvitr = 0
        mrotu = 0
        bterr = 0
        bimpa = 0
        for record in SeqIO.parse(alnfile,'fasta'):
            organism = record.description[:5].upper()
            geneid   = record.description.split()[-1]
            geneseq  = record.seq
            sys.stdout.write('\t'+organism+' - '+geneid+':\t')

            if (organism == 'AMELL'):
                  result = SearchFasta(geneid, geneseq, genepath+'/AMELL/'+amell_pep, genepath+'/AMELL/'+amell_rna)
                  amell += 1
            elif (organism == 'AFLOR'):
                  result = SearchFasta(geneid, geneseq, genepath+'/AFLOR/'+aflor_pep, genepath+'/AFLOR/'+aflor_rna)                  
                  aflor += 1
            elif (organism == 'NVITR'):
                  result = SearchFasta(geneid, geneseq, genepath+'/NVITR/'+nvitr_pep, genepath+'/NVITR/'+nvitr_rna)                  
                  nvitr += 1
            elif (organism == 'MROTU'):
                  result = SearchFasta(geneid, geneseq, genepath+'/MROTU/'+mrotu_pep, genepath+'/MROTU/'+mrotu_rna)                  
                  mrotu += 1
            elif (organism == 'BTERR'):
                  result = SearchGenbank(geneid, geneseq, genepath+'/BTERR/'+bterr_rna)
                  bterr += 1
            elif (organism == 'BIMPA'):
                  result = SearchGenbank(geneid, geneseq, genepath+'/BIMPA/'+bimpa_rna)
                  bimpa += 1
            else:
                sys.stdout.write('UNKNOWN ORGANISM!')
            sys.stdout.write('\n')

            if (result == None):  # Check missing genes file
                sys.stderr.write("\tLooking for gene among missing alignments")
                result = SearchFasta(geneid, geneseq, genepath+'/'+missing_pep, genepath+'/'+missing_rna)

            # Write output
            if (result == None):  # Couldn't find sequence automatically
                sys.stderr.write("\tError with %s\n" % organism)
                errfile.write(filename+"\t"+organism+"\t"+geneid+"\n")  
            else:
                pepfile.write('>%s\n%s\n' % (record.description, result['pepsequence']))
                dnafile.write('>%s\n%s\n' % (record.description, str(result['rnasequence']).upper().replace('U','T')))

                idsfile.write(record.description+'\n')
                idsfile.write('\tOrthoDB:\t'+filename+'\n')
                if ('rnafile' in result):
                    idsfile.write('\tCDS file:\t'+result['rnafile']+'\t\n')
                    idsfile.write('\tProtein file:\t'+result['pepfile']+'\n')
                elif ('gbfilename' in result):
                    idsfile.write('\tGenbank file:\t'+result['gbfilename']+'\n')
                else: 
                    pass
                idsfile.write('\tCDS:\t\t'+result['rnadescription']+'\n')
                idsfile.write('\tProtein:\t'+result['pepdescription']+'\n')
            #

        #

        statfile.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (filename, len(geneseq), amell, aflor, nvitr, mrotu, bterr, bimpa))

        if (amell == 1 and nvitr == 1 and bterr == 1 and bimpa == 1):
            one_each += 1

        if (amell == 0 or nvitr == 0 or bterr == 0 or bimpa == 0):
            missing_seq += 1

        alnfile.close()
        i += 1
    #
#
errfile.close()
pepfile.close()
dnafile.close()
idsfile.close()

sys.stdout.write('\n')
sys.stdout.write('%d alignments in total\n' % (i-1))
sys.stdout.write('%d alignments with exactly one sequence from each species\n' % one_each)
sys.stdout.write('%d alignments missing some species\n' % missing_seq)


end = time.time()
print "Total time taken: "+str(end-start)+" seconds"  