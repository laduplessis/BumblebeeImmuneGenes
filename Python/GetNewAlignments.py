import sys, os, time, math, shutil, string
from fnmatch import fnmatch
from optparse import OptionParser
from GeneticCode import parseCodons, translateSeq
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from subprocess import Popen, PIPE


# Produce alignments that can be used with PAML
#     - Extract sequences for each group from the files created with ExtractDNA.py
#     - Skip groups that have missing sequences for some species
#     - Check translation again with OrthoDB alignments
#     - If there is more than one sequence of a species:
#           * Perform a new alignment on the translated sequences (AA) with the chosen MSA program for 
#             every possible combination of sequences so that the MSA has only one sequence per species
#           * Use Gblocks (with strict settings) to restrict the alignments to the best matching areas only
#           * Select the alignment with the longest well-aligned parts (as determined by Gblocks)
#           * If there are more alignments with disjoint sets of sequences select alignments by longest well aligned parts until 
#             the maximum number is selected
#       Else:
#           * Perform a new alignment on the translated sequences with the chosen MSA program
#     - Trim the alignment:
#           * Using Gblocks
#     - Skip groups where the trimmed alignment is shorter than the minimum length
#     - Save trimmed cDNA alignments in Phylip format so it can be used with PAML
#
# Need to have MSA programs and Gblocks installled.  Check that programs are in your path.  
#     (For prographmsa the path is hardcoded for the moment)
#     (Some problems with Muscle and Prank)
#
# Output
#     - Alignments.csv:  Exactly like Alignments.csv from ExtractDNA.py, but only for groups selected here
#     - Skipped.csv:     Groups which have been skipped for whatever reason
#     - Colsremoved.csv: List of columns removed from each of the alignments
#     For each group that wasn't skipped:
#           - <group>.aa.fa:                          Protein sequences without gaps
#           - <group>.dna.fa:                         Corresponding cDNA sequences
#           - <group>.<msaprogram>.fas:               New alignment (AA, ids remapped)
#           - <group>.<msaprogram>.dna.fas            New alignment (cDNA, ids remapped)
#           - <group>.<msaprogram>.dna.fas-gb         Trimmed alignment from Gblocks
#           - <group>.<msaprogram>.dna.fas-gb.txt     Gblocks statistics
#           - <group>.<msaprogram>.dna.phylip         Trimmed cDNA alignment in Phylip, ready for Paml
#           - <group>.<msaprogram>.phylip             Trimmed AA alignment in phylip
#           - <group>.map                             Mapping of ids to the Phylip files
#           - <group>.ids                             Identifiers of the original sequences used
#           - <group>.combinations.tar                Temporary files from finding the best combination of sequences to use
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
                  help = "Directory containing sequences in Fasta format (from ExtractDNA.py) [default = %default]")

parser.add_option("-a","--alignpath",
                  dest = "alignpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing OrthoDB alignments in Fasta format [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "filename to store output in [required]")

parser.add_option("-l","--minlength",
                  dest = "minlength",
                  default = "100",
                  metavar = "integer",
                  help = "Minimum length for an alignment (in codons, after trimming) [default = %default]")

parser.add_option("-s","--speciesmin",
                  dest = "speciesmin",
                  default = "1,1,1,1,1,1",
                  metavar = "comma separated integer list",
                  help = "Minimum number of species in each alignment to be considered (Amell, Aflor, Nvitr, Mrotu, Bimpa, Bterr) [default = %default]")

parser.add_option("-S","--speciesmax",
                  dest = "speciesmax",
                  default = "100,100,100,100,100,100",
                  metavar = "comma separated integer list",
                  help = "Maximum number of species in each alignment to be considered (Amell, Aflor, Nvitr, Mrotu, Bimpa, Bterr) [default = %default]")

parser.add_option("-m","--msaprogram",
                  dest = "msaprogram",
                  default = "mafft",
                  metavar = "mafft/muscle/probcons/prank/prographmsa",
                  help = "MSA program to use [default = %default]")

parser.add_option("-g","--gblockspars",
                  dest = "gblockspars",
                  default = "",
                  metavar = "parameters",
                  help = "Parameters for Gblocks (strict/relaxed/none) - use Gblocks if this is specified [default = %default]")

parser.add_option("-c","--combalignpath",
                  dest = "combalignpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing result from GetNewAlignment to use for combinations [default = %default]")

(options,args) = parser.parse_args()

inputpath    = os.path.abspath(options.inputpath)+'/'
alignpath    = os.path.abspath(options.alignpath)+'/'
outputpath   = os.path.abspath(options.outputpath)+'/'
minlength    = int(options.minlength)
speciesmin   = map(int, options.speciesmin.split(','))
speciesmax   = map(int, options.speciesmax.split(','))
msaprogram   = options.msaprogram.lower()
gblockspars  = options.gblockspars
if (options.combalignpath != ""):
      combalignpath = os.path.abspath(options.combalignpath)+'/'
else:
      combalignpath = options.combalignpath

protdist = '/Users/louis/Documents/Packages/phylip-3.69/exe/protdist'
if (msaprogram == 'prank'):
      msaprogram = 'prank.best'

speciesnames = ['Amell','Aflor','Nvitr','Mrotu','Bterr','Bimpa']
StopCodons   = ['TAA','TAG','TGA']


################################################################################################################################  

def GetDNAAlignment(alignseq, dnaseq):

      dnaalignseq = ''
      codons = list(parseCodons(dnaseq))
      i = 0
      for aa in alignseq:
            if (aa == '-'):
                  dnaalignseq += '---'
            else:
                  dnaalignseq += codons[i]
                  i += 1
            #
      #

      return dnaalignseq

#      


def GetSequences(fname, seqid):
      for record in SeqIO.parse(fname, 'fasta'):
            if (record.description == seqid):
                  return record.seq
#


def ExtractOrthoDBAlignment(name, skipfile):

      orthodbin = open(alignpath+name)
      dnafile   = open(inputpath+'DNAsequences.fa')
      aafile    = open(inputpath+'AAsequences.fa')

      idmap     = dict()
      alignment = dict()
      speciesnr = [0]*len(speciesnames)
      for record in SeqIO.parse(orthodbin,'fasta'):
            sys.stdout.write('\t'+record.description+'\t')
            orthodbseq = record.seq
            dnaseq     = GetSequences(inputpath+'DNAsequences.fa', record.description)
            aaseq      = GetSequences(inputpath+'AAsequences.fa', record.description)
            

            
            # Check proteins are equal
            if (str(orthodbseq).replace('-','') == str(aaseq)):
                  sys.stdout.write('matched proteins...')
            elif (str(orthodbseq).replace('-','').replace('X','') == str(aaseq)):
                  sys.stdout.write('matched proteins (saved)...')
                  orthodbseq = Seq.Seq(str(orthodbseq).replace('X','-'))
            else:
                  sys.stdout.write('protein MISMATCH!\n')
                  return None

            # Check translation
            if (str(dnaseq.translate()) == str(aaseq)):
                  if (len(dnaseq) > 3*len(aaseq)):
                        dnaseq = dnaseq[:3*len(aaseq)]
                  sys.stdout.write('translation ok...\n')
            else:
                  sys.stdout.write('translation FAIL!\n')
                  return None

            org = record.description[:5]
            for i in range(0,6):
                  if (org == speciesnames[i].upper()):
                        speciesnr[i] += 1      
                        idmap[org+str(speciesnr[i])] = record.description
                        break
                  #
            #       


            sequence = {'id'         : org+str(speciesnr[i]),
                        'description': record.description,
                        'orthodbseq' : str(orthodbseq.upper()), 
                        'orthodbdna' : str(GetDNAAlignment(orthodbseq, dnaseq).upper()),
                        'dnaseq'     : str(dnaseq.upper()),
                        'aaseq'      : str(aaseq.upper())}

            if (org in alignment):
                  alignment[org].append(sequence)
            else:
                  alignment[org] = [sequence]

      #
      orthodbin.close()
      dnafile.close()
      aafile.close()

      return (alignment, idmap)
#

################################################################################################################################  

def GetAllowedCombinations(group, alignpath):
      """Retrieve all of the combinations allowed for a group from another directory of alignments

      """

      allowed = []
      print group
      for filename in os.listdir(alignpath):
          if (fnmatch(filename,group+'*.ids')):
              seqs = []
              idfile = open(alignpath+filename,'r')
              for line in idfile:
                  seqs.append(line[line.find('\t'):].strip())
              idfile.close()
              allowed.append(seqs)
      #
      print "Allowed "+str(len(allowed))+" combinations"
      return allowed
#

def CheckCompatible(combination, alignment):
      """Check if alignment is compatible with the combination

      (combination is a subset of the alignment)
      """

      for seq in combination:
            found = False
            for org in alignment:
                  if (org['description'] == seq):
                        found = True
                        break

            if (not found):
                return False
      #
      return True
#



def GetCombinations(speciesid, tempalignment, alignment, alncombs, allowed=None):
      """Get all combinations of sequences from an alignment so that each combination has only one sequence per species

      (Recursive, could use too much memory if there are too many sequences!)
      """

      if (speciesid >= len(speciesnames)):
            for org in tempalignment:
                  sys.stdout.write('\t'+org['id'])

            if (allowed != None):
                  compatible = False
                  for comb in allowed:
                        if (CheckCompatible(comb,tempalignment)):
                              compatible = True
                              break
                  #
                  if (compatible):
                        sys.stdout.write('\tOK!\n')
                        alncombs.append(tempalignment)       
                  else:
                        sys.stdout.write("\tincompatible\n")

            else:
                  sys.stdout.write('\n')
                  alncombs.append(tempalignment)

      else:
            org = speciesnames[speciesid].upper()

            if (org in alignment):
                  for sequence in alignment[org]:
                        GetCombinations(speciesid + 1, tempalignment + [sequence], alignment, alncombs, allowed)
            else:
                  GetCombinations(speciesid + 1, tempalignment, alignment, alncombs, allowed)
      #
#

def runMSASimple(group):
      """Run specified MSA program on the group
      
      Run MSA on AA sequences, no input guide tree
      
      Group is found from global inputpath
      Program to use from global msaprogram
      Output stored in global outputpath (after creating appropriate subdirectory
      
      """
      
      # Check if input file exists
      if not (os.path.exists(outputpath+group+'.aa.fa')):
            sys.stderr.write("\tWarning! Family "+str(group)+" not found!\n")                 
            return

      print "\n\tStarting MSA computation for group "+group
      start = time.time()

      # Run MSA program
      if (msaprogram == 'mafft'):
            handle = Popen("mafft-linsi %s > %s" % (outputpath+group+'.aa.fa',outputpath+group+'.mafft.fas'), stdout=PIPE, stderr=PIPE, shell=True)     
      elif (msaprogram == 'muscle'):            
            handle = Popen("muscle -in %s -out %s" % (outputpath+group+'.aa.fa',outputpath+group+'.muscle.fas'), 
                        stdout=PIPE, stderr=PIPE, shell=True)
      elif (msaprogram == 'probcons'):
            handle = Popen("probcons %s > %s" % (outputpath+group+'.aa.fa',outputpath+group+'.probcons.fas'),
                        stdout=PIPE, stderr=PIPE, shell=True)
      elif (msaprogram == 'prank.best'):
            handle = Popen("prank -d=%s -o=%s" % (outputpath+group+'.aa.fa',outputpath+group+'.prank'),
                        stdout=PIPE, stderr=PIPE, shell=True)
            #try:
            #      print os.path.exists(outputpath+group+".prank.best.fas")
            #      shutil.move(outputpath+group+".prank.best.fas",outputpath+group+".prank.fas")
            #except:
            #      sys.stderr.write("Warning! Could not rename PRANK output file!\n")      
      elif (msaprogram == 'prographmsa'):
            handle = Popen("/Users/louis/Documents/Packages/ProGraphMSA/PrographMSA+TR.sh --fasta %s -o %s" % (outputpath+group+'.aa.fa',outputpath+group+'.prographmsa.fas'),
                        stdout=PIPE, stderr=PIPE, shell=True)
      else:
            sys.stderr.write("\tWarning! Unknown MSA program "+msaprogram+"!\n")
            
      # Read and process errors from MSA program
      out = handle.stdout.read()

      err = handle.stderr.read()
      if (err != ""):
            sys.stderr.write("\tWarning! Errors encountered!\n")
            #sys.stderr.write(err)
      

      end = time.time()
      print "\tMSA computed ("+str(end-start)+" seconds)"
      print
##

def GetNewAlignment(group, alignment):
      """Write ungapped AA and DNA sequences to a file and run MSA program

      Returns Alignment object
      """
      dnafile      = open(outputpath+'/'+group+'.dna.fa','w')
      aafile       = open(outputpath+'/'+group+'.aa.fa','w')
      for record in alignment:
            dnafile.write   ('>%s\n%s\n' % (record['id'], record['dnaseq']))
            aafile.write    ('>%s\n%s\n' % (record['id'], record['aaseq']))
      dnafile.close()
      aafile.close()

      runMSASimple(group)
##


################################################################################################################################  


def ArchiveFiles(outputpath, group):
      """Archive and remove temporary alignments (all combinations)

      """

      sys.stdout.write("\tArchiving files...\n")

      here = os.getcwd()
      os.chdir(outputpath)

      pattern = group+'_*'

      callstring = "tar cvf "+group+".combinations.tar "+pattern
      handle = Popen(callstring, stdout=None, stderr=PIPE, shell=True)
        
      # Read and process errors from compressing
      err = handle.stderr.read()
      if (err != ""):
            sys.stderr.write(err)

      i = 0
      for filename in os.listdir(outputpath):
          if (fnmatch(filename,pattern)):
                  try:
                        os.remove(filename)
                        i += 1
                  except Exception:
                        sys.stdout.write("\tERROR deleting %s\n" % filename)
      sys.stdout.write("\tDeleted %d files\n" % i)

      os.chdir(here)
#

def CompressFiles(outputpath, pattern):
      sys.stdout.write("\tCompressing files...\n")

      here = os.getcwd()
      os.chdir(outputpath)

      
      callstring = "gzip "+pattern
      handle = Popen(callstring, stdout=None, stderr=None, shell=True)
        
      # Read and process errors from compressing
      #err = handle.stderr.read()
      #if (err != ""):
      #      sys.stderr.write(err)
#


################################################################################################################################


def runGblocks(fname, pars):
      """Trim alignment using Gblocks

      Returns Alignment object
      """

      print "\n\tStarting Gblocks computation for "+fname[fname.rfind('/')+1:]
      start = time.time()
      
      handle = Popen("Gblocks %s %s " % (fname,pars), stdout=PIPE, stderr=PIPE, shell=True)

      # Read and process errors from MSA program
      out = handle.stdout.read()

      err = handle.stderr.read()
      if (err != ""):
            sys.stderr.write("\tWarning! Errors encountered!\n")
            #sys.stderr.write(err)


      end = time.time()
      print "\tGblocks finished ("+str(end-start)+" seconds)\n"

      statfile = open(fname+'-gb.txt','r')
      line = statfile.readline()
      while (line.find("New number of positions") < 0):
            line = statfile.readline()
      trimlength = int(line[line.find(':')+1:line.find('(')])
      return trimlength
##



################################################################################################################################  
start = time.time()

if (not os.path.exists(outputpath)):
      os.mkdir(outputpath)
#

statfile = open(outputpath+'Alignments.csv','w')
statfile.write('Alignment\tLength\tAMELL\tAFLOR\tNVITR\tMROTU\tBTERR\tBIMPA\n')
skipfile = open(outputpath+'Skipped.csv','w')
skipfile.write('Alignment\tLength\tAMELL\tAFLOR\tNVITR\tMROTU\tBTERR\tBIMPA\n')
colsfile = open(outputpath+'Colsremoved.csv','w')
colsfile.write('Alignment\tNr\tRemoved\n')
alnfile  = open(inputpath+'Alignments.csv','r')
alnfile.readline()

for line in alnfile:
      parts   = line.split()
      name    = parts[0]
      length  = int(parts[1])
      species = map(int,parts[2:])

      sys.stdout.write('\nExtracting '+name+': ')

      ok = True
      for i in range(0,len(speciesnames)):
            if (species[i] < speciesmin[i] or species[i] > speciesmax[i]):
                  sys.stdout.write('Wrong number of %s sequences...skipping\n' % speciesnames[i])
                  skipfile.write(line)
                  ok = False
                  break
      if (not ok):
            continue
      else:
            sys.stdout.write('\n')


      # Alignment has all the required species, extract and check the sequences
      (alignment, idmap) = ExtractOrthoDBAlignment(name, skipfile)

      if (alignment == None):
             sys.stdout.write('\tPROBLEM WITH ALIGNMENT...skipping\n')
             skipfile.write(line)
      else:

            ######################
            # Get best alignment #
            ######################
            alncombs = []
            if (combalignpath != ""):
                  allowed  = GetAllowedCombinations(name[:name.find('.')], combalignpath)
                  GetCombinations(0, [], alignment, alncombs, allowed)
            else:
                  GetCombinations(0, [], alignment, alncombs)

            # Get group name and get number of groups to extract
            orthodbgroup = name[:name.find('.')]
            nrseqs = []
            for org in alignment:
                  nrseqs.append(len(alignment[org]))
            nrperms = min(nrseqs)


            # Save id map
            mapfile = open(outputpath+orthodbgroup+'.map','w')
            for seqid in idmap:
                  mapfile.write(seqid+'\t'+idmap[seqid]+"\n")
            mapfile.close()


            # Get scores for all combinations
            sys.stdout.write('Finding best combination(s) of gene duplications...\n')
            scores = []
            for i in range(0,len(alncombs)):
                  GetNewAlignment(orthodbgroup+"_%d" % i, alncombs[i])
                  strict = '-t=p -p=t -b3=8  -b4=10 -b5=n'
                  trimlength = runGblocks(outputpath+orthodbgroup+"_%d." % i +msaprogram+".fas", strict)
                  scores.append(trimlength)

            # Sort by score and get ids of alignments to save 
            removed  = []
            selected = []
            for comb in sorted(range(0,len(scores)), key=lambda x : scores[x], reverse=True):
                  #sys.stdout.write('%d\t%d\t' % (comb, scores[comb]))
                  #for s in alncombs[comb]:
                  #      sys.stdout.write(s['id']+' ')
                  #sys.stdout.write('\n')

                  unique = True
                  for seq in alncombs[comb]:
                        if (seq['id'] in removed):
                          unique = False
                          break

                  if (unique):
                        for seq in alncombs[comb]:
                            removed.append(seq['id'])   
                        selected.append(comb)
                        if (len(selected) == nrperms):
                              break
                        

            # Save selected alignments
            ind = 1
            for comb in selected:

                  ##############
                  # Copy files #
                  ##############

                  if (len(selected) == 1):
                        outgroup = orthodbgroup
                  else:
                        outgroup = orthodbgroup+'-%d' % ind
                        ind += 1

                  shutil.copy(outputpath+orthodbgroup+'_%d.' % (comb) + msaprogram+'.fas', 
                              outputpath+outgroup+'.'+msaprogram+'.fas')
                  shutil.copy(outputpath+orthodbgroup+'_%d.aa.fa' % (comb), 
                              outputpath+outgroup+'.'+msaprogram+'.aa.fa')
                  shutil.copy(outputpath+orthodbgroup+'_%d.dna.fa' % (comb), 
                              outputpath+outgroup+'.'+msaprogram+'.dna.fa')

                  # Read untrimmed AA msa
                  msa = AlignIO.read(outputpath+outgroup+'.'+msaprogram+'.fas', 'fasta')


                  ######################################
                  # Propagate to DNA and save used ids #
                  ######################################
                  dnaseqs = list(SeqIO.parse(outputpath+outgroup+'.'+msaprogram+'.dna.fa','fasta'))
                  idfile  = open(outputpath+outgroup+'.'+msaprogram+'.ids','w')
                  for dnaseq in dnaseqs:
                        found = False
                        for aaseq in msa:
                              if (dnaseq.description == aaseq.description):
                                    dnaseq.seq = Seq.Seq(GetDNAAlignment(str(aaseq.seq), str(dnaseq.seq)))
                                    idfile.write(dnaseq.id[:5]+'\t'+idmap[dnaseq.id]+'\n')
                                    dnaseq.id  = dnaseq.id[:5]
                                    found = True
                                    break;
                        if (not found):
                              sys.stdout.write("Can't find "+dnaseq.description+" in "+orthodbgroup+"!\n")
                              sys.exit()
                        #
                  idfile.close()
                  dnamsa = Align.MultipleSeqAlignment(dnaseqs)
                  AlignIO.write(dnamsa,outputpath+outgroup+'.'+msaprogram+'.dna.fas',"fasta")


                  ##################
                  # Trim alignment #
                  ##################
                  # (on DNA MSA this time)
                  if (gblockspars == "strict"):
                        strict = '-t=c -p=t -b3=8  -b4=10 -b5=n'
                        runGblocks(outputpath+outgroup+'.'+msaprogram+'.dna.fas', strict)
                        dnamsa = AlignIO.read(outputpath+outgroup+'.'+msaprogram+'.dna.fas-gb', 'fasta')

                  elif (gblockspars == "relaxed"): 
                        n = math.floor(len(dnamsa)/2+1)
                        relaxed = '-t=c -p=t -b2=%d -b3=10 -b4=5  -b5=h' % n
                        runGblocks(outputpath+outgroup+'.'+msaprogram+'.dna.fas', relaxed)
                        dnamsa = AlignIO.read(outputpath+outgroup+'.'+msaprogram+'.dna.fas-gb', 'fasta')

                  else:
                        sys.stdout.write("\tUnknown Gblocks parameters, not trimming\n")


                  #############################
                  # Check length of alignment #
                  ############################# 
                  if (len(dnamsa[0]) / 3 < minlength):
                        sys.stdout.write('Alignment too short...skipping\n')
                        skipfile.write(outgroup+'\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % \
                             (len(dnamsa[0]) / 3, species[0], species[1], species[2], species[3], species[4], species[5]))
                  else:
                        statfile.write(outgroup+'\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % \
                             (len(dnamsa[0]) / 3, species[0], species[1], species[2], species[3], species[4], species[5]))


                        ################################
                        # Save Phylip cDNA output file #
                        ################################
                        # (PAML needs "I" on first line of phylip output file)
                        AlignIO.write(dnamsa,outputpath+outgroup+'.'+msaprogram+'.dna.phylip_temp',"phylip")
                        infile  = open(outputpath+outgroup+'.'+msaprogram+'.dna.phylip_temp','r')
                        outfile = open(outputpath+outgroup+'.'+msaprogram+'.dna.phylip','w')
                        i = 0
                        for line in infile:
                              if (i == 0):
                                    outfile.write(line[:-1]+' I\n')
                              else:
                                    outfile.write(line)
                              i += 1
                        outfile.close()
                        infile.close()
                        os.remove(outputpath+outgroup+'.'+msaprogram+'.dna.phylip_temp')


                        #########################################################
                        # Translate to AA alignment and save Phylip output file #
                        #########################################################
                        aaseqs = []
                        for dnaseq in dnamsa:
                              aaseq = SeqRecord.SeqRecord(Seq.Seq(translateSeq(str(dnaseq.seq))), id = dnaseq.id[:5])
                              aaseqs.append(aaseq)
                        aamsa = Align.MultipleSeqAlignment(aaseqs)
                        AlignIO.write(aamsa, outputpath+outgroup+'.'+msaprogram+'.phylip',"phylip")

            ##

            ###########################
            # Archive temporary files #
            ###########################

            ArchiveFiles(outputpath, orthodbgroup)

     #
skipfile.close()
statfile.close()
colsfile.close()
alnfile.close()

# Compress archived files
CompressFiles(outputpath,'*.combinations.tar')




end = time.time()

print "Total time taken: "+str(end-start)+" seconds"