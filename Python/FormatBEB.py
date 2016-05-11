
import sys, os, time, string, stat, shutil

from Bio import AlignIO
from fnmatch import fnmatch
from optparse import OptionParser



# Format BEB into BEB values for each of the sequences in the alignment
#
#
################################################################################################################################
# Parameters
################################################################################################################################

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-b","--bebpath",
                  dest = "bebpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing result from ExtractBEB.py [default = %default]")

parser.add_option("-i","--msapath",
                  dest = "msapath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing the MSAs and Gblocks results [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "filename to store output in [required]")

parser.add_option("-m","--msaprogram",
                  dest = "msaprogram",
                  default = "mafft",
                  metavar = "mafft/muscle/probcons/prank/prographmsa",
                  help = "MSA program to use [default = %default]")


(options,args) = parser.parse_args()

msapath = os.path.abspath(options.msapath)
bebpath = os.path.abspath(options.bebpath)
outputpath = os.path.abspath(options.outputpath)
msaprogram = options.msaprogram.lower()

################################################################################################################################  

def GetGBlocks(group):

      gblocks = open(msapath+'/'+group+'.'+msaprogram+'.dna.fas-gb.txt','r')

      for line in gblocks:
      	if (line[:7] == 'Flanks:'):
      		flanks = map(int,line[7:].replace('[','').replace(']','').strip().split())
      gblocks.close()

      return flanks
#

def ParseAlign(fname, format):
      align = dict()
      msa   = AlignIO.read(fname,format)

      for record in msa:
            align[record.id[:5]] = str(record.seq)

      return align
#
 
def GetTrimSeq(sequence, flanks):

      trimseq = '';
      for i in range(0, len(flanks), 2):
            start = (flanks[i]-1)/3
            end   = (flanks[i+1]/3)
            trimseq += sequence[start:end]

      return trimseq
#


def MapIdx(trimidx, flanks):
      """From the trimmed indices to the original msa indices

      """
      seqlen = 0
      for i in range(0, len(flanks), 2):
            start = (flanks[i]-1)/3
            end   = (flanks[i+1]/3)
            
            seqlen += (end - start)
            if (seqlen > trimidx):
                  seqlen -= (end - start)
                  return start + (trimidx - seqlen)

      return -1
#

def GetTrimIds(flanks):
      
      trimids = []
      for i in range(0, len(flanks), 2):
            start = (flanks[i]-1)/3
            end   = (flanks[i+1]/3)
            
            trimids += range(start,end)

      return trimids
#

def GetBEB(group, model):

      lines   = []
      bebfile = open(bebpath+'/'+group+'.'+model+'.BEB')

      for line in bebfile:
            lines.append(line.split())
      bebfile.close()

      return lines
#

def ExtendBEB(trim_beb, trim_ids, orig_msa, outputfname):

      for record in orig_msa:
            outfile = open(outputfname+'.'+record,'w')

            seq = orig_msa[record]

            curridx = 1
            trimidx = 1
            outfile.write('\t'.join(trim_beb[0])+'\n')
            for i in range(0, len(seq)):

                  if (seq[i] != '-'):
                        outfile.write("%i\t%s" % (curridx, seq[i]))
                        curridx += 1

                        if (i in trim_ids):                  
                             outfile.write('\t'+'\t'.join(trim_beb[trimidx][2:])+'\n')
                             trimidx += 1
                        else:
                             outfile.write('%s\n'%('\t-1'*(len(trim_beb[0])-2)))
                  #
            outfile.close()
      #
#


def ExtendBEBAlign(trim_beb, trim_ids, orig_msa, outputfname):


      # Output alignment
      alnout = open(outputfname, 'w')
      alnout.write('Site\tGblocks')
      index  = dict()
      for record in orig_msa:
            alnout.write('\t'+record+'_Site\t'+record+'_AA')
            index[record] = 0

      alnout.write('\t'+'\t'.join(trim_beb[0][2:])+'\n')
      trimidx = 1
      for i in range(0, len(orig_msa[record])):

            if (i in trim_ids):                          
                  alnout.write('%d\tG' % (i+1))
            else:
                  alnout.write('%d\tX' % (i+1))

            for org in orig_msa:
                  if orig_msa[org][i] != '-':
                        index[org] += 1
                        ind = index[org]
                  else:
                        ind = -1

                  alnout.write('\t%d\t%c' % (ind, orig_msa[org][i]))

            if (i in trim_ids):                  
                  alnout.write('\t'+'\t'.join(trim_beb[trimidx][2:])+'\n')
                  trimidx += 1
            else:
                  alnout.write('%s\n'%('\t-1'*(len(trim_beb[0])-2)))

      alnout.close()
#




################################################################################################################################  

start = time.time()
#models = ['M8', 'BS.Tree1']

if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

for filename in sorted(os.listdir(bebpath), reverse=True):
      if (fnmatch(filename,'*.BEB')):
            group = filename[:filename.find('.')]
            model = filename[filename.find('.')+1:filename.rfind('.')]

            flanks   = GetGBlocks(group)
            trim_ids = GetTrimIds(flanks)

            orig_msa = ParseAlign(msapath+'/'+group+'.'+msaprogram+'.fas', 'fasta')
            trim_msa = ParseAlign(msapath+'/'+group+'.'+msaprogram+'.phylip', 'phylip')    

            sys.stdout.write(group + '\t' + model +'\n')
            trim_beb = GetBEB(group, model)

            ExtendBEB(trim_beb, trim_ids, orig_msa, outputpath+'/'+group+'.'+model+'.BEB')
            ExtendBEBAlign(trim_beb, trim_ids, orig_msa, outputpath+'/'+group+'.'+msaprogram+'.'+model+'.BEB')

            # Output sequences
            seqout  = open(outputpath+'/'+group+'.fa','w')
            for record in orig_msa:
                  seqout.write('>'+record+'\n'+orig_msa[record].replace('-','')+'\n')
            seqout.close()

            



            
      #


end = time.time()
print "Total time taken: "+str(end-start)+" seconds"
