import sys, os, time, math, shutil
from fnmatch import fnmatch
from optparse import OptionParser
from subprocess import Popen, PIPE

from Bio import AlignIO

# Trim alignments with Gblocks
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
                  help = "Input (DNA) sequences in Fasta format (or path) [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "path to store output files in [required]")

parser.add_option("-g","--gblockspars",
                  dest = "gblockspars",
                  default = "strict",
                  metavar = "parameters",
                  help = "Parameters for Gblocks (strict/relaxed) - use Gblocks if this is specified [default = %default]")

parser.add_option("-p","--pattern",
                  dest = "pattern",
                  default = "*.fas",
                  metavar = "",
                  help = "Pattern to search for (in quotes) [default = %default]")

(options,args) = parser.parse_args()

inputpath    = os.path.abspath(options.inputpath)+'/'
outputpath   = os.path.abspath(options.outputpath)+'/'
pattern      = options.pattern
gblockspars  = options.gblockspars



gblocks = "gblocks"

################################################################################################################################


def runGblocks(fname, pars):
      """Trim alignment using Gblocks

      Returns Alignment object
      """

      print "\n\tStarting Gblocks computation for "+fname[fname.rfind('/')+1:]
      start = time.time()
      
      handle = Popen("%s %s %s " % (gblocks, fname, pars), stdout=PIPE, stderr=PIPE, shell=True)

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

statfile = open(outputpath+"gblocks_"+gblockspars+".csv", "w")
statfile.write("Alignment\tRaw Length\tTrimmed Length\tRetained\n")

for filename in os.listdir(inputpath):
    if (fnmatch(filename,pattern)):

      sys.stdout.write("Trimming "+filename+"\n")
      msa = AlignIO.read(inputpath+filename, 'fasta')

      if (gblockspars == "strict"):
            pars = '-t=c -p=t -b3=8  -b4=10 -b5=n'
      elif (gblockspars == "relaxed"): 
            n = math.floor(len(msa)/2+1)
            pars = '-t=c -p=t -b2=%d -b3=10 -b4=5  -b5=h' % n
      else:
          sys.stdout.write("\tUnknown Gblocks parameters! Exiting...\n")
          sys.exit()

      trimlength = runGblocks(inputpath+filename, pars)
      statfile.write("%s\t%d\t%d\t%.3f\n" % (filename, len(msa[0]), trimlength, float(trimlength)/len(msa[0])))

      shutil.move(inputpath+filename+'-gb', outputpath+filename[:filename.rfind('.')]+".gb_"+gblockspars+".fas")
      shutil.move(inputpath+filename+'-gb.txt', outputpath+filename[:filename.rfind('.')]+".gb_"+gblockspars+".txt")


statfile.close()

end = time.time()
print "Total time taken: "+str(end-start)+" seconds" 

