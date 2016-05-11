import sys, os, time, string, stat, shutil
from fnmatch import fnmatch
from subprocess import Popen, PIPE
from optparse import OptionParser




# Run PhyML to get branch length estimates (topology not optimized, uses JTT model)
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
                  help = "Directory containing orthodb alignments in phylip format (*.phy) (from ExtractSequences.py) [default = %default]")

parser.add_option("-t","--treefile",
                  dest = "treefile",
                  default = "",
                  metavar = "filename",
                  help = "Starting tree topology [default = %default]")

parser.add_option("-m","--msaprogram",
                  dest = "msaprogram",
                  default = "mafft",
                  metavar = "mafft/muscle/probcons/prank/prographmsa",
                  help = "MSA program to use [default = %default]")

parser.add_option("-p","--phyml",
                  dest = "phyml",
                  default = "phyml",
                  metavar = "filename",
                  help = "File name of the phyml binaries [default = %default]")

(options,args) = parser.parse_args()

inputpath  = os.path.abspath(options.inputpath)
treefile   = options.treefile
msaprogram = options.msaprogram.lower()
phyml      = options.phyml

################################################################################################################################  

def Phyml(inputfile):
 
  print "Starting PhyML computation for group "+inputfile
  start = time.time()

  # Run codonPhyML
  if (treefile == ""):
      callstring = "%s -i %s -d aa -m LG -b 0 -o tlr --quiet" % (phyml, inputpath+'/'+inputfile)
  else:
      callstring = "%s -i %s -d aa -m LG -b 0 -o lr --quiet -u %s " % (phyml, inputpath+'/'+inputfile, os.path.abspath(treefile))
  
  print callstring
  handle = Popen(callstring, stdout=None, stderr=PIPE, shell=True)
        
  err = handle.stderr.read()
  if (err != ""):
        sys.stderr.write("Warning! Errors encountered!\n")    
        sys.stderr.write(err)

  end = time.time()
  print "PhyML finished in ("+str(end-start)+" seconds)"
  print

##


################################################################################################################################  

start = time.time()

for filename in os.listdir(inputpath):
      if (fnmatch(filename,'*.'+msaprogram+'.phylip')):
            Phyml(filename)
      #
#

end = time.time()
print "Total time taken: "+str(end-start)+" seconds"


