import sys, os
from fnmatch import fnmatch
from optparse import OptionParser
from Bio import SeqIO

#
# Extract all alignments that contain sequences from all the species in the organism list
# Only sequences from species in the list are extracted
# Restricted alignments are saved in a new directory in fasta format
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

parser.add_option("-l","--orglist",
                  dest = "orglist",
                  default = "",
                  metavar = "path",
                  help = "path to file containing list of organisms (lines starting with # skipped [optional]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path + filename",
                  help = "filename to store output in [required]")

parser.add_option("-p","--pattern",
                  dest = "pattern",
                  default = "*.aln",
                  metavar = "filename",
                  help = "Pattern for alignment filenames to search for (in quotes) [default = %default]")

(options,args) = parser.parse_args()

inputpath  = os.path.abspath(options.inputpath)
outputpath = os.path.abspath(options.outputpath)
orglist    = os.path.abspath(options.orglist)
pattern    = options.pattern
removegaps = True

################################################################################################################################    

# Load organism list if specified
orgfile = open(orglist,'r')
organisms = dict()
for line in orgfile:
    if (not line.isspace()):
        line = line.strip()        
        if (line[0] != '#'):
          organisms[line] = []
sys.stdout.write("Organisms: "+str(organisms.keys())+"\n")
##

# Read alignments
allfiles = dict()
sys.stdout.write("Reading alignments...\n")
for filename in os.listdir(inputpath):
    if (fnmatch(filename,pattern)):
        allfiles[filename] = 0
        alnfile = open(inputpath+'/'+filename,'r')
        for record in SeqIO.parse(alnfile,'fasta'):
            for org in organisms:
                if (org in record.description):
                    organisms[org].append(filename)
        alnfile.close()
##

# Print output and find which files contain all organisms
extractfiles = []
for org in organisms:
    sys.stdout.write(org+"\t%d genes\n" % len(organisms[org]))
sys.stdout.write('\nAlignment\t')
for org in organisms:
    sys.stdout.write('\t'+org)
sys.stdout.write('\n')

for filename in sorted(allfiles):
    sys.stdout.write(filename)
    for org in organisms:
        if (filename in organisms[org]):
            sys.stdout.write('\t%d' % organisms[org].count(filename))
            allfiles[filename] += 1
        else:
            sys.stdout.write('\t-')
    if (allfiles[filename] > 0):
        sys.stdout.write('\textracted!')
        extractfiles.append(filename)
    sys.stdout.write('\n')
##

# Extract alignments containing all organisms and only extract sequences
# or organisms in organism list, and remove gaps from alignments
# (will realign anyway)
sys.stdout.write("\nExtracting alignments... (%d alignments)\n" % len(extractfiles))
if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)
for filename in extractfiles:
    infile  = open(inputpath+'/'+filename,'r')
    outfile = open(outputpath+'/'+filename,'w')
    for record in SeqIO.parse(infile,'fasta'):
        for org in organisms:
            if (org in record.description):

                parts = record.description.split()
                record.id = parts[-1]+'_'+parts[0]
                record.description = parts[1]
 
                SeqIO.write(record, outfile, 'fasta')
    infile.close()
    outfile.close()
##    




