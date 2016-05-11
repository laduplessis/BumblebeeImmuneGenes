import sys, os, time, string, stat, shutil

from fnmatch import fnmatch
from optparse import OptionParser



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
                  help = "Directory containing result from ExtractParameters.py [default = %default]")

parser.add_option("-b","--bestpath",
                  dest = "bestpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing result from ExtractBest.py [default = %default]")

#parser.add_option("-m","--model",
#                  dest = "model",
#                  default = "M8,Branchsite",
#                  metavar = "",
#                  help = "(sites/branchsite/clade) or comma separated list [default = %default]")

parser.add_option("-m","--msaprogram",
                  dest = "msaprogram",
                  default = "mafft",
                  metavar = "mafft/muscle/probcons/prank/prographmsa",
                  help = "MSA program to use [default = %default]")

#parser.add_option("-w","--wmax",
#                  dest = "wmax",
#                  default = "1000",
#                  metavar = "",
#                  help = "Maximum value for w, disregard runs where w bigger than this value [default = %default]")

(options,args) = parser.parse_args()

inputpath   = os.path.abspath(options.inputpath)
bestpath    = os.path.abspath(options.bestpath)
msaprogram = options.msaprogram.lower()

models = ['m8','branchsite']
#models      = options.model.lower().split(',')
#wmax        = float(options.wmax)


################################################################################################################################  

def GetBEBM8(infname, outfname, number):

	infile = open(infname, 'r')
	i = 0
	step1 = False
	step2 = False
	for line in infile:
		if ('BEB' in line):
			i += 1
			if (i == number):
				step1   = True
				classes = int(line[line.find('for')+3:line.find('classes')])
				
				outfile = open(outfname, 'w')
				outfile.write('Site\tAA\t')
				for j in range(0,classes):
					outfile.write('p%d\t' % j)
				outfile.write('mean_w\tstd_w\n')

				continue

		if (step1 and line.strip() != ""):
			parts = line.split()
			try:
				site = int(parts[0])
				step2 = True

				outfile.write('%d\t%s\t' %(site, parts[1]))
				for j in range(0,classes):
					outfile.write('%s\t' % parts[j+2])
				outfile.write(parts[-3]+'\t'+parts[-1]+'\n')	
				
			except:
				if (step2):
					break

	infile.close()
#

def GetBEBBranchSite(infname, outfname, number):

	infile = open(infname, 'r')
	i = 0
	step1 = False
	step2 = False
	for line in infile:
		if ('BEB' in line):
			i += 1
			if (i == number):
				step1   = True
				classes = int(line[line.find('for')+3:line.find('classes')])
				
				outfile = open(outfname, 'w')
				outfile.write('Site\tAA\t')
				for j in range(0,classes):
					outfile.write('p%d\t' % j)
				outfile.write('\n')

				continue

		if (step1 and line.strip() != ""):
			parts = line.split()
			try:
				site = int(parts[0])
				step2 = True

				outfile.write('%d\t%s\t' %(site, parts[1]))
				for j in range(0,classes):
					outfile.write('%s\t' % parts[j+2])
				outfile.write('\n')
				
			except:
				if (step2):
					break

	infile.close()
#


################################################################################################################################  

start = time.time()

outputpath = bestpath+'/BEB'
if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

for model in models:

      if (model == 'm8'):
      	bestfile = open(bestpath+'/Sites-Best.M8.w')
      	bestfile.readline()
      	for line in bestfile:
      		parts = line.split()
      		GetBEBM8(inputpath+'/'+parts[0]+'.'+msaprogram+'_Sites-'+parts[-1]+'/rst', outputpath+'/'+parts[0]+'.M8.BEB',2)
      	bestfile.close()

      elif (model == 'branchsite'):

	      for filename in os.listdir(bestpath):
	            if (fnmatch(filename,'BranchSiteAlt-Best*')):
					nr = int(filename[filename.find('Tree')+4:filename.rfind('.')])

					bestfile = open(bestpath+'/'+filename)
					bestfile.readline()
					for line in bestfile:
						parts = line.split()
						GetBEBBranchSite(inputpath+'/'+parts[0]+'.'+msaprogram+'_BranchSiteAlt-'+parts[-1]+'/rst', outputpath+'/'+parts[0]+'.BS.Tree%d.BEB' % nr,nr)
					bestfile.close()

      #


end = time.time()
print "Total time taken: "+str(end-start)+" seconds"
