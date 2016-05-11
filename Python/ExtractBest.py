import sys, os, time, string, stat, shutil

from scipy.stats import chi2
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

parser.add_option("-m","--model",
                  dest = "model",
                  default = "Sites,Branchsite,Clade",
                  metavar = "",
                  help = "(sites/branchsite/clade) or comma separated list [default = %default]")

#parser.add_option("-w","--wmax",
#                  dest = "wmax",
#                  default = "1000",
#                  metavar = "",
#                  help = "Maximum value for w, disregard runs where w bigger than this value [default = %default]")

(options,args) = parser.parse_args()

inputpath   = os.path.abspath(options.inputpath)
models      = options.model.lower().split(',')
#wmax        = float(options.wmax)





################################################################################################################################  

def GetBest(pattern):

      runslk   = []
      runsdata = []
      runids   = []

      for filename in os.listdir(inputpath):
            if (fnmatch(filename,pattern) and 'Best' not in filename and 'Skipped' not in filename):
                  sys.stdout.write("Extracting best for "+filename+'\n')

                  data    = []
                  lkhoods = []

                  lkfile = open(inputpath+'/'+filename, 'r')
                  header = lkfile.readline()
                  for line in lkfile:
                        data.append(line)
                        #wpvals = map(float, line.split()[2:])
                        #if (max(wpvals) >= wmax):
                        #      sys.stdout.write('\t%s: w out of range!\n' % line[:line.find('\t')])
                        #      lkhoods.append(-999999)
                        #else:
                        lkhoods.append(float(line.split()[1]))
                  lkfile.close()
                  runslk.append(lkhoods)
                  runsdata.append(data)
                  runids.append(filename[filename.find('-')+1:filename.find('.',filename.find('-')+3)])
            #
      #

      if (len(runslk) == 0):
            return
      
      outfile  = open(inputpath+'/'+pattern.replace('*','-Best'),'w')
      #skipfile = open(inputpath+'/'+pattern.replace('*','-Skipped'),'w')
      outfile.write(header.strip()+'\tBest\n')
      #skipfile.write(header)
      for i in range(0,len(runslk[0])):
            maxidx = 0
            maxlk  = runslk[0][i]
            sys.stdout.write('%d\t%.3f' % (i, maxlk))
            for j in range(1,len(runslk)):
                  sys.stdout.write("\t%.3f" % runslk[j][i])
                  if (runslk[j][i] > maxlk):
                        maxlk  = runslk[j][i]
                        maxidx = j
            sys.stdout.write('\t - %.3f\t%d' % (maxlk, maxidx))

            # Write output
            #if (maxlk == -999999):
            #      sys.stdout.write('\t- w out of range! DELETED!\n')
            #      skipfile.write(runsdata[maxidx][i])
            #else: 
            sys.stdout.write('\n')
            outfile.write(runsdata[maxidx][i].strip()+'\t'+runids[maxidx]+'\n')
      #
      outfile.close()
      #skipfile.close()

##




################################################################################################################################  

start = time.time()


for model in models:

      if (model == 'sites'):
            GetBest('Sites*.M0.w')
            GetBest('Sites*.M1.w')
            GetBest('Sites*.M2.w')
            GetBest('Sites*.M3.w')
            GetBest('Sites*.M7.w')
            GetBest('Sites*.M8.w')
      elif (model == 'branchsite'):
            GetBest('BranchSiteAlt*.Tree1.w')
            GetBest('BranchSiteAlt*.Tree2.w')
            GetBest('BranchSiteAlt*.Tree3.w')
            GetBest('BranchSiteAlt*.Tree4.w')
            GetBest('BranchSiteAlt*.Tree5.w')
      elif (model == 'clade'):
            GetBest('CladeC*.Tree1.w')
            GetBest('CladeC*.Tree2.w')
            GetBest('CladeD*.Tree1.w')
            GetBest('CladeD*.Tree2.w')

      #
#

end = time.time()
print "Total time taken: "+str(end-start)+" seconds"