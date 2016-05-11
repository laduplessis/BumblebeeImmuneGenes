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
                  help = "Directory containing orthodb alignments in phylip format (*.phylip) and trees (from ExtractSequences.py) [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing template PAML control files (and where to write new ones) [default = %default]")

parser.add_option("-m","--model",
                  dest = "model",
                  default = "Sites,M8A,Branchsite,Clade",
                  metavar = "",
                  help = "(sites/branchsite/clade) or comma separated list [default = %default]")


(options,args) = parser.parse_args()

inputpath   = os.path.abspath(options.inputpath)
outputpath  = os.path.abspath(options.outputpath)
models      = options.model.lower().split(',')





################################################################################################################################  





def SitesModel(pattern):

      # Create files (for overwriting)
      M0file = open(outputpath+'/'+pattern+'.M0.w','w')
      M1file = open(outputpath+'/'+pattern+'.M1.w','w')
      M2file = open(outputpath+'/'+pattern+'.M2.w','w')
      M3file = open(outputpath+'/'+pattern+'.M3.w','w')
      M7file = open(outputpath+'/'+pattern+'.M7.w','w')
      M8file = open(outputpath+'/'+pattern+'.M8.w','w')

      M0file.write('Alignment\tLk\tw0\tp0\tdNlen\tdSlen\n')
      M1file.write('Alignment\tLk\tw0\tw1\tp0\tp1\n')
      M2file.write('Alignment\tLk\tw0\tw1\tw2\tp0\tp1\tp2\n')
      M3file.write('Alignment\tLk\tw0\tw1\tw2\tp0\tp1\tp2\n')
      M7file.write('Alignment\tLk\tw0\tw1\tw2\tw3\tw4\tw5\tw6\tw7\tw8\tw9\tp0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\n')
      M8file.write('Alignment\tLk\tw0\tw1\tw2\tw3\tw4\tw5\tw6\tw7\tw8\tw9\tw10\tp0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\tp10\n')

      M0file.close()
      M1file.close()
      M2file.close()
      M3file.close()
      M7file.close()
      M8file.close()

      for filename in os.listdir(inputpath):
            if (fnmatch(filename,'*_'+pattern)):
                  name = filename[:filename.find('.')]
                  sys.stdout.write(name+'\n')

                  pamlfile = open(inputpath+'/'+filename+'/'+filename+'.out','r')

                  model = -1
                  for line in pamlfile:
                        #print line
                        if (len(line) >= 6 and line[:5] == 'Model' and line[5] != ':'):
                              if (model >= 0):
                                    outfile = open(outputpath+'/'+pattern+'.M%d.w' % model,'a')
                                    outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                                    for i in w:
                                          outfile.write('\t%8.5f' % i)
                                    for i in p:
                                          outfile.write('\t%8.5f' % i)
                                    if (dNlen >= 0):
                                          outfile.write('\t%8.5f' % dNlen)
                                    if (dSlen >= 0): 
                                          outfile.write('\t%8.5f' % dSlen)
                                    outfile.write('\n')
                                    outfile.close()
                              model = int(line[5:line.find(':')])
                              dNlen = -1
                              dSlen = -1

                        if (model == 0):
                              if (len(line) >= 3 and line[:3] == 'lnL'):
                                    lnL = float(line.split()[4])
                              if (len(line) >= 5 and line[:5] == 'omega'):
                                    w = [float(line[line.find('=')+1:])]
                                    p = [1]
                              if (len(line) >= 5 and line[:18] == 'tree length for dN'):
                                    dNlen = float(line[line.find(':')+1:]) 
                              if (len(line) >= 5 and line[:18] == 'tree length for dS'):
                                    dSlen = float(line[line.find(':')+1:]) 
                        elif (model > 0):
                              if (len(line) >= 3 and line[:3] == 'lnL'):
                                    lnL = float(line.split()[4])
                              if (len(line) >= 2 and line[:2] == 'p:'):
                                    p = map(float, line.split()[1:])
                              if (len(line) >= 2 and line[:2] == 'w:'):
                                    w = map(float, line.split()[1:])
                        else:
                              pass
                  pamlfile.close()
                  if (model >= 0):
                        outfile = open(outputpath+'/'+pattern+'.M%d.w' % model,'a')
                        outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                        for i in w:
                              outfile.write('\t%8.5f' % i)
                        for i in p:
                              outfile.write('\t%8.5f' % i)
                        if (dNlen >= 0):
                              outfile.write('\t%8.5f' % dNlen)
                        if (dSlen >= 0): 
                              outfile.write('\t%8.5f' % dSlen)
                        outfile.write('\n')
                        outfile.close()
      #
#


def M8AModel():

      # Create files (for overwriting)
      outfile = open(outputpath+'/Sites.M8A-1.000.w','w')
      outfile.write('Alignment\tLk\tw0\tw1\tw2\tw3\tw4\tw5\tw6\tw7\tw8\tw9\tw10\tp0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\tp10\n')

      for filename in os.listdir(inputpath):
            if (fnmatch(filename,'*_M8A-1.000')):
                  name = filename[:filename.find('.')]
                  sys.stdout.write(name+'\n')

                  pamlfile = open(inputpath+'/'+filename+'/'+filename+'.out','r')

                  model = -1
                  for line in pamlfile:

                        if (len(line) >= 3 and line[:3] == 'lnL'):
                              lnL = float(line.split()[4])
                        if (len(line) >= 2 and line[:2] == 'p:'):
                              p = map(float, line.split()[1:])
                        if (len(line) >= 2 and line[:2] == 'w:'):
                              w = map(float, line.split()[1:])
                  pamlfile.close()
                  outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                  for i in w:
                        outfile.write('\t%8.5f' % i)
                  for i in p:
                        outfile.write('\t%8.5f' % i)
                  outfile.write('\n')
      outfile.close()         
#

def BranchSite(pattern):

      trees = []
      for filename in os.listdir(inputpath):
            if (fnmatch(filename,'*_'+pattern)):
                  name = filename[:filename.find('.')]
                  sys.stdout.write(name+'\n')

                  pamlfile = open(inputpath+'/'+filename+'/'+filename+'.out','r')

                  tree = 0
                  for line in pamlfile:
                        if (len(line) >= 4 and line[:4] == 'TREE'):
                              if (tree > 0):
                                    if (tree in trees):
                                        outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','a')
                                    else:
                                        outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','w')
                                        outfile.write('Alignment\tLk\tw0\tw1\tw2\tp0\tp1\tp2a\tp2b\n')
                                        trees.append(tree)
                                    outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                                    outfile.write('\t%8.5f\t%8.5f\t%8.5f' % (bw[0], bw[1], fw[2]))
                                    for i in p:
                                          outfile.write('\t%8.5f' % i)
                                    outfile.write('\n')
                                    outfile.close()
                              tree += 1
                        if (len(line) >= 3 and line[:3] == 'lnL'):
                              lnL = float(line.split()[4])
                        if (len(line) >= 10 and line[:10] == 'proportion'):
                              p  = map(float, line.split()[1:])
                        if (len(line) >= 10 and line[:10] == 'background'):
                              bw = map(float, line.split()[2:])
                        if (len(line) >= 10 and line[:10] == 'foreground'):
                              fw = map(float, line.split()[2:])
                        
                  pamlfile.close()
                  if (tree > 0):
                        if (tree in trees):
                            outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','a')
                        else:
                            outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','w')
                            outfile.write('Alignment\tLk\tw0\tw1\tw2\tp0\tp1\tp2a\tp2b\n')
                            trees.append(tree)
                        outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                        outfile.write('\t%8.5f\t%8.5f\t%8.5f' % (bw[0], bw[1], fw[2]))
                        for i in p:
                              outfile.write('\t%8.5f' % i)
                        outfile.write('\n')
                        outfile.close()

      #
#


def Clade(pattern):

      trees = []
      for filename in os.listdir(inputpath):
            if (fnmatch(filename,'*_'+pattern)):
                  name = filename[:filename.find('.')]
                  sys.stdout.write(name+'\n')

                  pamlfile = open(inputpath+'/'+filename+'/'+filename+'.out','r')

                  tree = 0
                  for line in pamlfile:
                        if (len(line) >= 4 and line[:4] == 'TREE'):
                              if (tree > 0):
                                    if (tree in trees):
                                        outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','a')
                                    else:
                                        outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','w')
                                        outfile.write('Alignment\tLk\tw0\tw1')
                                        for i in range(0,len(w)):
                                                outfile.write('\tw%d' % (i+2))
                                        outfile.write('\tp0\tp1\tp2\n')
                                        trees.append(tree)

                                    outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                                    outfile.write('\t%8.5f\t%8.5f' % (w[0][0], w[0][1]))
                                    for i in w:
                                          outfile.write('\t%8.5f' % i[2])
                                    for i in p:
                                          outfile.write('\t%8.5f' % i)
                                    outfile.write('\n')
                                    outfile.close()
                              tree += 1
                              w = []

                        if (len(line) >= 3 and line[:3] == 'lnL'):
                              lnL = float(line.split()[4])
                        if (len(line) >= 10 and line[:10] == 'proportion'):
                              p  = map(float, line.split()[1:])
                        if (len(line) >= 10 and line[:10] == 'branch typ'):
                              w.append(map(float, line.split()[3:]))
                        
                  pamlfile.close()
                  if (tree > 0):
                        if (tree in trees):
                            outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','a')
                        else:
                            outfile = open(outputpath+'/'+pattern+'.Tree'+str(tree)+'.w','w')
                            outfile.write('Alignment\tLk\tw0\tw1')
                            for i in range(0,len(w)):
                                    outfile.write('\tw%d' % (i+2))
                            outfile.write('\tp0\tp1\tp2\n')
                            trees.append(tree)

                        outfile.write(name.ljust(20)+'\t%.5f' % lnL)
                        outfile.write('\t%8.5f\t%8.5f' % (w[0][0], w[0][1]))
                        for i in w:
                              outfile.write('\t%8.5f' % i[2])
                        for i in p:
                              outfile.write('\t%8.5f' % i)
                        outfile.write('\n')
                        outfile.close()
      #
#

def GetOmegaVals(pattern):

      for filename in os.listdir(inputpath):
            if (fnmatch(filename,'*.dna.phylip')):
                  group = filename[:-11]
                  break

      print group

      wvals = []
      for filename in os.listdir(inputpath):
            if (fnmatch(filename, group+'_'+pattern+'*')):
                  wvals.append(filename[filename.find('-')+1:])


      print wvals
      return wvals
#



################################################################################################################################  

start = time.time()

if (not os.path.exists(outputpath)):
      os.mkdir(outputpath)


for model in models:

      if (model == 'sites'):
            wvals = GetOmegaVals('Sites')
            for w in wvals:
                  SitesModel('Sites-'+w)
      elif (model == 'm8a'):
            M8AModel()
      elif (model == 'branchsite'):
            wvals = GetOmegaVals('BranchSiteAlt')
            for w in wvals:
                  BranchSite('BranchSiteAlt-'+w)
            BranchSite('BranchSiteNull-1.000')
      elif (model == 'clade'):
            wvals = GetOmegaVals('CladeC')
            for w in wvals:
                Clade('CladeC-'+w)
                Clade('CladeD-'+w)
      #
#

end = time.time()
print "Total time taken: "+str(end-start)+" seconds"