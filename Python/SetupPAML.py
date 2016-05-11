import sys, os, time, string, stat, shutil
import numpy as np
from fnmatch import fnmatch
from optparse import OptionParser




# Setup control files for different PAML models and runs
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
                  help = "Directory containing MSAs in phylip format (*.phylip) and trees [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory to write output to [default = %default]")

parser.add_option("-c","--ctlpath",
                  dest = "ctlpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory containing template PAML control files [default = %default]")

parser.add_option("-m","--msaprogram",
                  dest = "msaprogram",
                  default = "mafft",
                  metavar = "mafft/muscle/probcons/prank/prographmsa",
                  help = "MSA program to use [default = %default]")

parser.add_option("-s","--species",
                  dest = "species",
                  default = "6",
                  metavar = "integer",
                  help = "Number of species to use (6,5,4) [default = %default]")

parser.add_option("-r","--randseed",
                  dest = "randseed",
                  default = "-1",
                  metavar = "Int",
                  help = "Seed for random number generator [default = %default]")

parser.add_option("-R","--randstarts",
                  dest = "randstarts",
                  default = "0",
                  metavar = "Int",
                  help = "Number of random starting points for w [default = %default]")


(options,args) = parser.parse_args()

inputpath  = os.path.abspath(options.inputpath)
outputpath = os.path.abspath(options.outputpath)
ctlpath    = os.path.abspath(options.ctlpath)
msaprogram = options.msaprogram.lower()
species    = int(options.species)
randseed   = int(options.randseed)
randstarts = int(options.randstarts)


################################################################################################################################  

def MakeBranchSiteTrees6Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    mrotubr = parts[1]
    nvitrbr = parts[3]
    bimpabr = parts[5]
    bterrbr = parts[7]
    aflorbr = parts[10]
    amellbr = parts[12]

    apisbr    = parts[13]
    bombusbr  = parts[8]
    outlierbr = parts[14] 

    trees = []

    # Tree 1: 3 internal branches under selection
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr+" '#1' ", outlierbr+" '#1' ")
    trees.append(tree)

    # Tree 2: 2 internal branches between Apis and Bombus under selection
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr+" '#1' ", outlierbr)
    trees.append(tree)

    # Tree 3: Branch connecting Megachile and Nasonia under selection
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr, bterrbr, bombusbr, aflorbr, amellbr, apisbr, outlierbr+" '#1' ")
    trees.append(tree)

    # Tree 4: Only Bombus under selection
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr, outlierbr)
    trees.append(tree)

    treefile = open(outputpath+'/'+filename+'_branchsitetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#

def MakeBranchSiteTrees5Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    mrotubr = parts[1]
    bimpabr = parts[3]
    bterrbr = parts[5]
    aflorbr = parts[8]
    amellbr = parts[10]

    apisbr    = parts[11]
    bombusbr  = parts[6]

    trees = []

    # Tree 1: 3 internal branches under selection
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr+" '#1' ", bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr+" '#1' ")
    trees.append(tree)

    # Tree 2: 2 internal branches between Apis and Bombus under selection
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr, bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr+" '#1' ")
    trees.append(tree)

    # Tree 3: Branch connecting Megachile and Nasonia under selection
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr+" '#1' ", bimpabr, bterrbr, bombusbr, aflorbr, amellbr, apisbr)
    trees.append(tree)

    # Tree 4: Only Bombus under selection
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr, bimpabr, bterrbr, bombusbr+" '#1' ", aflorbr, amellbr, apisbr)
    trees.append(tree)

    # Tree 5: Only Apis under selection
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr, bimpabr, bterrbr, bombusbr, aflorbr, amellbr, apisbr+" '#1' ")
    trees.append(tree)


    treefile = open(outputpath+'/'+filename+'_branchsitetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#


def MakeBranchSiteTrees4Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    bimpabr = parts[1]
    bterrbr = parts[3]
    aflorbr = parts[5]
    amellbr = parts[7]
    interbr = parts[8]

    trees = []

    # Tree 1: Internal branch under selection
    tree = "(BIMPA:%s,BTERR:%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (bimpabr, bterrbr, aflorbr, amellbr, interbr+" '#1' ")
    trees.append(tree)

    treefile = open(outputpath+'/'+filename+'_branchsitetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#



def MakeCladeTrees6Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    mrotubr = parts[1]
    nvitrbr = parts[3]
    bimpabr = parts[5]
    bterrbr = parts[7]
    aflorbr = parts[10]
    amellbr = parts[12]

    apisbr    = parts[13]
    bombusbr  = parts[8]
    outlierbr = parts[14] 

    trees = []

    # Tree 1: 2 clades
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr+" '#1' ", bterrbr+" '#1' ", bombusbr+" '#1' ", aflorbr+" '#1' ", amellbr+" '#1' ", apisbr+" '#1' ", outlierbr)
    trees.append(tree)

    # Tree 2: 3 clades
    tree = "(MROTU:%s,NVITR:%s,((BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s):%s);\n" % \
            (mrotubr, nvitrbr, bimpabr+" '#1' ", bterrbr+" '#1' ", bombusbr+" '#1' ", aflorbr+" '#2' ", amellbr+" '#2' ", apisbr+" '#2' ", outlierbr)
    trees.append(tree)

    treefile = open(outputpath+'/'+filename+'_cladetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#

def MakeCladeTrees5Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    mrotubr = parts[1]
    bimpabr = parts[3]
    bterrbr = parts[5]
    aflorbr = parts[8]
    amellbr = parts[10]

    apisbr    = parts[11]
    bombusbr  = parts[6]

    trees = []

    # Tree 1: 2 clades
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr, bimpabr+" '#1' ", bterrbr+" '#1' ", bombusbr+" '#1' ", aflorbr+" '#1' ", amellbr+" '#1' ", apisbr+" '#1' ")
    trees.append(tree)

    # Tree 2: 3 clades
    tree = "(MROTU:%s,(BIMPA:%s,BTERR:%s):%s,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (mrotubr, bimpabr+" '#1' ", bterrbr+" '#1' ", bombusbr+" '#1' ", aflorbr+" '#2' ", amellbr+" '#2' ", apisbr+" '#2' ")
    trees.append(tree)

    treefile = open(outputpath+'/'+filename+'_cladetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#

def MakeCladeTrees4Species(filename):
    """Create trees with switches for PAML analyses for the group in question (for branch-site model)

    """

    # Read tree
    treefile = open(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt','r')
    tree = treefile.readline()
    treefile.close()

    parts       = tree[:-2].replace(',',':').replace(')','').split(':')
    bimpabr = parts[1]
    bterrbr = parts[3]
    aflorbr = parts[5]
    amellbr = parts[7]
    interbr = float(parts[8])

    trees = []

    # Tree 1: 3 internal branches under selection
    tree = "((BIMPA:%s,BTERR:%s):%.12f,(AFLOR:%s,AMELL:%s):%s);\n" % \
            (bimpabr, bterrbr, interbr/2, aflorbr+" '#1' ", amellbr+" '#1' ", "%.12f '#1' " % (interbr/2))
    trees.append(tree)

    treefile = open(outputpath+'/'+filename+'_cladetrees.txt','w')
    for tree in trees:
        treefile.write(tree+'\n')
    treefile.close()
#

def MakeBranchSiteTrees(filename):

    if (species == 6):
        MakeBranchSiteTrees6Species(filename)
    elif (species == 5):
        MakeBranchSiteTrees5Species(filename)
    elif (species == 4):
        MakeBranchSiteTrees4Species(filename)
    else:
        sys.stdout.write('Invalid number of species\n')
        sys.exit()
#
    

def MakeCladeTrees(filename):

    if (species == 6):
        MakeCladeTrees6Species(filename)
    elif (species == 5):
        MakeCladeTrees5Species(filename)
    elif (species == 4):
        MakeCladeTrees4Species(filename)
    else:
        sys.stdout.write('Invalid number of species\n')
        sys.exit()
#
    




def GetControlFiles(outputpath):

    ctlfiles = []
    for filename in os.listdir(outputpath):
        if (fnmatch(filename, '*.ctl')):
            ctlfiles.append(filename)
    return ctlfiles
#


def MakeSitesCTLFile(filename, ctlfile, winit):

    # Filename is '*.'+msaprogram+'.dna.phylip'
    group   = filename[:-11]
    runname = ctlfile[:ctlfile.rfind('.')]+'-%.3f' % winit
    runpath = outputpath+'/'+group+'_'+runname+'/'
    if (not os.path.exists(runpath)):
        os.mkdir(runpath)

    infile  = open(ctlpath+'/'+ctlfile,'r')
    outfile = open(runpath+group+'_'+runname+'.ctl','w')

    i = 0
    for line in infile:
        if (i == 0):
            outfile.write(line[:line.find('=')+1]+' ../'+filename+'\n')
        elif (i == 1):
            outfile.write(line[:line.find('=')+1]+' ../'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt\n')
        elif (i == 2):
            outfile.write(line[:line.find('=')+1]+' '+group+'_'+runname+'.out\n')
        elif (line[:line.find('=')].strip() == 'omega'):
            outfile.write(line[:line.find('=')+1]+' '+str(winit)+'   '+line[line.find('*'):])
        else:
            outfile.write(line.rstrip()+'\n')
        i += 1

    infile.close()
    outfile.close()

    script = open(outputpath+'/'+runname+'.sh','a')
    script.write('cd '+group+'_'+runname+'\ncodeml '+group+'_'+runname+'.ctl\ncd ..\n')
    script.close()

    brutusscript = open(outputpath+'/'+runname+'.brutus.sh','a')
    brutusscript.write('cd '+group+'_'+runname+"\nbsub -W%d:0 -R 'rusage[mem=1024]' -o %s" % (8, group+'_'+runname+'.brutus.out'))
    brutusscript.write(' ~/paml44/bin/codeml '+group+'_'+runname+'.ctl\ncd ..\n')
    brutusscript.close()
#


def MakeBranchCTLFile(filename, ctlfile, cladebranch, winit):

    # Filename is '*.'+msaprogram+'.dna.phylip'
    group   = filename[:-11]
    runname = ctlfile[:ctlfile.rfind('.')]+'-%.3f' % winit
    runpath = outputpath+'/'+group+'_'+runname+'/'
    if (not os.path.exists(runpath)):
        os.mkdir(runpath)

    infile  = open(ctlpath+'/'+ctlfile,'r')
    outfile = open(runpath+group+'_'+runname+'.ctl','w')

    i = 0
    for line in infile:
        if (i == 0):
            outfile.write(line[:line.find('=')+1]+' ../'+filename+'\n')
        elif (i == 1):
            if (cladebranch == 'branchsite'):
                outfile.write(line[:line.find('=')+1]+' ../'+filename+'_branchsitetrees.txt\n')
            else:
                outfile.write(line[:line.find('=')+1]+' ../'+filename+'_cladetrees.txt\n')
        elif (i == 2):
            outfile.write(line[:line.find('=')+1]+' '+group+'_'+runname+'.out\n')
        elif (line[:line.find('=')].strip() == 'omega'):
            outfile.write(line[:line.find('=')+1]+' '+str(winit)+'   '+line[line.find('*'):])
        else:
            outfile.write(line.rstrip()+'\n')
        i += 1

        
    infile.close()
    outfile.close()

    script = open(outputpath+'/'+runname+'.sh','a')
    script.write('cd '+group+'_'+runname+'\ncodeml '+group+'_'+runname+'.ctl\ncd ..\n')
    script.close()

    brutusscript = open(outputpath+'/'+runname+'.brutus.sh','a')
    brutusscript.write('cd '+group+'_'+runname+"\nbsub -W%d:0 -R 'rusage[mem=1024]' -o %s" % (1, group+'_'+runname+'.brutus.out'))
    brutusscript.write(' ~/paml44/bin/codeml '+group+'_'+runname+'.ctl\ncd ..\n')
    brutusscript.close()
#





################################################################################################################################  

start = time.time()

#ctlfiles = GetControlFiles(outputpath)
if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

if (randseed > 0):
    sys.stdout.write("Seeded random number generator with %d\n" % randseed)
    np.random.seed(randseed)
#
randomegas = np.round(np.random.rand(randstarts)*4, 2)

wvals = [0.5, 1, 2]
if (randstarts > 0):
    wvals = wvals+list(randomegas)



# Do different starting values for models where omega is not fixed (some may get stuck in local optima)
for filename in os.listdir(inputpath):
      if (fnmatch(filename,'*.'+msaprogram+'.dna.phylip')):
            sys.stdout.write(filename+'\n')
            shutil.copy(inputpath+'/'+filename, outputpath+'/'+filename)
            shutil.copy(inputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt', 
                        outputpath+'/'+filename[:filename.find('dna')]+'phylip_phyml_tree.txt')

            MakeBranchSiteTrees(filename)
            MakeCladeTrees(filename)

            MakeSitesCTLFile(filename,  'M8A.ctl', 1)
            MakeBranchCTLFile(filename, 'BranchSiteNull.ctl', 'branchsite',1)
            for w in wvals:
                MakeSitesCTLFile(filename,  'Sites.ctl', w)    
                MakeBranchCTLFile(filename, 'BranchSiteAlt.ctl',  'branchsite', w)
                MakeBranchCTLFile(filename, 'CladeC.ctl', 'clade', w)
                MakeBranchCTLFile(filename, 'CladeD.ctl', 'clade', w)
            #
      #
#

end = time.time()
print "Total time taken: "+str(end-start)+" seconds"