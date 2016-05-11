import sys, os, time, string, stat, shutil

#from Bio import AlignIO
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
                  help = "Directory containing the BEB and InterPro results [default = %default]")

parser.add_option("-d","--databases",
                  dest = "databases",
                  default = "SMART,Pfam",
                  metavar = "Comma separated list",
                  help = "Comma separated list of dbs to extract domains from [default = %default]")


(options,args) = parser.parse_args()

inputpath = os.path.abspath(options.inputpath)
databases = map(string.lower, options.databases.split(','))

################################################################################################################################  


def GetIPRDomains(fname, databases):

      domains = []
      iprfile = open(fname, 'r')
      for line in iprfile:
            parts = line.split('\t')
            if (parts[3].lower() in databases):
                  domain = dict()
                  domain['name']  = parts[5]
                  domain['db']    = parts[3]
                  domain['id']    = parts[4]
                  domain['start'] = parts[6]
                  domain['end']   = parts[7]
                  domains.append(domain)

      iprfile.close()
      return domains
#

################################################################################################################################  

start = time.time()

for filename in os.listdir(inputpath+'/Raw/'):
      if (fnmatch(filename,'*.tsv.*')):
            (group, org, a, b) = filename.split('.')

            domains = GetIPRDomains(inputpath+'/Raw/'+filename, databases)
            domfile = open(inputpath+'/'+group+'.'+org+'.domains','w')
            domfile.write('\t'.join(['Name', 'Database', 'Id', 'Start', 'End'])+'\n')
            for d in domains:
                  domfile.write('\t'.join(['"'+d['name']+'"', d['db'], d['id'], d['start'], d['end']])+'\n')
            domfile.close()

            # bebfile = open(inputpath+'/'+'.'.join([group, model, 'BEB', org]), 'r')
            # outfile = open(inputpath+'/'+'.'.join([group, model, org, 'IPR', 'BEB']),'w')
            # header = bebfile.readline()
            # outfile.write(header.strip()+'\tDomain\n')
            # for line in bebfile:
            #       site = int(line[:line.find('\t')])

            #       seqdoms = []
            #       for d in domains:
                       
            #             if (site >= int(d['start']) and site <= int(d['end'])):
            #                   seqdoms.append(d['name'])

            #       if (seqdoms == []):
            #             seqdoms = ['None']
            #       outfile.write(line.strip()+'\t"'+','.join(seqdoms)+'"\n')
            # bebfile.close()
            # outfile.close()
      #
#


end = time.time()
print "Total time taken: "+str(end-start)+" seconds"
