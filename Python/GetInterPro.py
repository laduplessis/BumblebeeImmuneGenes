import base64, platform, os, suds, sys, time, urllib2
import logging
from suds import WebFault
from suds.client import Client

from Bio import SeqIO
from fnmatch import fnmatch
from optparse import OptionParser



# Get results of InterPro Scan for AA sequences 
#
# A lot of functions based on http://www.ebi.ac.uk/Tools/webservices/download_clients/python/suds/iprscan5_suds.py
#
#
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
                  help = "Directory containing result  [default = %default]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "filename",
                  help = "Directory to store interpro results  [default = %default]")

#parser.add_option("-c","--crc",
#                  dest = "crc",
#                  action = "store_true",
#                  help = "Enable InterProScan Matches look-up (faster)  [default = %default]")

parser.add_option("-g","--goterms",
                  dest = "goterms",
                  action = "store_true",
                  help = "Enable inclusion of GO terms  [default = %default]")

parser.add_option("-e","--email",
                  dest = "email",
                  default = "",
                  metavar = "filename",
                  help = "email address  [default = %default]")





(options,args) = parser.parse_args()

wsdlUrl     = 'http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl'
inputpath   = os.path.abspath(options.inputpath)
outputpath  = os.path.abspath(options.outputpath)
#crc         = options.crc
goterms     = options.goterms
email       = options.email




################################################################################################################################

# Debug level
debugLevel = 0

# Debug print
def printDebugMessage(functionName, message, level):
    if(level <= debugLevel):
        print >>sys.stderr, '[' + functionName + '] ' + message

# Submit job
def serviceRun(email, title, params):
    printDebugMessage('serviceRun', 'Begin', 1)
    jobid = server.run(email=email, title=title, parameters=params)
    printDebugMessage('serviceRun', 'End', 1)
    return jobid

# Get job status
def serviceCheckStatus(jobId):
    printDebugMessage('serviceCheckStatus', 'jobId: ' + jobId, 1)
    result = server.getStatus(jobId = jobId)
    return result

# Get available result types for job
def serviceGetResultTypes(jobId):
    printDebugMessage('serviceGetResultTypes', 'Begin', 1)
    result = server.getResultTypes(jobId=jobId)
    printDebugMessage('serviceGetResultTypes', 'End', 1)
    return result['type']

# Get result
def serviceGetResult(jobId, type):
    printDebugMessage('serviceGetResult', 'Begin', 1)
    printDebugMessage('serviceGetResult', 'jobId: ' + jobId, 1)
    printDebugMessage('serviceGetResult', 'type: ' + type, 1)
    resultBase64 = server.getResult(jobId=jobId, type=type)
    result = base64.decodestring(resultBase64)
    printDebugMessage('serviceGetResult', 'End', 1)
    return result    


# Client-side poll
def clientPoll(jobId):
    printDebugMessage('clientPoll', 'Begin', 1)
    result = 'PENDING'
    sys.stdout.write("\t\t")
    while result == 'RUNNING' or result == 'PENDING':
        result = serviceCheckStatus(jobId)
        sys.stdout.write(result+ "...")
        if result == 'RUNNING' or result == 'PENDING':
            time.sleep(15)
    sys.stdout.write("\n")
    printDebugMessage('clientPoll', 'End', 1)


# Get result for a jobid
def getResult(jobId, group, species):
    printDebugMessage('getResult', 'Begin', 1)
    printDebugMessage('getResult', 'jobId: ' + jobId, 1)

    # Check status and wait if necessary
    clientPoll(jobId)
    # Get available result types
    resultTypes = serviceGetResultTypes(jobId)
    for resultType in resultTypes:
		# Get the result
		result = serviceGetResult(jobId, resultType['identifier'])
		# Derive the filename for the result
		filename = outputpath + '/' + group + '.'+ species + '.'+ resultType['identifier'] + '.'+ resultType['fileSuffix']



		# Write a result file
		fh = open(filename, 'w');
		fh.write(result)
		fh.close()
		#sys.stdout.write('\t\t'+filename+'\n')

    printDebugMessage('getResult', 'End', 1)



################################################################################################################################

start = time.time()

# Set up connection
client = Client(wsdlUrl)
server = client.service

# Set the client user-agent.
clientRevision = '$Revision: 2762 $'
clientVersion = '0'
if len(clientRevision) > 11:
    clientVersion = clientRevision[11:-2] 
userAgent = 'EBI-Sample-Client/%s (%s; Python %s; %s) suds/%s Python-urllib/%s' % (
    clientVersion, os.path.basename( __file__ ),
    platform.python_version(), platform.system(),
    suds.__version__, urllib2.__version__
)
sys.stderr.write('userAgent: ' + userAgent+'\n')
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)

# Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
if os.environ.has_key('http_proxy'):
    http_proxy_conf = os.environ['http_proxy'].replace('http://', '')
    proxy = dict(http=http_proxy_conf)
    client.set_options(proxy=proxy)
elif os.environ.has_key('HTTP_PROXY'):
    http_proxy_conf = os.environ['HTTP_PROXY'].replace('http://', '')
    proxy = dict(http=http_proxy_conf)
    client.set_options(proxy=proxy)


# Set global parameters
params = {} 
#if crc:
#	params['nocrc'] = 0
#else:#
#	params['nocrc'] = 1
   
if goterms:
    params['goterms'] = 1
else:
    params['goterms'] = 0
# Add the other options (if defined)
#if options.appl:
#    params['appl'] = {'string':options.appl}



if (not os.path.exists(outputpath)):
	os.mkdir(outputpath)

for filename in sorted(os.listdir(inputpath), reverse=True):
    if (fnmatch(filename,'*.fa')):
  		sys.stdout.write(filename[:-3]+'\n')
  		group = filename[:-3]

  		for seq in SeqIO.parse(inputpath+'/'+filename,'fasta'):
			sys.stdout.write('\t'+seq.id+'\t'+str(seq.seq[:25])+'...  \t')

			species  = seq.id 
			params['sequence'] = str(seq.seq)
			jobid = serviceRun(email, group+'_'+species, params)
			sys.stdout.write("JobID="+jobid+"\n")

			time.sleep(5)
			getResult(jobid, group, species)
  	#
#


end = time.time()
print "Total time taken: "+str(end-start)+" seconds"

