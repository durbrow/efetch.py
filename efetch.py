#!python
"""
    Retrieve sequence data as FASTA from NCBI using eutils
    
    The following example will fetch the sequence of chromosome 8 of GRCh37:
        python efetch.py fasta CM000670.1

    See:
        http://www.ncbi.nlm.nih.gov/books/NBK25499/
        http://www.ncbi.nlm.nih.gov/books/NBK25500/
        http://www.ncbi.nlm.nih.gov/books/NBK25501/
"""

import httplib
import string
import json
import sys
import xml.parsers.expat

def _query(function, params):
    if 'db' not in params: params['db'] = 'nuccore' # set default database
    return '/entrez/eutils/{}.fcgi?'.format(function) + '&'.join(['{}={}'.format(k, v) for k, v in params.iteritems()])

def _lines(resp):
    buffer = resp.read(32768)
    at = 0
    while len(buffer) != 0:
        next = buffer.find('\n', at)
        if next < 0:
            more = resp.read(32768)
            if len(more) == 0:
                rslt = buffer[at:]
                buffer = ''
                at = 0
                yield rslt
            else:
                buffer = buffer[at:] + more
                at = 0
        else:
            rslt = buffer[at:next]
            at = next + 1
            yield rslt


class SRARun:
    def __init__(self, properties):
        self.properties = properties
    
    def accession(self):
        return str(self.properties['acc'])
    
    def totalSpots(self):
        return long(self.properties['total_spots'])
    
    def totalBases(self):
        return long(self.properties['total_bases'])
    
    def isPublic(self):
        return self.properties['load_done'] == 'true' and self.properties['is_public'] == 'true'
    
    def __repr__(self):
        return self.properties.__repr__()
        
    @classmethod
    def extractRuns(cls, xmlStr):
        runs = []
        def startElement(name, attrs):
            if name == "Run":
                runs.append(SRARun(attrs))
        
        xmlStr = xmlStr.replace('&lt;', '<')
        xmlStr = xmlStr.replace('&gt;', '>')
        xmlStr = xmlStr.replace('&amp;', '&')
        p = xml.parsers.expat.ParserCreate('utf-8')
        p.StartElementHandler = startElement
        p.Parse(("<Runs>" + xmlStr + "</Runs>").encode('utf-8'))
        return runs
    
def SRARunList(term):
    """
        tries to resolve a given search term to gi's
        then tries to get SRR's for those gi's
    """
    conn = httplib.HTTPConnection("eutils.ncbi.nlm.nih.gov")
    conn.request('GET', _query('esearch', dict(retmode='json', term=term, db='sra')))
    resp = conn.getresponse()
    try:
        idlist = json.load(resp)['esearchresult']['idlist']
    except:
        return []
    conn.request('GET', _query('esummary', dict(db='sra', retmode='json', id=','.join(idlist))))
    resp = conn.getresponse()
    try:
        rslt = json.load(resp)['result']
        runlist = []
        for uid in rslt['uids']:
            runlist.extend(SRARun.extractRuns(rslt[uid]['runs']))
        return runlist
    except KeyError:
        return []


def FASTA(acc):
    """
        tries to resolve the given accession into gi's
        then tries to fetch the FASTA of the first gi
        returns the defline and a Generator on the lines of FASTA
    """
    conn = httplib.HTTPConnection("eutils.ncbi.nlm.nih.gov")
    conn.request('GET', _query('esearch', dict(retmode='json', term=acc)))
    resp = conn.getresponse()
    try:
        idlist = json.load(resp)['esearchresult']['idlist']
    except:
        idlist = []
    if len(idlist) == 0:
        sys.stderr.write("Nothing was found for '{}'\n".format(acc))
        return None, None
    if len(idlist) > 1:
        sys.stderr.write("More than one ID was found for '{}'\n".format(acc))
    conn.request('GET', _query('efetch', dict(retmode='text', rettype='fasta', id=idlist[0])))
    resp = conn.getresponse()
    gen = _lines(resp)
    defline = gen.next()
    if defline.startswith('>'):
        return defline, gen
    sys.stderr.write("Unexpected output from eutils:\n{}\n".format(defline))
    return None, None


if __name__ == '__main__':
    import sys
    try:
        verb = sys.argv[1]
        query = sys.argv[2:]
    except:
        pass
    if verb == "fasta":
        for arg in query:
            defline, lines = FASTA(arg)
            if defline:
                print(defline)
                for line in lines:
                    print(line)
    if verb == "runs":
        for arg in query:
            for run in SRARunList(arg):
                print(run.accession())
