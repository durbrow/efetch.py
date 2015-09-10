#!python
"""
    Retrieve sequence data as FASTA from NCBI using eutils
    
    The following example will fetch the sequence of chromosome 8 of GRCh37:
        python efetch.py CM000670.1

    See:
        http://www.ncbi.nlm.nih.gov/books/NBK25499/
        http://www.ncbi.nlm.nih.gov/books/NBK25500/
        http://www.ncbi.nlm.nih.gov/books/NBK25501/
"""

import httplib
import json
import sys

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
    for arg in sys.argv[1:]:
        defline, lines = FASTA(arg)
        print(defline)
        for line in lines:
            print(line)
