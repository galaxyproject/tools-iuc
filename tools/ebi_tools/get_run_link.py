#!/usr/bin/env python
# ======================================================================
# Script derived from the EB-eye (REST) Python-3 client available at
# http://www.ebi.ac.uk/Tools/webservices/services/eb-eye_rest
# and distributed under the Apache License
# ======================================================================
# Load libraries
import platform
import os
import urllib
import argparse
from gzip import GzipFile
from xmltramp2 import xmltramp
import re
# python2
from StringIO import StringIO
import urllib2
# python3
# import urllib.request as urllib2


# Service base URL
baseUrl = 'http://www.ebi.ac.uk/ebisearch/ws/rest'

# Debug level
debugLevel = 0


# Debug print
def printDebugMessage(functionName, message, level):
    if(level <= debugLevel):
        print ('[' + functionName + '] ' + message)


# User-agent for request.
def getUserAgent():
    printDebugMessage('getUserAgent', 'Begin', 11)
    urllib_agent = 'Python-urllib/%s' % urllib2.__version__
    clientRevision = '$Revision: 2468 $'
    clientVersion = '0'
    if len(clientRevision) > 11:
        clientVersion = clientRevision[11:-2]
    user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        clientVersion, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        urllib_agent
    )
    printDebugMessage('getUserAgent', 'user_agent: ' + user_agent, 12)
    printDebugMessage('getUserAgent', 'End', 11)
    return user_agent


# Wrapper for a REST (HTTP GET) request
def restRequest(url):
    printDebugMessage('restRequest', 'Begin', 11)
    printDebugMessage('restRequest', 'url: ' + url, 11)
    # python 2
    url = urllib.quote(url, safe="%/:=&?~#+!$,;'@()*[]")
    # python 3
    # url = urllib.request.quote(url, safe="%/:=&?~#+!$,;'@()*[]")

    try:
        user_agent = getUserAgent()
        http_headers = {
            'User-Agent': user_agent,
            'Accept-Encoding': 'gzip'
        }
        req = urllib2.Request(url, None, http_headers)
        resp = urllib2.urlopen(req)
        # python2
        encoding = resp.info().getheader('Content-Encoding')
        # python3
        # encoding = resp.info().__getitem__('Content-Encoding')
        result = None
        if encoding is None or encoding == 'identity':
            # python2
            result = resp.read()
            # python3
            # result = str(resp.read(), 'utf-8')
        elif encoding == 'gzip':
            result = resp.read()
            printDebugMessage('restRequest', 'result: ' + str(result), 21)
            # python2
            gz = GzipFile(
                fileobj=StringIO(result),
                mode="r")
            result = gz.read()
            # python3
            # result = str(gzip.decompress(result), 'utf-8')
        else:
            raise Exception('Unsupported Content-Encoding')
        resp.close()
    except urllib2.HTTPError as ex:
        raise ex
    printDebugMessage('restRequest', 'result: ' + result, 11)
    printDebugMessage('restRequest', 'End', 11)
    return result


# Get run link
def get_run_link(run_id):
    printDebugMessage('getEntries', 'Begin', 1)
    requestUrl = baseUrl + '/metagenomics_runs/entry/' + run_id + '?fieldurl=true'
    printDebugMessage('getEntries', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    entries = doc['entries']['entry':]
    fieldURL = ''
    for entry in entries:
        for fieldurl in entry['fieldURLs']['fieldURL':]:
            fieldURL += str(fieldurl)
    printDebugMessage('getEntries', 'End', 1)
    return fieldURL


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_id', required=True)
    args = parser.parse_args()

    url = get_run_link(args.run_id)
    p = re.compile('http')
    url = p.sub('https', url)
    print url
