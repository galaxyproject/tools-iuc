#!/usr/bin/env python
# ======================================================================
# Script derived from the EB-eye (REST) Python client available at
# http://www.ebi.ac.uk/Tools/webservices/services/eb-eye_rest
# and distributed under the Apache License
# ======================================================================
# Load libraries

from __future__ import print_function

import os
import platform
import re
from gzip import GzipFile
from io import BytesIO
from optparse import OptionParser

from six import text_type
from six.moves.urllib.error import HTTPError
from six.moves.urllib.parse import quote
from six.moves.urllib.request import Request, urlopen
from xmltramp2 import xmltramp

# Service base URL
baseUrl = 'https://www.ebi.ac.uk/ebisearch/ws/rest'

# Debug level
debugLevel = 0


# Debug print
def printDebugMessage(functionName, message, level):
    if(level <= debugLevel):
        print('[' + functionName + '] ' + message)


# User-agent for request.
def getUserAgent():
    printDebugMessage('getUserAgent', 'Begin', 11)
    urllib_agent = 'Python-urllib/%s' % platform.python_version()
    clientRevision = '$Revision: 2468 $'
    clientVersion = '0'
    if len(clientRevision) > 11:
        clientVersion = clientRevision[11:-2]
    user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        clientVersion, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        urllib_agent
    )
    printDebugMessage('getUserAgent', 'user_agent: %s' % user_agent, 12)
    printDebugMessage('getUserAgent', 'End', 11)
    return user_agent


# Wrapper for a REST (HTTP GET) request
def restRequest(url):
    printDebugMessage('restRequest', 'Begin', 11)
    printDebugMessage('restRequest', 'url: %s' % url, 11)
    url = quote(url, safe="%/:=&?~#+!$,;'@()*[]")

    try:
        user_agent = getUserAgent()
        http_headers = {
            'User-Agent': user_agent,
            'Accept-Encoding': 'gzip'
        }
        req = Request(url, None, http_headers)
        resp = urlopen(req)
        encoding = resp.info().get('Content-Encoding')
        result = None
        if encoding is None or encoding == 'identity':
            result = text_type(resp.read(), 'utf-8')
        elif encoding == 'gzip':
            result = resp.read()
            printDebugMessage('restRequest', 'result: %s' % result, 21)
            gz = GzipFile(
                fileobj=BytesIO(result),
                mode="r")
            result = text_type(gz.read(), 'utf-8')
        else:
            raise Exception('Unsupported Content-Encoding')
        resp.close()
    except HTTPError as ex:
        raise ex
    printDebugMessage('restRequest', 'result: %s' % result, 11)
    printDebugMessage('restRequest', 'End', 11)
    return result


def hasSubdomains(domainInfo):
    for dir in domainInfo._dir:
        if dir._name == 'subdomains':
            return True
    return False


def extractUsefulFields(fieldInfos):
    searchable = []
    retrievable = []

    for fieldInfo in fieldInfos:
        if fieldInfo('id') == "$facets":
            continue

        options = fieldInfo['options']['option':]
        for option in options:
            if option("name") == "searchable" and str(option) == "true":
                searchable.append(fieldInfo('id'))
            if option("name") == "retrievable" and str(option) == "true":
                retrievable.append(fieldInfo('id'))
    return searchable, retrievable


def extractLowerLevelDomains(domainInfo, domains):
    if hasSubdomains(domainInfo):
        subdomains = domainInfo['subdomains']['domain':]
        for subdomain in subdomains:
            domains = extractLowerLevelDomains( subdomain, domains)
    else:
        searchable, retrievable = extractUsefulFields(
            domainInfo['fieldInfos']['fieldInfo':])

        domain_id = domainInfo('id')
        domains.setdefault(domain_id, {})
        domains[domain_id]["name"] = domainInfo('name')
        domains[domain_id]["searchable_fields"] = sorted(searchable)
        domains[domain_id]["retrievable_fields"] = sorted(retrievable)
    return domains


# Get domain Hierarchy
def getDomainHierarchy():
    requestUrl = baseUrl + '/allebi'
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    allebi = doc['domains']['domain']
    lower_level_domains = extractLowerLevelDomains(allebi, {})
    printDebugMessage('getDomainHierarchy', 'End', 1)
    return lower_level_domains


# Check if a databaseInfo matches a database name.
def is_database(dbInfo, dbName):
    printDebugMessage('is_database', 'Begin', 11)
    retVal = False
    if str(dbInfo.name) == dbName:
        retVal = True
    else:
        for dbAlias in dbInfo.aliasList:
            if str(dbAlias) == dbName:
                retVal = True
    printDebugMessage('is_database', 'retVal: %s' % retVal, 11)
    printDebugMessage('is_database', 'End', 11)
    return retVal


# Get number of results
def getNumberOfResults(domain, query):
    printDebugMessage('getNumberOfResults', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '?query=' + query + '&size=0'
    printDebugMessage('getNumberOfResults', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    numberOfResults = int(str(doc['hitCount']))
    printDebugMessage('getNumberOfResults', 'End', 1)
    return numberOfResults


def makeRequest(requestUrl):
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    entries = doc['entries']['entry':]
    formatted_output = printEntries(entries)
    return formatted_output


# Get search results
def getResults(domain, query, fields):
    numberOfResults = getNumberOfResults(domain, query)
    maximum_size = 100
    quotient = numberOfResults // maximum_size
    start = 0

    printDebugMessage('getResults', 'Begin', 1)
    request_output = "%s\tlink\n" % (fields.replace(",", "\t"))
    for i in range(quotient):
        start = maximum_size * i
        requestUrl = baseUrl + '/' + domain + '?query=' + query
        requestUrl += '&fields=' + fields + '&size=' + str(maximum_size)
        requestUrl += '&start=' + str(start) + '&fieldurl=true'
        request_output += makeRequest(requestUrl)

    if (numberOfResults % 100) > 0:
        start = maximum_size * quotient
        remainder = numberOfResults - start
        requestUrl = baseUrl + '/' + domain + '?query=' + query
        requestUrl += '&fields=' + fields + '&size=' + str(remainder)
        requestUrl += '&start=' + str(start) + '&fieldurl=true'
        request_output += makeRequest(requestUrl)

    print(request_output)


def printEntries(entries):
    output = ""
    printDebugMessage('printEntries', 'Begin', 1)
    for entry in entries:
        sep = ""
        for field in entry['fields']['field':]:
            output += "%s" % (sep)
            fields = field['values']['value':]
            if len(fields) > 0:
                sub_sep = ""
                for value in field['values']['value':]:
                    output += "%s%s" % (sub_sep, value)
                    sub_sep = ","
            sep = "\t"

        if hasFieldUrls(entry):
            output += "%s" % (sep)
            sub_sep = ""
            for fieldurl in entry['fieldURLs']['fieldURL':]:
                output += "%s%s" % (sub_sep, str(fieldurl))
                sub_sep = ","
            sep = "\t"
        if hasViewUrls(entry):
            output += "%s" % (sep)
            sub_sep = ""
            for viewurl in entry['viewURLs']['viewURL':]:
                output += "%s%s" % (sub_sep, str(viewurl))
                sub_sep = ","
        output += "\n"
    printDebugMessage('printEntries', 'End', 1)
    return output


def hasFieldUrls(entry):
    for dir in entry._dir:
        if dir._name == 'fieldURLs':
            return True
    return False


def hasViewUrls(entry):
    for dir in entry._dir:
        if dir._name == 'viewURLs':
            return True
    return False


def getRunLink(run_id):
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
    p = re.compile('http')
    fieldURL = p.sub('https', fieldURL)
    print(fieldURL)


if __name__ == '__main__':
    # Usage message
    usage = """
      %prog getDomainHierarchy
      %prog getResults <domain> <query> <fields>
      %prog getRunLink <runId>
      """

    description = "Tools to query and download data from several EMBL-EBI databases"
    description += "The searching tools are using the EB-eye search engine. "
    description += "http://www.ebi.ac.uk/ebisearch/"
    # Process command-line options
    parser = OptionParser(
        usage=usage,
        description=description,
        version='1.0')
    (options, args) = parser.parse_args()

    # No arguments, print usage
    if len(args) < 1:
        parser.print_help()

    # Get domain hierarchy
    elif args[0] == 'getDomainHierarchy':
        getDomainHierarchy()

    # Get search results
    elif args[0] == 'getResults':
        if len(args) < 4:
            print('domain, query and fields should be given.')
        else:
            getResults(args[1], args[2], args[3])

    # Get run link results
    elif args[0] == 'getRunLink':
        if len(args) < 2:
            print('run id should be given.')
        else:
            getRunLink(args[1])

    # Unknown argument combination, display usage
    else:
        print('Error: unrecognised argument combination')
        parser.print_help()
