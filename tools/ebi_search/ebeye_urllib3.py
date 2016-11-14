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
from optparse import OptionParser
from gzip import GzipFile
from xmltramp2 import xmltramp
# python2
from StringIO import StringIO
import urllib2
# python3
# import urllib.request as urllib2

# Output level
outputLevel = 1
# Debug level
debugLevel = 0

# Service base URL
baseUrl = 'http://www.ebi.ac.uk/ebisearch/ws/rest'


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


def hasSubdomains(domainInfo):
    for dir in domainInfo._dir:
        if dir._name == 'subdomains':
            return True
    return False


def printDomains(domainInfo, indent):
    printDebugMessage('printDomains', 'Begin', 1)
    print (indent + domainInfo('id') + ':' + domainInfo('name'))
    if hasSubdomains(domainInfo):
        subdomains = domainInfo['subdomains']['domain':]
        for subdomain in subdomains:
            printDomains(subdomain, indent + '\t')
    printDebugMessage('printDomains', 'End', 1)


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
    printDebugMessage('getDomainHierarchy', 'Begin', 1)
    requestUrl = baseUrl + '/allebi'
    printDebugMessage('getDomainHierarchy', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    allebi = doc['domains']['domain']
    # printDomains(allebi, '')
    lower_level_domains = extractLowerLevelDomains(allebi, {})
    # domainInfoList = doc['result.domains.domain':]
    # for domainInfo in domainInfoList:
    #    print (domainInfo.id)
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
    printDebugMessage('is_database', 'retVal: ' + str(retVal), 11)
    printDebugMessage('is_database', 'End', 11)
    return retVal


# Get domain details
def getDomainDetails(domain):
    printDebugMessage('getDomainDetails', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain
    printDebugMessage('getDomainDetails', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    domainInfo = doc['domains']['domain']
    printDomainDetails(domainInfo)
    printDebugMessage('getDomainDetails', 'End', 1)


def printDomainDetails(domainInfo):
    printDebugMessage('printDomainDetails', 'Begin', 1)
    print (domainInfo('name') + ' (' + domainInfo('id') + ')')
    if hasSubdomains(domainInfo):
        subdomains = domainInfo['subdomains']['domain':]
        for subdomain in subdomains:
            printDomainDetails(subdomain)
    else:
        indexInfos = domainInfo['indexInfos']['indexInfo':]
        for indexInfo in indexInfos:
            print (indexInfo('name') + ': ' + str(indexInfo))
        print
        fieldInfos = domainInfo['fieldInfos']['fieldInfo':]
        print ('field_id\tsearchable\tretrievable\tsortable\tfacet')
        fieldStr = ''
        for fieldInfo in fieldInfos:
            fieldStr = fieldInfo('id') + '\t'
            options = fieldInfo['options']['option':]
            for option in options:
                fieldStr += str(option) + '\t'
            print (fieldStr)
        print
    printDebugMessage('printDomainDetails', 'End', 1)


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


def printNumber(num):
    printDebugMessage('printNumber', 'Begin', 1)
    print(num)
    printDebugMessage('printNumber', 'End', 1)


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
    quotient = numberOfResults / maximum_size
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


def printFacets(facets):
    printDebugMessage('printFacets', 'Begin', 1)
    for facet in facets:
        print (facet('label') + ': ' + facet('id'))
        for facetValue in facet['facetValues']['facetValue':]:
            printFacetValue(facetValue, 0)
        print
    printDebugMessage('printFacets', 'End', 1)


def printFacetValue(facetValue, depth=0):
    printDebugMessage('printFacetValue', 'Begin', 1)

    print ("%s%s (%s):%s" % (
        '\t' * depth,
        str(facetValue['label']),
        str(facetValue['value']),
        str(facetValue['count'])))

    if hasFacetValueChildren(facetValue):
        for child in facetValue['children']['facetValue':]:
            printFacetValue(child, depth + 1)

    printDebugMessage('printFacetValue', 'End', 1)


def hasFacetValueChildren(facetValue):
    for dir in facetValue._dir:
        if dir._name == 'children':
            return True
    return False


# Get search results with facets
def getFacetedResults(
        domain, query, fields, size='', start='', fieldurl='', viewurl='',
        sortfield='', order='', sort='', facetcount='10', facetfields='',
        facets='', facetsdepth=''):
    printDebugMessage('getFacetedResults', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '?query=' + query
    requestUrl += '&fields=' + fields + '&size=' + size + '&start=' + start
    requestUrl += '&fieldurl=' + fieldurl + '&viewurl=' + viewurl
    requestUrl += '&sortfield=' + sortfield + '&order=' + order
    requestUrl += '&sort=' + sort + '&facetcount=' + facetcount
    requestUrl += '&facetfields=' + facetfields + '&facets=' + facets
    requestUrl += '&facetsdepth=' + facetsdepth
    printDebugMessage('getFacetedResults', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    entries = doc['entries']['entry':]
    printEntries(entries)
    print
    facets = doc['facets']['facet':]
    printFacets(facets)
    printDebugMessage('getFacetedResults', 'End', 1)


# Get entry details
def getEntries(domain, entryids, fields, fieldurl='', viewurl=''):
    printDebugMessage('getEntries', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/entry/' + entryids
    requestUrl += '?fields=' + fields + '&fieldurl=' + fieldurl
    requestUrl += '&viewurl=' + viewurl
    print(requestUrl)
    printDebugMessage('getEntries', requestUrl, 2)
    request_output = makeRequest(requestUrl)
    print(request_output)
    printDebugMessage('getEntries', 'End', 1)


# Get domain ids referenced in a domain
def getDomainsReferencedInDomain(domain):
    printDebugMessage('getDomainsReferencedInDomain', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/xref/'
    printDebugMessage('getDomainsReferencedInDomain', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    for domain in doc['domains']['domain':]:
        print (domain('id'))
    printDebugMessage('getDomainsReferencedInDomain', 'End', 1)


# Get domain ids referenced in an entry
def getDomainsReferencedInEntry(domain, entryid):
    printDebugMessage('getDomainsReferencedInEntry', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/entry/' + entryid + '/xref/'
    printDebugMessage('getDomainsReferencedInEntry', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    for domain in doc['domains']['domain':]:
        print (domain('id'))
    printDebugMessage('getDomainsReferencedInEntry', 'End', 1)


# Get cross-references
def getReferencedEntries(
        domain, entryids, referenceddomain, fields, size='', start='',
        fieldurl='', viewurl='', facetcount='', facetfields='', facets=''):
    printDebugMessage('getReferencedEntries', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/entry/' + entryids
    requestUrl += '/xref/' + referenceddomain + '?fields=' + fields
    requestUrl += '&size=' + size + '&start=' + start + '&fieldurl=' + fieldurl
    requestUrl += '&viewurl=' + viewurl + '&facetcount=' + facetcount
    requestUrl += '&facetfields=' + facetfields + '&facets=' + facets
    printDebugMessage('getReferencedEntries', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    entries = doc['entries']['entry':]
    for entry in entries:
        printEntries(entry['references']['reference':])
        if hasReferenceFacet(entry):
            printFacets(entry['referenceFacets']['referenceFacet':])
        print
    printDebugMessage('getEntries', 'End', 1)


def hasReferenceFacet(entry):
    for dir in entry._dir:
        if dir._name == 'referenceFacets':
            return True
    return False


# Get top terms
def getTopTerms(domain, field, size='', excludes='', excludesets=''):
    printDebugMessage('getTopTerms', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/topterms/' + field
    requestUrl += '?size=' + size + '&excludes=' + excludes
    requestUrl += '&excludesets=' + excludesets
    printDebugMessage('getTopTerms', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    topterms = doc['topTerms']['term':]
    for term in topterms:
        printTerm(term)
        print
    printDebugMessage('getTopTerms', 'End', 1)


def printTerm(term):
    printDebugMessage('printTerm', 'Begin', 1)
    print (str(term['text']) + ': ' + str(term['docFreq']))
    print
    printDebugMessage('printTerm', 'End', 1)


# Get similar documents to a given one
def getMoreLikeThis(
        domain, entryid, targetDomain, fields, size='', start='', fieldurl='',
        viewurl='', mltfields='', mintermfreq='', mindocfreq='',
        maxqueryterm='', excludes='', excludesets=''):
    printDebugMessage('getMoreLikeThis', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/entry/' + entryid
    requestUrl += '/morelikethis/' + targetDomain + '?size=' + size
    requestUrl += '&start=' + start + '&fields=' + fields
    requestUrl += '&fieldurl=' + fieldurl + '&viewurl=' + viewurl
    requestUrl += '&mltfields=' + mltfields + '&mintermfreq=' + mintermfreq
    requestUrl += '&mindocfreq=' + maxqueryterm + '&maxqueryterm=' + mindocfreq
    requestUrl += '&excludes=' + excludes + '&excludesets=' + excludesets
    printDebugMessage('getMoreLikeThis', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    entries = doc['entries']['entry':]
    printEntries(entries)
    printDebugMessage('getMoreLikeThis', 'End', 1)


# Get suggestions
def getAutoComplete(domain, term):
    printDebugMessage('getAutoComplete', 'Begin', 1)
    requestUrl = baseUrl + '/' + domain + '/autocomplete?term=' + term
    printDebugMessage('getAutoComplete', requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    suggestions = doc['suggestions']['suggestion':]
    printSuggestions(suggestions)
    printDebugMessage('getAutoComplete', 'End', 1)


def printSuggestions(suggestions):
    printDebugMessage('printSuggestions', 'Begin', 1)
    for suggetion in suggestions:
        print (str(suggetion['suggestion']))
    print
    printDebugMessage('printSuggestions', 'End', 1)


if __name__ == '__main__':
    # Usage message
    usage = """
      %prog getDomainHierarchy
      %prog getDomainDetails  <domain>

      %prog getNumberOfResults <domain> <query>
      %prog getResults <domain> <query> <fields>
      [OPTIONS: --size | --start | --fieldurl | --viewurl | --sortfield |
      --order | --sort ]
      %prog getFacetedResults <domain> <query> <fields>
      [OPTIONS: --size | --start | --fieldurl | --viewurl | --sortfield |
      --order | --sort | --facetcount | --facetfields | --facets |
      --facetsdepth ]

      %prog getEntries <domain> <entryids> <fields>
      [OPTIONS: --fieldurl | --viewurl]

      %prog getDomainsReferencedInDomain <domain>
      %prog getDomainsReferencedInEntry  <domain> <entryid>
      %prog getReferencedEntries <domain> <entryids> <referencedDomain>
      <fields>
      [OPTIONS: --size | --start | --fieldurl | --viewurl | --facetcount |
      --facetfields | --facets]

      %prog getTopTerms <domain> <field>
      [OPTIONS: --size | --excludes | --excludesets]

      %prog getAutoComplete   <domain> <term>

      %prog getMoreLikeThis   <domain> <entryid> <fields>
      [OPTIONS: --size | --start | --fieldurl | --viewurl | --mltfields |
      --mintermfreq | --mindocfreq | --maxqueryterm | --excludes |
      --excludesets]
      %prog getExtendedMoreLikeThis <domain> <entryid> <targetDomain> <fields>
      [OPTIONS: --size | --start | --fieldurl | --viewurl | --mltfields |
      --mintermfreq | --mindocfreq | --maxqueryterm | --excludes |
      --excludesets]"""

    description = "Search at EMBL-EBI using the EB-eye search engine. "
    description += "For more information on EB-eye refer to"
    description += " http://www.ebi.ac.uk/ebisearch/"
    version = "$Id: dbfetch_urllib2.py 2468 2013-01-25 14:01:01Z hpm $"
    # Process command-line options
    parser = OptionParser(
        usage=usage,
        description=description,
        version=version)
    parser.add_option('--size', help='number of entries to retrieve')
    parser.add_option('--start', help='index of the first entry in results')
    parser.add_option('--fieldurl', help='whether field links are included')
    parser.add_option('--viewurl', help='whether view links are included')
    parser.add_option('--sortfield', help='field id to sort')
    parser.add_option('--order', help='sort in ascending/descending order')
    parser.add_option('--sort', help='sort criteria')
    parser.add_option(
        '--facetcount',
        help='number of facet values to retrieve')
    parser.add_option(
        '--facetfields',
        help='field ids associated with facets to retrieve')
    parser.add_option('--facets', help='list of selected facets')
    parser.add_option('--facetsdepth', help='depth in hierarchical facet')
    parser.add_option(
        '--mltfields',
        help='field ids  to be used for generating a morelikethis query')
    parser.add_option(
        '--mintermfreq',
        help='frequency below which terms will be ignored in the base '
        'document')
    parser.add_option(
        '--mindocfreq',
        help='frequency at which words will be ignored which do not occur in '
        'at least this many documents')
    parser.add_option(
        '--maxqueryterm',
        help='maximum number of query terms that will be included in any'
        ' generated query')
    parser.add_option('--excludes', help='terms to be excluded')
    parser.add_option('--excludesets', help='stop word sets to be excluded')

    parser.add_option(
        '--quiet',
        action='store_true',
        help='decrease output level')
    parser.add_option(
        '--verbose',
        action='store_true',
        help='increase output level')
    parser.add_option(
        '--baseUrl',
        default=baseUrl,
        help='base URL for EBI Search')
    parser.add_option(
        '--debugLevel',
        type='int',
        default=debugLevel,
        help='debug output level')
    (options, args) = parser.parse_args()

    # Increase output level.
    if options.verbose:
        outputLevel += 1

    # Decrease output level.
    if options.quiet:
        outputLevel -= 1

    # Debug level.
    if options.debugLevel:
        debugLevel = options.debugLevel

    # Base URL for service.
    if options.baseUrl:
        baseUrl = options.baseUrl

    # No arguments, print usage
    if len(args) < 1:
        parser.print_help()
    # Get domain hierarchy
    elif args[0] == 'getDomainHierarchy':
        getDomainHierarchy()
    # Get domain details
    elif args[0] == 'getDomainDetails':
        if len(args) < 2:
            print ('domain should be given.')
        else:
            getDomainDetails(args[1])

    # Get number of results
    elif args[0] == 'getNumberOfResults':
        if len(args) != 3:
            print ('domain and query should be given.')
        else:
            getNumberOfResults(args[1], args[2])

    # Get search results
    elif args[0] == 'getResults':
        if len(args) < 4:
            print ('domain, query and fields should be given.')
        else:
            size = options.size if options.size else ''
            start = options.start if options.start else ''
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            sortfield = options.sortfield if options.sortfield else ''
            order = options.order if options.order else ''
            sort = options.sort if options.sort else ''
            getResults(args[1], args[2], args[3])
    # Get search results with facets
    elif args[0] == 'getFacetedResults':
        if len(args) < 4:
            print ('domain, query and fields should be given.')
        else:
            size = options.size if options.size else ''
            start = options.start if options.start else ''
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            sortfield = options.sortfield if options.sortfield else ''
            order = options.order if options.order else ''
            sort = options.sort if options.sort else ''
            facetcount = options.facetcount if options.facetcount else ''
            facetfields = options.facetfields if options.facetfields else ''
            facets = options.facets if options.facets else ''
            facetsdepth = options.facetsdepth if options.facetsdepth else ''
            getFacetedResults(
                args[1],
                args[2],
                args[3],
                size,
                start,
                fieldurl,
                viewurl,
                sortfield,
                order,
                sort,
                facetcount,
                facetfields,
                facets,
                facetsdepth)
    # Get entry details.
    elif args[0] == 'getEntries':
        if len(args) < 4:
            print ('domain, entry ids and fields should be given.')
        else:
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            getEntries(args[1], args[2], args[3], fieldurl, viewurl)
    # Get domain ids referenced in a domain
    elif args[0] == 'getDomainsReferencedInDomain':
        if len(args) < 2:
            print ('domain sholud be given.')
        else:
            getDomainsReferencedInDomain(args[1])
    # Get domain ids referenced in an entry
    elif args[0] == 'getDomainsReferencedInEntry':
        if len(args) < 3:
            print ('domain and entry id should be given.')
        else:
            getDomainsReferencedInEntry(args[1], args[2])
    # Get cross-references
    elif args[0] == 'getReferencedEntries':
        if len(args) < 5:
            print (
                'domain, entryids, referencedDomain and fields should be '
                'given.')
        else:
            size = options.size if options.size else ''
            start = options.start if options.start else ''
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            facetcount = options.facetcount if options.facetcount else ''
            facetfields = options.facetfields if options.facetfields else ''
            facets = options.facets if options.facets else ''
            getReferencedEntries(
                args[1],
                args[2],
                args[3],
                args[4],
                size,
                start,
                fieldurl,
                viewurl,
                facetcount,
                facetfields,
                facets)
    # Get top terms
    elif args[0] == 'getTopTerms':
        if len(args) < 3:
            print ('domain and field should be given.')
        else:
            size = options.size if options.size else ''
            excludes = options.excludes if options.excludes else ''
            excludesets = options.excludesets if options.excludesets else ''
            getTopTerms(args[1], args[2], size, excludes, excludesets)
    # Get similar documents to a given one
    elif args[0] == 'getMoreLikeThis':
        if len(args) < 4:
            print ('domain, entryid and fields should be given.')
        else:
            size = options.size if options.size else ''
            start = options.start if options.start else ''
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            mltfields = options.mltfields if options.mltfields else ''
            mintermfreq = options.mintermfreq if options.mintermfreq else ''
            mindocfreq = options.mindocfreq if options.mindocfreq else ''
            maxqueryterm = options.maxqueryterm if options.maxqueryterm else ''
            excludes = options.excludes if options.excludes else ''
            excludesets = options.excludesets if options.excludesets else ''
            getMoreLikeThis(
                args[1],
                args[2],
                args[1],
                args[3],
                size,
                start,
                fieldurl,
                viewurl,
                mltfields,
                mintermfreq,
                mindocfreq,
                maxqueryterm,
                excludes,
                excludesets)
    # Get similar documents to a given one
    elif args[0] == 'getExtendedMoreLikeThis':
        if len(args) < 5:
            print ('domain, entryid, targetDomain and fields should be given.')
        else:
            size = options.size if options.size else ''
            start = options.start if options.start else ''
            fieldurl = options.fieldurl if options.fieldurl else ''
            viewurl = options.viewurl if options.viewurl else ''
            mltfields = options.mltfields if options.mltfields else ''
            mintermfreq = options.mintermfreq if options.mintermfreq else ''
            mindocfreq = options.mindocfreq if options.mindocfreq else ''
            maxqueryterm = options.maxqueryterm if options.maxqueryterm else ''
            excludes = options.excludes if options.excludes else ''
            excludesets = options.excludesets if options.excludesets else ''
            getMoreLikeThis(
                args[1],
                args[2],
                args[3],
                args[4],
                size,
                start,
                fieldurl,
                viewurl,
                mltfields,
                mintermfreq,
                mindocfreq,
                maxqueryterm,
                excludes,
                excludesets)

    elif args[0] == 'getAutoComplete':
        if len(args) < 3:
            print ('domain and term should be given.')
        else:
            getAutoComplete(args[1], args[2])
    # Unknown argument combination, display usage
    else:
        print ('Error: unrecognised argument combination')
        parser.print_help()
