import json
import os
from io import StringIO

from Bio import Entrez

Entrez.tool = "GalaxyEutils_1_0"
BATCH_SIZE = 200


class Client(object):

    def __init__(self, history_file=None, user_email=None, admin_email=None):
        self.using_history = False
        self.using_parsedids = False

        if user_email is not None and admin_email is not None:
            Entrez.email = ';'.join((admin_email, user_email))
        elif user_email is not None:
            Entrez.email = user_email
        elif admin_email is not None:
            Entrez.email = admin_email
        else:
            Entrez.email = os.environ.get('NCBI_EUTILS_CONTACT', None)

        if Entrez.email is None:
            raise Exception("Cannot continue without an email; please set "
                            "administrator email in NCBI_EUTILS_CONTACT")

        if history_file is not None:
            with open(history_file, 'r') as handle:
                data = json.loads(handle.read())
                # esearch
                if 'QueryKey' in data:
                    self.query_key = data['QueryKey']
                    self.webenv = data['WebEnv']
                    self.query_keys = []
                    self.query_keys += [data['QueryKey']]
                    self.using_history = True
                elif 'query_key' in data:
                    self.query_key = data['query_key']
                    self.webenv = data['WebEnv']
                    self.query_keys = []
                    self.query_keys += [data['query_key']]
                    self.using_history = True
                elif 'esearchresult' in data:
                    self.query_key = data['esearchresult']['querykey']
                    self.webenv = data['esearchresult']['webenv']
                    self.query_keys = []
                    self.query_keys += [data['esearchresult']['querykey']]
                    self.using_history = True
                # elink
                elif 'linksets' in data:
                    # elink for cmd=neighbor_history
                    if 'linksetdbhistories' in data['linksets'][0]:
                        self.webenv = data['linksets'][0]['webenv']
                        self.query_key = data['linksets'][0]['linksetdbhistories'][0]['querykey']
                        self.using_history = True
                    # elink for cmd=neighbor|neighbor_score
                    elif 'linksetdbs' in data['linksets'][0]:
                        self.using_parsedids = True
                        # elink for neighbor
                        if isinstance(data['linksets'][0]['linksetdbs'][0]['links'][0], str):
                            self.idstr = ','.join(data['linksets'][0]['linksetdbs'][0]['links'])
                        # elink for neighbor_score
                        else:
                            self.idstr = ','.join(map(lambda x: x['id'], data['linksets'][0]['linksetdbs'][0]['links']))
                    if 'linksetdbhistories' in data['linksets'][0]:
                        self.webenv = data['linksets'][0]['webenv']
                        self.query_keys = []
                        for query in data['linksets'][0]['linksetdbhistories']:
                            if 'querykey' in query:
                                self.query_keys += [query['querykey']]
                else:
                    print("No match")
                    print(data)

    def get_history(self):
        if self.using_history:
            return {
                'query_key': self.query_key,
                'WebEnv': self.webenv,
            }
        elif self.using_parsedids:
            return {
                'id': self.idstr,
            }
        else:
            return {}

    def get_histories(self):
        histories = []
        for key in self.query_keys:
            histories += [{'WebEnv': self.webenv, 'query_key': key}]
        return histories

    def post(self, database, **payload):
        return json.dumps(Entrez.read(Entrez.epost(database, **payload)), indent=4)

    def fetch(self, db, ftype=None, **payload):
        os.makedirs("downloads")

        if 'id' in payload:
            summary = self.id_summary(db, payload['id'])
        elif 'WebEnv' not in payload or 'query_key' not in payload:
            summary = self.history_summary(db)
        else:
            summary = payload

        count = len(summary)
        payload['retmax'] = BATCH_SIZE

        # This may be bad. I'm not sure yet. I think it will be ... but UGH.
        for i in range(0, count, BATCH_SIZE):
            payload['retstart'] = i
            file_path = os.path.join('downloads', 'EFetch Results Chunk %s.%s' % (i, ftype))
            with open(file_path, 'w') as handle:
                handle.write(Entrez.efetch(db, **payload).read())

    def id_summary(self, db, id_list):
        payload = {
            'db': db,
            'id': id_list,
        }
        return Entrez.read(Entrez.esummary(**payload))

    def history_summary(self, db):
        if not self.using_history:
            raise Exception("History must be available for this method")

        payload = {
            'db': db,
            'query_key': self.query_key,
            'WebEnv': self.webenv,
        }
        return Entrez.read(Entrez.esummary(**payload))

    def summary(self, **payload):
        return Entrez.esummary(**payload).read()

    def link(self, **payload):
        return Entrez.elink(**payload).read()

    def extract_history_from_xml_file(self, xml_file):
        history = {}
        with open(xml_file, 'r') as handle:
            xml_str = handle.read()
            history = self.extract_history_from_xml(xml_str)
        return history

    def extract_history_from_xml(self, xml_str):
        try:
            parsed_data = Entrez.read(StringIO(xml_str))
            history = {}
            gotit = 0

            # New code doesn't work for esearch input to elink - Parsing esearch output (reading an xml history) does not work as an elink input payload, which needs 'QueryKey'.  Notably, if parsing elink output as input to elink, conversion of xml 'QueryKey' to 'query_key' is needed for some reason.  Also Notably, efetch returned results using the 'QueryKey' key
            # For esearch xml history results
            if 'QueryKey' in parsed_data:
                history['query_key'] = parsed_data['QueryKey']
                gotit += 1
            if 'WebEnv' in parsed_data:
                history['WebEnv'] = parsed_data['WebEnv']
                gotit += 1
            # For elink xml history results
            if gotit < 2:
                if 'LinkSetDbHistory' in parsed_data[0]:
                    if 'QueryKey' in parsed_data[0]['LinkSetDbHistory'][0]:
                        history['query_key'] = parsed_data[0]['LinkSetDbHistory'][0]['QueryKey']
                        gotit += 1
                if 'WebEnv' in parsed_data[0]:
                    history['WebEnv'] = parsed_data[0]['WebEnv']
                    gotit += 1
                if gotit < 2:
                    raise Exception("Could not find WebEnv in xml response")
        except Exception as e:
            print("Error parsing...")
            print(xml_str)
            raise e

        return history

    def extract_histories_from_xml_file(self, xml_file):
        histories = []
        with open(xml_file, 'r') as handle:
            xml_str = handle.read()
            histories = self.extract_histories_from_xml(xml_str)
        return histories

    def extract_histories_from_xml(self, xml_str):
        try:
            parsed_data = Entrez.read(StringIO(xml_str))
            histories = []
            gotit = 0

            # New code doesn't work for esearch input to elink - Parsing esearch output (reading an xml history) does not work as an elink input payload, which needs 'QueryKey'.  Notably, if parsing elink output as input to elink, conversion of xml 'QueryKey' to 'query_key' is needed for some reason.  Also Notably, efetch returned results using the 'QueryKey' key
            # For esearch xml history results
            if 'QueryKey' in parsed_data:
                tmp_hist = {}
                tmp_hist['query_key'] = parsed_data['QueryKey']
                gotit += 1
            if 'WebEnv' in parsed_data:
                tmp_hist['WebEnv'] = parsed_data['WebEnv']
                gotit += 1
            if gotit == 2:
                histories += [tmp_hist]
            # For elink xml history results
            else:
                gotenv = 0
                if 'LinkSetDbHistory' in parsed_data[0]:
                    for query in parsed_data[0]['LinkSetDbHistory']:
                        tmp_hist = {}
                        if 'WebEnv' in parsed_data[0]:
                            tmp_hist['WebEnv'] = parsed_data[0]['WebEnv']
                        if 'QueryKey' in query:
                            tmp_hist['query_key'] = query['QueryKey']
                            histories += [tmp_hist]
                            gotit += 1
                if gotit == 0 and gotenv == 0:
                    raise Exception("Could not find WebEnv in xml response")
        except Exception as e:
            print("Error parsing...")
            print(xml_str)
            raise e

        return histories

    def search(self, **payload):
        return Entrez.esearch(**payload).read()

    def info(self, **kwargs):
        return Entrez.einfo(**kwargs).read()

    def gquery(self, **kwargs):
        return Entrez.egquery(**kwargs).read()

    def citmatch(self, **kwargs):
        return Entrez.ecitmatch(**kwargs).read()

    @classmethod
    def jsonstring2jsondata(cls, json_str):
        json_handle = StringIO(json_str)
        json_data = json.loads(json_handle.read())
        return json_data

    @classmethod
    def jsonfile2UIlist(cls, json_file):
        merged_ids = []
        with open(json_file, 'r') as handle:
            json_data = json.loads(handle.read())
            for id in cls.jsondata2UIlist(json_data):
                merged_ids += [id]
        return merged_ids

    @classmethod
    def jsondata2UIlist(cls, json_data):
        merged_ids = []

        # Always prioritize the result links as opposed to the search links
        # elink - retrieves linked IDs for cmd=neighbor|neighbor_score only
        if 'linksets' in json_data:
            for lnk in json_data['linksets'][0]['linksetdbs']:
                if 'links' in lnk:
                    for id in lnk['links']:
                        # elink for neighbor
                        if isinstance(id, str):
                            merged_ids.append(id)
                        # elink for neighbor_score
                        else:
                            merged_ids.append(id['id'])
        # esearch
        elif 'esearchresult' in json_data:
            for id in json_data['esearchresult']['idlist']:
                merged_ids += [id]

        return merged_ids

    @classmethod
    def xmlfile2UIlist(cls, xml_file):
        merged_ids = []
        with open(xml_file, 'r') as handle:
            xml_data = Entrez.read(handle)
            for id in cls.xmldata2UIlist(xml_data):
                merged_ids += [id]
        return merged_ids

    @classmethod
    def xmlstring2UIlist(cls, xml_str):
        merged_ids = []
        xml_data = Entrez.read(StringIO(xml_str))
        for id in cls.xmldata2UIlist(xml_data):
            merged_ids += [id]
        return merged_ids

    @classmethod
    def xmldata2UIlist(cls, xml_data):
        merged_ids = []

        try:
            # Always prioritize the result links as opposed to the search links
            # elink - retrieves linked IDs for cmd=neighbor|neighbor_score only
            if 'LinkSetDb' in xml_data[0]:
                for lnk in xml_data[0]['LinkSetDb'][0]['Link']:
                    # elink for neighbor
                    if isinstance(lnk, str):
                        merged_ids.append(lnk)
                    # elink for neighbor_score
                    else:
                        merged_ids.append(lnk['Id'])
            # esearch
            elif 'IdList' in xml_data:
                for id in xml_data['IdList']:
                    merged_ids += [id]
        # If it was not elink output, we will end up here
        except Exception:
            # esearch
            if 'IdList' in xml_data:
                for id in xml_data['IdList']:
                    merged_ids += [id]

        return merged_ids

    @classmethod
    def parse_ids(cls, id_list, id, history_file, xml_file, json_file):
        """Parse IDs passed on --cli or in a file passed to the cli
        """
        merged_ids = []
        if id is not None:
            for pid in id.replace('__cn__', ',').replace('\n', ',').split(','):
                if pid is not None and len(pid) > 0:
                    merged_ids.append(pid)

        if id_list is not None:
            with open(id_list, 'r') as handle:
                merged_ids += [x.strip() for x in handle.readlines()]

        if xml_file is not None:
            tmp_ids = cls.xmlfile2UIlist(xml_file)
            for id in tmp_ids:
                merged_ids += [id]

        if json_file is not None:
            tmp_ids = cls.jsonfile2UIlist(json_file)
            for id in tmp_ids:
                merged_ids += [id]

        return merged_ids

    @classmethod
    def getVersion(cls):
        """Return the biopython version
        """
        import Bio
        return Bio.__version__
