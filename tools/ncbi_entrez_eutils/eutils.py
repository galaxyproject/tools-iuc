import os
import json
import StringIO
from Bio import Entrez
Entrez.email = os.environ.get('NCBI_EUTILS_CONTACT', None)
Entrez.tool = "GalaxyEutils"
BATCH_SIZE = 200


class Client(object):

    def __init__(self, history_file=None):
        self.using_history = False
        if Entrez.email is None:
            raise Exception("Cannot continue without an email; please set "
                            "administrator email in NCBI_EUTILS_CONTACT")

        if history_file is not None:
            with open(history_file, 'r') as handle:
                data = json.loads(handle.read())
                self.query_key = data['QueryKey']
                self.webenv = data['WebEnv']
                self.using_history = True

    def get_history(self):
        if not self.using_history:
            return {}
        else:
            return {
                'query_key': self.query_key,
                'WebEnv': self.webenv,
            }

    def post(self, database, **payload):
        return json.dumps(Entrez.read(Entrez.epost(database, **payload)), indent=4)

    def fetch(self, db, whole=False, **payload):
        if whole:
            if 'id' in payload:
                summary = self.id_summary(db, payload['id'])
            else:
                summary = self.history_summary(db)

            count = len(summary)

            payload['retmax'] = BATCH_SIZE

            # Print the first one
            print Entrez.efetch(db, **payload).read()
            # Then write subsequent to files for <discover datasets>
            for i in range(BATCH_SIZE, count, BATCH_SIZE):
                payload['retstart'] = i
                # TODO: output multiple files??? Collection?
                with open('%s.out' % i, 'w') as handle:
                    handle.write(Entrez.efetch(db, **payload).read())
        else:
            print Entrez.efetch(db, **payload).read()

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

    def extract_history(self, xml_data):
        parsed_data = Entrez.read(StringIO.StringIO(xml_data))
        history = {}
        for key in ('QueryKey', 'WebEnv'):
            if key in parsed_data:
                history[key] = parsed_data[key]

        return history

    def search(self, **payload):
        return Entrez.esearch(**payload).read()

    def info(self, **kwargs):
        return Entrez.einfo(**kwargs).read()

    @classmethod
    def parse_ids(cls, id_list, id, history_file):
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

        # Exception hanlded here for uniformity
        if len(merged_ids) == 0 and history_file is None:
            raise Exception("Must provide history file or IDs")

        return merged_ids
