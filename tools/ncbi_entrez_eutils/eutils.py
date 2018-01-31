import json
import os

from Bio import Entrez
from six import StringIO

Entrez.tool = "GalaxyEutils_1_0"
BATCH_SIZE = 200


class Client(object):

    def __init__(self, history_file=None, user_email=None, admin_email=None):
        self.using_history = False

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

    def fetch(self, db, ftype=None, **payload):
        os.makedirs("downloads")

        if 'id' in payload:
            summary = self.id_summary(db, payload['id'])
        else:
            summary = self.history_summary(db)

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

    def gquery(self, **kwargs):
        return Entrez.egquery(**kwargs).read()

    def citmatch(self, **kwargs):
        return Entrez.ecitmatch(**kwargs).read()

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
