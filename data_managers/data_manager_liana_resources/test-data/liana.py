from pathlib import Path
import pandas as pd

class Resource:
    @staticmethod
    def show_resources():
        return ["consensus", "cellphonedb", "cellchatdb", "connectomedb2020", "mouseconsensus"]
    @staticmethod
    def select_resource(resource_name):
        return pd.DataFrame({"ligand": ["L"], "receptor": ["R"], "resource": [resource_name]})
    @staticmethod
    def get_hcop_orthologs(target_organism="mouse", min_evidence=3, **kwargs):
        return pd.DataFrame({"human_symbol": ["A"], f"{target_organism}_symbol": ["a"]})
    @staticmethod
    def get_metalinks(biospecimen_location=None, **kwargs):
        return pd.DataFrame({"type": ["lr"], "biospecimen_location": [str(biospecimen_location)]})
resource = Resource()
