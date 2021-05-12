"""Download the pdb structure for a list of pdb_ids"""

# import argparse
import re
import requests
import time
import urllib.request


class Download(object):
    def __init__(self):
        """
        A class to download pdb structures
        """
        self.base_url = "https://www.ebi.ac.uk/pdbe/"
        self.api_base = self.base_url + "api/"
        self.secondary_structure_url = self.api_base + 'pdb/entry/secondary_structure/'

    def make_request(self, url, mode, pdb_id):
        if mode == "get":
            response = requests.get(url=url + pdb_id)

        for _ in range(10):
            response = requests.get(url=url + pdb_id)
            if response.status_code == 200:
                return response.json()
            else:
                time.sleep(5)
        return None

    def get_pdb(self, pdb_id=None):
        if not pdb_id:
            print("no pdb_id")
            return None

        if pdb_id:
            # pdbl = PDBList()
            # pdbl.retrieve_pdb_file(pdb_id)
            urllib.request.urlretrieve('http://files.rcsb.org/download/' + pdb_id + '.pdb', pdb_id + '.pdb')

    def read_list(self, input):
        with open(input, "rt") as f:
            for line in f:
                line = line.split()
                pdb_searcher = re.compile(r'^[0-9][a-zA-Z0-9]{3}')
                pdb_id = line[0]
                if pdb_searcher.match(pdb_id):
                    if len(pdb_id) > 4:
                        pdb_id = pdb_id[0:-1]
                    self.get_pdb(pdb_id=pdb_id)
                else:
                    continue


if __name__ == "__main__":
    pdb_list = "pdb_list.txt"
    download = Download()
    download.read_list(pdb_list)
