# import re
# import argparse
import requests
import time
from Bio.PDB import *

base_url = "https://www.ebi.ac.uk/pdbe/"
api_base = base_url + "api/"
secondary_structure_url = api_base + 'pdb/entry/secondary_structure/'


def make_request(url, mode, pdb_id):
    if mode == "get":
        response = requests.get(url=url + pdb_id)

    for _ in range(10):
        response = requests.get(url=url + pdb_id)
        if response.status_code == 200:
            return response.json()
        else:
            time.sleep(5)
    return None


def get_pdb(pdb_id=None):
    if not pdb_id:
        print("no pdb_id")
        return None

    if pdb_id:
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_id)

if __name__ == "__main__":
    with open("5pdb.txt", "rt") as f:
        for line in f:
            line = line.split()
            pdb_id=line[0]
            get_pdb(pdb_id=pdb_id)
