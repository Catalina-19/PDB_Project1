# import re
import requests

base_url = "https://www.ebi.ac.uk/pdbe/"
api_base = base_url + "api/"
secondary_structure_url = api_base + 'pdb/entry/secondary_structure/'


def make_request(url, mode, pdb_id):
    if mode == "get":
        response = requests.get(url=url + pdb_id)

    if response.status_code == 200:
        return response.json()
    else:
        print("[No data retrieved - %s] %s" % (response.status_code, response.text))
    return None


def get_secondary_structure_ranges(pdb_id=None, pdb_list=None):
    if not pdb_id:
        print("no pdb_id")
        return None

    if pdb_id:
        data = make_request(secondary_structure_url, "get", pdb_id)

    for entry_id in data.keys():
        entry = data[entry_id]
        molecules = entry["molecules"]

        for i in range(len(molecules)):
            chains = molecules[i]["chains"]

            for j in range(len(chains)):
                secondary_structure = chains[j]["secondary_structure"]
                helices = secondary_structure["helices"]
                strands = secondary_structure["strands"]
                helix_list = []
                strand_list = []
                for k in range(len(helices)):
                    start = helices[k]["start"]["residue_number"]
                    end = helices[k]["end"]["residue_number"]
                    helix_list.append("%s-%s" % (start, end))

                for l in range(len(strands)):
                    start = strands[l]["start"]["residue_number"]
                    end = strands[l]["end"]["residue_number"]
                    strand_list.append("%s-%s" % (start, end))

                report = "%s chain %s has " % (entry_id, chains[j]["chain_id"])
                if len(helix_list) > 0:
                    report += "helices at residue ranges %s " % str(helix_list)
                else:
                    report += "no helices "
                report += "and "
                if len(strand_list) > 0:
                    report += "strands at %s" % str(strand_list)
                else:
                    "no strands"
                print(report)
    return None


f = open("5pdb.txt", "rt")
for line in f:
    get_secondary_structure_ranges(pdb_id=line)

f.close()
