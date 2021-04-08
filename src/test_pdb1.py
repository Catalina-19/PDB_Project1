# import re
import requests
import time
import argparse


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
    # This should be changed to account for the cases that the server is busy or there is problem to
    # establish connection. something like:
    """
    for _ in range(10):
        response = requests.get(url=url + pdb_id)
        if response.status_code == 200:
            return response.json()
        else:
            time.sleep(5)
    
    # If after 10 loop didn't get anything
    return None
    """


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

# Add this line so you can use your code both as a module (to be imported to other codes) and as a
# stand alone program. If you call the code as a program (python PDB_project1) the code after if __name__ ...
# will be executed. if you import with out if __name__ ... line, the code will be executed at once.
if __name__ == "__main__":
    # do not hard code file names, parameters etc. Use argparse (I imported for you) to read the parameters from
    # commad line. Read more about argparse it is a good one
    """
    A example from my stuff
    parser = argparse.ArgumentParser(description='Active Site Design')
    parser.add_argument('conf', type=str, help='Configuration file.')
    parser.add_argument('-debug', action='store_true', help='Debug mode.')

    args = parser.parse_args()
    """

    # use with cluase. It is recommended and takes care of closing the file
    with open("5pdb.txt", "rt") as f:
        # trim the line. Dont trust the input given by the user
        for line in f:
            line = line.split()
            # if the id is the first field (it should be 0)
            pdb_id = line[0]
            get_secondary_structure_ranges(pdb_id=pdb_id)

    """
    Also your code is suppose to get the pdb files not the secondary structure annotation.
    """