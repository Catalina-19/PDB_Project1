# Import necessary packages
from prody import *
from pylab import *

import os

# TODO Complete this either with native prody class or manual
# standard_aa_manual = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
# 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# First, read the text file in the directory and create an empty list (pdb_names)
for file in os.listdir('.'):
    if file.endswith('.txt'):
        with open(file, "rt") as f:
            pdb_names = list()
            # Read all the lines in the text file and append them to the pdb_names list
            for line in f:
                line = line.split()
                pdb_names.append(line[0])
            # Read each pdb file, get each structure and get the Hierarchical View of every one. Besides, define all
            # the standard amino acids in a variable.
            for pdb_id in pdb_names:
                structure = parsePDB(pdb_id)
                hv = structure.getHierView()
                standard_aa = structure.protein.stdaa
                # To iterate over every chain
                for chain in hv:
                    # Get all residues listed into a variable, useful later
                    residues = list(chain)
                    for index in range(len(residues) - 3):
                        res_i = residues[index]
                        res_j = residues[index + 1]
                        res_k = residues[index + 2]
                        # Check that all residues are one of the standard residues
                        # Compute distance between atom C of res_i and atom N of res_k
                        if res_i in standard_aa and res_j in standard_aa and res_k in standard_aa:
                            distance_CA_N = calcDistance(residues[index].getAtom("CA"), residues[index + 2].getAtom("N"))
                            # If distance < cutoff=4.0
                            # Compute Phi and psi
                            if distance_CA_N < 4.0:
                                print("%s: Distance between the atoms C-alpha of residue %s and N of residue %s is %.2f"
                                      % (chain, residues[index], residues[index + 2], distance_CA_N))
                                print("%s: the Psi of the residue %s is %.2f" % (
                                    chain, residues[index], calcPsi(residues[index])))
                                print("%s: the Phi of the residue %s is %.2f" % (
                                    chain, residues[index], calcPhi(residues[index])))
                        else:
                            continue
