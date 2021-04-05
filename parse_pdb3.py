#Import necessary packages
from prody import *
from pylab import *

import os

#First, read all PDB files in the directory
for file in os.listdir('.'):
    if file.endswith('.pdb'):
        with open(file, "rt") as f:
            #Read the structure and get the Hierarchical View
            structure = parsePDB(file)
            hv = structure.getHierView()
            #To iterate over every chain
            for chain in hv:
                #Get all residues listed into a variable, useful later
                residues = list(chain)
                #Set count to 0, in order to parse the residues in a window of 3
                count = 0
                for _ in range(len(chain) // 3):
                    #Iterate over all the residues in the chain
                    for i, residue in enumerate(residues[count:count+1]):
                        #Take into account only amino acids, not waters or others. Avoid the terminal residues (first and last, in which Phi and Psi angles cannot be calculated).
                        if residue.getResname() != "HOH" and residue.getResname() != "HEM" and residue.getResname() != "CD" and residue.getResname() != "ACT" and residues[count + 3].getResname() != "HOH" and residue.getResnum() != 1:
                            try:
                                #Compute the Phi and the Psi for the desired residues.
                                print("%s: the Phi of the residue %s is %.2f" % (chain, residue, calcPhi(residue)))
                                print("%s: the Psi of the residue %s is %.2f" % (chain, residue, calcPsi(residue)))
                            except:
                                continue
                    #Sum 3 to the count, to calculate over the next 3 residues.
                    count += 3
            for chain in hv:
                residues = list(chain)
                count = 0
                for _ in range(len(chain) // 3):
                    for i, residue in enumerate(residues[count:count + 1]):
                        if residue.getResname() != "HOH" and residue.getResname() != "HEM" and residue.getResname() != "CD" and residue.getResname() != "ACT" and residues[count + 3].getResname() != "HOH":
                            try:
                                #Calculate the distance between the C-alphas of the residues two residues away.
                                print("%s: Distance between C-alpha of residues %s and %s is %.2f" % (chain, residues[count], residues[count + 2],calcDistance(residues[count].getAtom("CA"), residues[count + 2].getAtom("CA"))))
                            except:
                                continue
                    count += 3