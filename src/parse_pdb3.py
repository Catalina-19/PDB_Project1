#Import necessary packages
from prody import *
from pylab import *

import os
# TODO Complete this either with native prody class or manual
standard_aa = ['ALA', 'PHE']

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

                for index in range(len(residues) - 3):
                    res_i = residues[index]
                    res_j = residues[index + 1]
                    res_k = residues[index + 2]

                    if res_i not in standard_aa and res_j not in standard_aa and res_k not in standard_aa:
                        continue
                    # Compute distance between atom C of res_i and atom N of res_k
                    #   if distance < cutoff=4.0
                    #       Compute Phi and psi

                for _ in range(len(chain)-3):
                    #Iterate over all the residues in the chain
                    for i, residue in enumerate(residues[count:count+1]):
                        #Take into account only amino acids, not waters or others.
                        if residue.getResname() != "HOH" and residue.getResname() != "HEM" and residue.getResname() != "CD" and residue.getResname() != "ACT" and residues[count + 3].getResname() != "HOH":
                            try:
                                #Calculate the distance between the C-alphas of the residues two residues away.
                                print("%s: Distance between C-alpha of residues %s and %s is %.2f" % (chain, residues[count], residues[count + 2], calcDistance(residues[count].getAtom("CA"), residues[count + 2].getAtom("CA"))))
                                #Compute the Phi and the Psi of the residues in which the C-alpha distance is shorter than 3.5 Angstrom..Avoid the terminal residues(first and last, in which Phi and Psi angles cannot be calculated).
                                #if calcDistance(residues[count].getAtom("C"), residues[count + 2].getAtom("N")) < 3.5 and residue.getResnum() != 1:
                                if calcDistance(residues[count].getAtom("CA"), residues[count + 2].getAtom("CA")) < 3.5 and residue.getResnum() != 1:
                                    print("%s: the Phi of the residue %s is %.2f" % (chain, residues[count], calcPhi(residues[count])))
                                    print("%s: the Psi of the residue %s is %.2f" % (chain, residues[count], calcPsi(residues[count])))
                            except:
                                continue
                    #Sum 1 to the count, to calculate over the next residue.
                    count += 1