from prody import *


class Analysis(object):
    def __init__(self):
        """
        A class to perform structural analysis
        """
        pass

    def read_pdb(self, file_object):
        """
        Reads a PDB file and returns a structure object
        :param file_object:
        :return: structure
        """
        structure = parsePDB(file_object)
        return structure

    def read_input(self, input):
        pdb_names = list()
        with open(input, 'r') as f:
            for line in f:
                line = line.split()
                pdb_names.append(line[0])
        return pdb_names

    def detect_three_residue_tight_loop_serial(self, pdb_list: str, cutoff: float = 4.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k
        :param pdb_list:
        :param cutoff:
        :return: None
        """
        # TODO
        # 1 Read the input file with pdb names into a list
        pdb_names = self.read_input(pdb_list)
        # 2 Read pdbs one by one
        for pdb_id in pdb_names:
           struct = self.read_pdb(pdb_id)

        # 3 Compute distance and phi - psi
        hv = struct.getHierView()
        standard_aa = struct.protein.stdaa
        for chain in hv:
            # Get all residues listed into a variable
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
                    if distance_CA_N < cutoff:
                        print(f"{pdb_id:>4}", f"{str(chain):>5}", f"{str(res_i):>5}", f"{str(res_j):>5}",
                              f"{str(res_k):>5}")
                        print("Distance:", end="")
                        print(f"{str(round(distance_CA_N, 2)):>8}")
                        # print("Distance: %.2f" % (distance_CA_N))
                        print("Phi:", end="")
                        print(f"{str(round(calcPhi(residues[index]), 2)):>15}",
                              f"{str(round(calcPhi(residues[index + 1]), 2)):>5}",
                              f"{str(round(calcPhi(residues[index + 2]), 2)):>5}")
                        # print("Phi: %.2f, %.2f, %.2f" % (calcPhi(residues[index]), calcPhi(residues[index+1]),
                        #                                 calcPhi(residues[index+2])))
                        print("Psi:", end="")
                        print(f"{str(round(calcPsi(residues[index]), 2)):>15}",
                              f"{str(round(calcPsi(residues[index + 1]), 2)):>5}",
                              f"{str(round(calcPsi(residues[index + 2]), 2)):>5}")
                        # print("Psi: %.2f, %.2f, %.2f" % (calcPsi(residues[index]), calcPsi(residues[index+1]),
                        #                                 calcPsi(residues[index+2])))

                else:
                    continue

    def detect_three_residue_tight_loop_MPI(self, pdb_list: list, cutoff: float = 4.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k. the computation is done in parallel with distributed
        memory scheme.

        :param pdb_list:
        :param cutoff:
        """
        pass

if __name__ == "__main__":
    pdb_list = "5pdb.txt"
    analysis = Analysis()
    analysis.detect_three_residue_tight_loop_serial(pdb_list)
