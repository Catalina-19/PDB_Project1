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
        pass

    def read_input(self, input):
        pdb_names = list()
        with open(input, 'r'):
            pass
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
        for pdb_name in pdb_names:
           struct = self.read_pdb(pdb_name)

        # 3 Compute distance and phi - psi

    def detect_three_residue_tight_loop_MPI(self, pdb_list: list, cutoff: float = 4.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phik,
        psi-i, psi-j, psi-k. the computation is done in parallel with distributed
        memory scheme.

        :param pdb_list:
        :param cutoff:
        """
        pass

if __name__ == "__main__":
    pdb_list = "input_file"
    analysis = Analysis()
    analysis.detect_three_residue_tight_loop_serial(pdb_list)
