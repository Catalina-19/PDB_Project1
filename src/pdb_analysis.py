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

    def detect_three_residue_tight_loop_serial(self, pdb_list: list, cutoff: float = 3.5):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k
        :param pdb_list:
        :param cutoff:
        :return: None
        """
        pass

    def detect_three_residue_tight_loop_MPI(self, pdb_list: list, cutoff: float = 3.5):
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
    pass
