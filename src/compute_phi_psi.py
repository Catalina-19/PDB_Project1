""" From the PDB ids of a text file, for each PDB, compute Phi and Psi angles of three residues when the distance
between the C, CA and N atoms of the first and third residues is lower than a cutoff """

from prody import *
from mpi4py import MPI
import re


# Disable ProDy output
confProDy(verbosity='info')


class Analysis(object):
    def __init__(self):
        """
        A class to perform structural analysis
        """
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

    def read_pdb(self, pdb_name):
        """
        Reads a PDB file and returns a structure object
        :param pdb_name:
        :return: structure
        """
        structure = parsePDB(pdb_name)
        return structure

    def read_input(self, input):
        """
        Joins all PDB ids into a list (pdb_names)
        :param input:
        :return: pdb_names
        """
        pdb_names = list()
        with open(input, 'r') as f:
            for line in f:
                line = line.split()
                pdb_searcher = re.compile(r'^[0-9][a-zA-Z0-9]{3}')
                pdb_id = line[0]
                if pdb_searcher.match(pdb_id):
                    pdb_names.append(pdb_id)
                else:
                    continue
        return pdb_names

    def detect_three_residue_tight_loop_serial(self, pdb_list: str, cutoff_1: float = 4.0, cutoff_2: float = 5.0,
                                               cutoff_3: float = 5.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k

        :param pdb_list:
        :param cutoff_1:
        :param cutoff_2:
        :param cutoff_3:
        :return: result
        """
        # TODO
        result = list()
        # 1 Read the input file with pdb names into a list
        # pdb_names = self.read_input(pdb_list)
        # 2 Read pdbs one by one
        for pdb_id in pdb_list:
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
                # resi: C - resk N
                # resi: CA - resk CA
                # resi: N - resk C
                if res_i in standard_aa and res_j in standard_aa and res_k in standard_aa:
                    distance_C_N = calcDistance(residues[index].getAtom("C"), residues[index + 2].getAtom("N"))
                    distance_CA_CA = calcDistance(residues[index].getAtom("CA"), residues[index + 2].getAtom("CA"))
                    distance_N_C = calcDistance(residues[index].getAtom("N"), residues[index + 2].getAtom("C"))
                    # If distance < cutoff=4.0
                    # Compute Phi and psi
                    if distance_C_N < cutoff_1 and distance_CA_CA < cutoff_2 and distance_N_C < cutoff_3:
                        try:
                            Phi_i = calcPhi(residues[index])
                            Phi_j = calcPhi(residues[index + 1])
                            Phi_k = calcPhi(residues[index + 2])
                            Psi_i = calcPsi(residues[index])
                            Psi_j = calcPsi(residues[index + 1])
                            Psi_k = calcPsi(residues[index + 2])
                            result.append(
                                '{:^5s}   {:^8s}   {:^9s}   {:^9s}   {:^9s}   {:^12.3f}   {:^14.3f}   {:^11.3f}   {:^13.3f}   {:^10.3f}   {:^11.3f}   {:^11.3f}   {:^11.3f}   {:^11.3f} \n'.format(
                                    pdb_id, str(chain), str(res_i), str(res_j), str(res_k), distance_C_N,
                                    distance_CA_CA, distance_N_C, Phi_i, Phi_j,
                                    Phi_k, Psi_i, Psi_j, Psi_k))
                        except Exception:
                            continue
                else:
                    continue
        return result

    def detect_three_residue_tight_loop_MPI(self, pdb_list: list, cutoff_1: float = 4.0, cutoff_2: float = 5.0,
                                            cutoff_3: float = 5.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k. the computation is done in parallel with distributed
        memory scheme.

        :param pdb_list:
        :param cutoff_1:
        :param cutoff_2:
        :param cutoff_3:
        """

        # (1)
        #   On Master process (rank 0)
        #       reads input file and get the pdb list
        #       divide the pdb list into a list of sub list
        #       send the sub lists to each clients (rank > 0). Keep the fist element to be processed by master process
        if self.rank == 0:
            data = self.read_input(pdb_list)
            sub_lists = list()
            ave, rem = divmod(len(data), self.size)
            start, end = 0, 0
            for i in range(self.size):
                start = end
                if i < rem:
                    end = start + ave + 1
                else:
                    end = start + ave
                sub_lists.append(data[start:end])

            for i in range(1, self.size):
                self.comm.send(sub_lists[i], dest=i)
            data = sub_lists[0]


            #   On Client processes (rank > 0)
            #       receive the sub list
        if self.rank > 0:
            data = self.comm.recv(source=0)

        # (2)
        #   On all process
        #       Read pdb File
        #       Compute the metrics, save them in a variable
        metrics = self.detect_three_residue_tight_loop_serial(data)

        # (3)
        #   On Master process
        #       Receive the metrics from other process
        #       Print them out
        #   On Clients
        #       Sent the metrics to the master node
        if self.rank == 0:
            results = list()
            results.extend(metrics)
            for i in range(1, self.size):
                results.extend(self.comm.recv(source=i))

        if self.rank > 0:
            self.comm.send(metrics, dest=0)

        # Print results
        if self.rank == 0:
            with open("output.txt", 'w') as f:
                f.write(
                '{:^5}   {:^8}   {:^9}   {:^9}   {:^9}   {:^13}  {:^13}   {:^12}   {:^8}   {:^8}   {:^8}   {:^8}   {:^8}   {:^8}'.format(
                    "PDB", "Chain", "Residue i", "Residue j", "Residue k", "Distance_C_N", "Distance_CA_CA",
                    "Distance_N_C", "Angle_Phi_i", "Angle_Phi_j", "Angle_Phi_k", "Angle_Psi_i", "Angle_Psi_j",
                    "Angle_Psi_k \n",
                ))
            for i in results:
                with open("output.txt", 'a') as f:
                    f.write(i)


if __name__ == "__main__":
    pdb_list = "cullpdb_pc50_res2.0_R0.25_d2021_03_25_chains17980.gz"
    analysis = Analysis()
    # analysis.detect_three_residue_tight_loop_serial(pdb_list)
    analysis.detect_three_residue_tight_loop_MPI(pdb_list)
