from prody import *
from mpi4py import MPI
# import time

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
                    distance_CA_N = calcDistance(residues[index].getAtom("C"), residues[index + 2].getAtom("N"))
                    # If distance < cutoff=4.0
                    # Compute Phi and psi
                    if distance_CA_N < cutoff:
                        try:
                            # I have changed it so it does not print anything unless wanted.
                            Phi_i = calcPhi(residues[index])
                            Phi_j = calcPhi(residues[index + 1])
                            Phi_k = calcPhi(residues[index + 2])
                            Psi_i = calcPsi(residues[index])
                            Psi_j = calcPsi(residues[index + 1])
                            Psi_k = calcPsi(residues[index + 2])
                            result.append('{:5s}   {:8s}   {:8s}   {:8s}   {:8s}   {:10.3f}   {:10.3f}   {:10.3f}   {:10.3f}   {:10.3f}   {:10.3f}   {:10.3f}'.format(pdb_id, str(chain), str(res_i), str(res_j), str(res_k), distance_CA_N, Phi_i, Phi_j, Phi_k, Psi_i, Psi_j, Psi_k))
                        except Exception:
                            continue

                    # If wanted to print nicely
                    # if distance_CA_N < cutoff:
                    #     print(f"{pdb_id:>4}", f"{str(chain):>5}", f"{str(res_i):>5}", f"{str(res_j):>5}",
                    #           f"{str(res_k):>5}")
                    #     print("Distance:", end="")
                    #     print(f"{str(round(distance_CA_N, 2)):>8}")
                    #     # print("Distance: %.2f" % (distance_CA_N))
                    #     print("Phi:", end="")
                    #     print(f"{str(round(calcPhi(residues[index]), 2)):>15}",
                    #           f"{str(round(calcPhi(residues[index + 1]), 2)):>5}",
                    #           f"{str(round(calcPhi(residues[index + 2]), 2)):>5}")
                    #     # print("Phi: %.2f, %.2f, %.2f" % (calcPhi(residues[index]), calcPhi(residues[index+1]),
                    #     #                                 calcPhi(residues[index+2])))
                    #     print("Psi:", end="")
                    #     print(f"{str(round(calcPsi(residues[index]), 2)):>15}",
                    #           f"{str(round(calcPsi(residues[index + 1]), 2)):>5}",
                    #           f"{str(round(calcPsi(residues[index + 2]), 2)):>5}")
                    #     # print("Psi: %.2f, %.2f, %.2f" % (calcPsi(residues[index]), calcPsi(residues[index+1]),
                    #     #                                 calcPsi(residues[index+2])))

                else:
                    continue
        return result

    def detect_three_residue_tight_loop_MPI(self, pdb_list: list, cutoff: float = 4.0):
        """
        Detects tight three residue-turns with the defined distance cutoff and print
        the pdb name, res-i, res-j, res-k, distance, phi-i, phi-j, phi-k,
        psi-i, psi-j, psi-k. the computation is done in parallel with distributed
        memory scheme.

        :param pdb_list:
        :param cutoff:
        """

        # (1)
        #   On Master process (rank 0)
        #       reads input file and get the pdb list
        #       divide the pdb list into a list of sub list
        #       send the sub lists to each clients (rank > 0). Keep the fist element to be processed by master process
        print('start on: ', self.rank)
        if self.rank == 0:
            data = self.read_input(pdb_list)
            #sub_lists = [data[x:x + len(data) // self.size] for x in range(0, len(data), len(data) // self.size)]
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


            print('Rank', self.rank, 'Starting sending', data)
            print('Rank', self.rank, 'Starting sending', sub_lists)

            for i in range(1, self.size):
                self.comm.send(sub_lists[i], dest=i)
            data = sub_lists[0]
            print('Rank', self.rank, 'finished sending')

            #   On Client processes (rank > 0)
            #       receive the sub list
        if self.rank > 0:
            print('Rank', self.rank, 'Starting reciveing')
            data = self.comm.recv(source=0)
            print('Rank', self.rank, 'Finished reciveing')

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

        # I have tried putting here a barrier in order to synchronise all processes, but it dit not help
        # comm.Barrier()
        # I have also tried to define a function to finish communication, but it did not result either
        # def finish_comm():

        if self.rank == 0:
            results = list()
            results.extend(metrics)
            #print(metrics)
            #print(results)
            for i in range(1, self.size):
                results.extend(self.comm.recv(source=i))

        if self.rank > 0:
            self.comm.send(metrics, dest=0)


        # finish_comm()
        if self.rank == 0:
            for i in results:
                print(i)


if __name__ == "__main__":
    pdb_list = "5pdb.txt"
    analysis = Analysis()
    # analysis.detect_three_residue_tight_loop_serial(pdb_list)
    analysis.detect_three_residue_tight_loop_MPI(pdb_list)
