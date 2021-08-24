""" From the PDB ids of a text file, for each PDB, compute Phi and Psi angles of three residues when the distance
between the C, CA and N atoms of the first and third residues is lower than a cutoff """
import math

# import nglview.color
# import simtk.unit
import pickle

from prody import *
from mpi4py import MPI
import re
# from pdb_download_structure import Download
import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
import MDAnalysis.analysis.align as align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import nglview as nv
import seaborn as sns
from MDAnalysis.analysis import rms

# Disable ProDy output
confProDy(verbosity='debug')


class Analysis(object):
    def __init__(self):
        """
        A class to perform structural analysis
        """
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.pdb_id_list = list()

    def read_cif(self, pdb_name):
        """
        Reads a PDB file and returns a structure object
        :param pdb_name:
        :return: structure
        """
        # structure = parsePDB(pdb_name)
        try:
            structure = parseMMCIF(pdb_name)
            return structure
        except:
            return None

    def read_input(self, input):
        """
        Joins all PDB ids into a list (pdb_names)
        :param input:
        :return: pdb_names
        """
        pdb_names = list()
        with open(input, 'r') as f:
            for line in f:
                # print(line)
                line = line.split()
                if len(line) >= 1:
                    # pdb_searcher = re.compile(r'^[0-9][a-zA-Z0-9]{3}')
                    pdb_names.append(line[0])
                # if pdb_searcher.match(pdb_id):
                #    pdb_names.append(pdb_id)
                else:
                    print('Could not read line {}'.format(line))
                #    continue
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
            # print(self.rank, pdb_id)
            struct = self.read_cif(pdb_id)
            if struct is None:
                print('Failed to open pdb {}'.format(pdb_id), flush=True)
                continue

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
                        try:
                            distance_C_N = calcDistance(residues[index].getAtom("C"), residues[index + 2].getAtom("N"))
                            distance_CA_CA = calcDistance(residues[index].getAtom("CA"),
                                                          residues[index + 2].getAtom("CA"))
                            distance_N_C = calcDistance(residues[index].getAtom("N"), residues[index + 2].getAtom("C"))
                        except Exception as a:
                            print(pdb_id, str(chain), str(res_i), str(res_j), str(res_k))
                            continue
                        # If distance < cutoff=4.0
                        # Compute Phi and psi
                        # print(self.rank, distance_C_N, distance_CA_CA, distance_N_C)
                        if distance_C_N < cutoff_1 and distance_CA_CA < cutoff_2 and distance_N_C < cutoff_3:
                            try:
                                Phi_i = calcPhi(residues[index])
                                Phi_j = calcPhi(residues[index + 1])
                                Phi_k = calcPhi(residues[index + 2])
                                Psi_i = calcPsi(residues[index])
                                Psi_j = calcPsi(residues[index + 1])
                                Psi_k = calcPsi(residues[index + 2])
                                string = '{:^5s}   {:^8s}   {:^9s}   {:^9s}   {:^9s}   {:^12.3f}   {:^14.3f}   {:^11.3f}   {:^13.3f}   {:^10.3f}   {:^11.3f}   {:^11.3f}   {:^11.3f}   {:^11.3f} \n'.format(
                                    pdb_id, str(chain), str(res_i), str(res_j), str(res_k), distance_C_N,
                                    distance_CA_CA, distance_N_C, Phi_i, Phi_j,
                                    Phi_k, Psi_i, Psi_j, Psi_k)
                                result.append(string)
                                print('###  ', string, flush=True)
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

            # print(sub_lists)
            for l in zip(sub_lists[0], sub_lists[1], sub_lists[2], sub_lists[3]):
                print(l)

            for i in range(1, self.size):
                self.comm.send(sub_lists[i], dest=i)
            data = sub_lists[0]

            #   On Client processes (rank > 0)
            #       receive the sub list
        if self.rank > 0:
            data = self.comm.recv(source=0)

        print(str([l for l in data]))

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

    def filtering_redundant_results(self):
        """
        Filter in order to remove the redundant (repeated) motifs
        :return:
        """
        with open("output.txt", 'rt') as f:
            filtered_list = list()
            unfiltered_pdb_list = list()
            unfiltered_res_i = list()
            for i in f:
                line = i.split()
                pdb_id = line[0]
                res_i = line[3] + ' ' + line[4]
                unfiltered_pdb_list.append(pdb_id)
                unfiltered_res_i.append(res_i)

            output_file = open("output.txt")

            for position, line in enumerate(output_file):
                if unfiltered_pdb_list[position] != unfiltered_pdb_list[position - 1] and unfiltered_pdb_list[position] != unfiltered_pdb_list[position - 2]:
                    filtered_list.append(line)

                elif unfiltered_pdb_list[position] == unfiltered_pdb_list[position - 1] and unfiltered_pdb_list[
                    position] == unfiltered_pdb_list[position - 2] and unfiltered_res_i[position] != unfiltered_res_i[position - 1] and unfiltered_res_i[position] != unfiltered_res_i[position - 2]:
                    filtered_list.append(line)

        with open("non_redundant_output.txt", 'w') as f:
            pass
        for i in filtered_list:
            with open("non_redundant_output.txt", 'a') as f:
                f.write(i)

    def filter_surface_results(self):
        """
        Compute the neighbors for each motif, in order to know whether they are on the surface
        :return:
        """
        neighbors_list = list()
        main_structure = parsePDB("1gk2")
        standard_aa = main_structure.protein.stdaa
        with open("correct_non_redundant.txt", 'r') as f:
            for i in f:
                line = i.split()
                if len(line[0]) == 3:
                    continue
                pdb_id = line[0]
                chain = line[2]
                # res_i = '{} {}'.format(line[3], line[4]).capitalize()
                # res_j = '{} {}'.format(line[5], line[6]).capitalize()
                # res_k = '{} {}'.format(line[7], line[8]).capitalize()
                try:
                    number_res_i = int(line[4])
                except:
                    long_number_res_i = line[4]
                    number_res_i = int(long_number_res_i[:-1])
                # number_res_j = int(line[6])
                try:
                    number_res_k = int(line[8])
                except:
                    long_number_res_k = line[8]
                    number_res_k = int(long_number_res_k[:-1])

                try:
                    structure = parsePDB(pdb_id)
                    # hv = structure.getHierView()
                    # standard_aa_list = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
                    contacts_residues = structure.select(
                        'protein and chain %s and within 5 of resnum %s:%s' % (chain, number_res_i, number_res_k))
                    contacts_list = list()
                    for i in contacts_residues:
                        if i.getResnum() not in contacts_list:
                            contacts_list.append(i.getResnum())
                    neighbors_list.append(len(contacts_list))
                    #
                    # u = mda.Universe(pdb_id + '.pdb')
                    # ag = u.select_atoms("resid %s:%s" % (number_res_i, number_res_k))
                    # neighbors = mda.lib.NeighborSearch.AtomNeighborSearch(atom_group=ag).search(atoms=ag,radius=20.0)
                    # count = 0
                    # for neighbor in neighbors:
                    #     if neighbor in standard_aa:
                    #         count += 1
                    # neighbors_list.append(count)
                except:
                    print("{} could not compute neighbors".format(pdb_id))
                    neighbors_list.append(-1)
        # Append the neighbors count to a new file with all the other information
        input_file = open("correct_non_redundant.txt", 'r')
        with open("neighbors_non_redundant.txt", 'w') as f:
            index = 0
            for line in input_file:
                if index == 0:
                    f.write('{}   {} \n'.format(line[:-1], "NEIGHBORS_COUNT"))
                elif index > 0:
                    f.write('{}    {:^15s} \n'.format(line[:-1], str(neighbors_list[index - 1])))
                index += 1

    def compare_phi_difference(self):
        # Compare the phi angles of each loop with respect to the main structure, the mutant 1gk2
        main_structure = parsePDB("1gk2")
        hv = main_structure.getHierView()
        main_serine = hv['A', 143]
        serine_phi = calcPhi(main_serine)
        serine_psi = calcPsi(main_serine)
        distances_list = list()
        with open("neighbors_non_redundant.txt", 'r') as f:
            for i in f:
                line = i.split()
                pdb_id = line[0]
                if len(pdb_id) > 3:
                    try:
                        # distance = float(math.sqrt( ((float(line[12]) - float(serine_phi))**2 + (float(line[15]) - float(serine_psi))**2))/2 )
                        distance = float(math.sqrt(((float(line[12]) - float(serine_phi)) ** 2 + (
                                    float(line[15]) - float(serine_psi)) ** 2)) / 2)
                        distances_list.append(distance)
                    except:
                        print("{} could not compare phi angles".format(pdb_id))
                        distances_list.append(-1)

        # Append the phi difference to the output in the last column
        input_file = open("neighbors_non_redundant.txt", 'r')
        with open("final_filtered_output.txt", 'w') as f:
            index = 0
            for line in input_file:
                if index == 0:
                    f.write('{}   {} \n'.format(line[:-1], "PHI_DIFFERENCE"))
                elif index > 0:
                    f.write('{}    {:^15.3f} \n'.format(line[:-1], distances_list[index - 1]))
                index += 1

    def pca_analysis_first_part(self, structure="1GK2-system.pdb", traj=["md01.nc", "md02.nc"]):
        # print("here starts function first part")
        u = mda.Universe(structure, traj, dt=2)
        r = mda.Universe(structure)

        # ALIGN. The trajectory should be aligned to the atom group for better results.
        aligner = align.AlignTraj(u, r, select='name CA', in_memory=True).run()
        # print("finshed alignment first part")

        # RUN PCA. It is performed with the CA atoms in order to reduce the computational requirements.
        calculated_pca = pca.PCA(u, select='name CA', n_components=None, align=True).run()
        # print("finished pca first part")

        # number_ag = len(atom_group)
        # print('There are {} backbone atoms in the analysis'.format(number_ag))

        # SAVE. Save the trajectory to a file
        atom_group = u.select_atoms('name CA')
        with mda.Writer("trajectory.xtc", atom_group.n_atoms) as W:
            for ts in u.trajectory:
                W.write(atom_group)
        print("finished writing trajectory first part")

        # Check the cumulated variance of every component
        # calculated_pca.cumulated_variance
        # print(pc.p_components.shape)

        # VISUALIZE VARIANCE
        # plt.ioff()
        plt.plot(calculated_pca.cumulated_variance[:100])
        plt.title('Cumulative variance of first 100 PCs')
        plt.xlabel('Principal component')
        plt.ylabel('Cumulative variance')
        plt.savefig('variance.png', bbox_inches='tight')

        # Save the PCA vectors
        pickle.dump(calculated_pca.p_components, open("pcomponents_wt.pkl", "wb"))

    def pca_analysis_second_part(self, structure="1GK2-system.pdb", traj=["md01.nc", "md02.nc"], number_of_pcs=2):
        u = mda.Universe(structure, traj, dt=2)
        r = mda.Universe(structure)

        # ALIGN. The trajectory should be aligned to the atom group for better results.
        aligner = align.AlignTraj(u, r, select='name CA', in_memory=True).run()
        print("finshed alignment first part", flush=True)

        # RUN PCA. It is performed with the CA atoms in order to reduce the computational requirements.
        calculated_pca = pca.PCA(u, select='name CA', n_components=None, align=False).run()
        print("finished pca second part", flush=True)

        atom_group = u.select_atoms('name CA')
        # Transform the atom group into the weight of each over the pc.
        transformed = calculated_pca.transform(atom_group, n_components=number_of_pcs)
        print("transformed second part", flush=True)
        # calculated_pca.transform()

        # PROJECTIONS. Project the original trajectory onto each of the pcs we are interested in to visualize the motion of each.
        for i in range(number_of_pcs):
            pc = calculated_pca.p_components[:, i]
            trans = transformed[:, i]
            trans.sort()
            projected_trajectory = np.outer(trans, pc) + calculated_pca.mean
            coordinates = projected_trajectory.reshape(len(trans), -1, 3)

            proj1 = mda.Merge(atom_group)
            proj1.load_new(coordinates, order="fac")
            # pca_atoms = u.select_atoms('name CA')
            with mda.Writer("mt1_projected_trajectory_pc{}.xtc".format(i + 1), proj1.atoms.n_atoms) as W:
                for ts in proj1.trajectory:
                    W.write(proj1.atoms)
        print("written trajectory second part")

    def project_mt_pc_onto_wt(self, wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], mt_structure="1GK2-mutant-system.pdb", mt_traj=["md01_mt.nc", "md02_mt.nc"], number_of_pcs=2):
        w = mda.Universe(wt_structure, wt_traj, dt=2)
        r_wt = mda.Universe(wt_structure)
        m = mda.Universe(mt_structure, mt_traj, dt=2)
        r_mt = mda.Universe(mt_structure)
        atom_group_wt = w.select_atoms('name CA')
        atom_group_mt = m.select_atoms('name CA')
        aligner_wt = align.AlignTraj(w, r_wt, select='name CA', in_memory=True).run()
        aligner_mt = align.AlignTraj(m, r_mt, select='name CA', in_memory=True).run()
        wt_calculated_pca = pca.PCA(w, select='name CA', n_components=None, align=False).run()
        mt_calculated_pca = pca.PCA(m, select='name CA', n_components=None, align=False).run()
        print("first part")

        transformed = mt_calculated_pca.transform(atom_group_mt, n_components=2)

        for i in range(2):
            pc = mt_calculated_pca.p_components[:, i]
            trans = transformed[:, i]
            projected = np.outer(trans, pc) + mt_calculated_pca.mean
            coordinates = projected.reshape(len(trans), -1, 3)

            proj1 = mda.Merge(atom_group_mt)
            proj1.load_new(coordinates, order="fac")
            with mda.Writer("wtmt1_projected_trajectory_pc{}.xtc".format(i + 1), proj1.atoms.n_atoms) as W:
                for ts in proj1.trajectory:
                    W.write(proj1.atoms)
        print("end")

    def pc1_vs_pc2_wt(self, wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], number_of_pcs=2):
        print("start")
        w = mda.Universe(wt_structure, wt_traj, dt=2)
        r_wt = mda.Universe(wt_structure)
        atom_group_wt = w.select_atoms('name CA')
        aligner_wt = align.AlignTraj(w, r_wt, select='name CA', in_memory=True).run()
        print("start pca")
        wt_calculated_pca = pca.PCA(w, select='name CA', n_components=None, align=False).run()
        print("finished pca")

        transformed_wt = wt_calculated_pca.transform(atom_group_wt, n_components=number_of_pcs)

        print("start df")
        df_wt = pd.DataFrame(transformed_wt, columns=['PC{}'.format(i+1) for i in range(2)])
        df_wt['Time (ps)'] = df_wt.index * w.trajectory.dt
        df_wt.to_csv('wt_pc1pc2_data.csv')

        print("first plot")
        g = sns.PairGrid(df_wt, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df_wt)))
        g.map(plt.scatter, marker=".")
        plt.savefig("wt_pc1vspc2.png")

        print("finished")

    def pc1_vs_pc2_mt1(self, wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], mt1_structure="mt1_1GK2-system.pdb", mt1_traj=["mt1_md01.nc", "mt1_md02.nc"], number_of_pcs=2):
        print("start")
        w = mda.Universe(wt_structure, wt_traj, dt=2)
        r_wt = mda.Universe(wt_structure)
        atom_group_wt = w.select_atoms('name CA')
        aligner_wt = align.AlignTraj(w, r_wt, select='name CA', in_memory=True).run()
        print("start pca")
        wt_calculated_pca = pca.PCA(w, select='name CA', n_components=None, align=False).run()
        print("finished pca")

        m1 = mda.Universe(mt1_structure, mt1_traj, dt=2)
        r_mt1 = mda.Universe(mt1_structure)
        atom_group_mt1 = m1.select_atoms('name CA')
        aligner_mt1 = align.AlignTraj(m1, r_mt1, select='name CA', in_memory=True).run()

        transformed_mt1 = wt_calculated_pca.transform(atom_group_mt1, n_components=number_of_pcs)

        df_mt1 = pd.DataFrame(transformed_mt1, columns=['PC{}'.format(i + 1) for i in range(2)])
        df_mt1['Time (ps)'] = df_mt1.index * m1.trajectory.dt
        df_mt1.to_csv('mt1_pc1pc2_data.csv')
        print("second plot")

        g = sns.PairGrid(df_mt1, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df_mt1)))
        g.map(plt.scatter, marker=".")
        plt.savefig("mt1_pc1vspc2.png")

        print("finished")

    def pc1_vs_pc2_mt2(self, wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], mt2_structure="mt2_1GK2-system.pdb", mt2_traj=["mt2_md01.nc", "mt2_md02.nc"], number_of_pcs=2):
        print("start")
        w = mda.Universe(wt_structure, wt_traj, dt=2)
        r_wt = mda.Universe(wt_structure)
        atom_group_wt = w.select_atoms('name CA')
        aligner_wt = align.AlignTraj(w, r_wt, select='name CA', in_memory=True).run()
        print("start pca")
        wt_calculated_pca = pca.PCA(w, select='name CA', n_components=None, align=False).run()
        print("finished pca")

        m2 = mda.Universe(mt2_structure, mt2_traj, dt=2)
        r_mt2 = mda.Universe(mt2_structure)
        atom_group_mt2 = m2.select_atoms('name CA')
        aligner_mt2 = align.AlignTraj(m2, r_mt2, select='name CA', in_memory=True).run()

        transformed_mt2 = wt_calculated_pca.transform(atom_group_mt2, n_components=number_of_pcs)

        print("start df")
        df_mt2 = pd.DataFrame(transformed_mt2, columns=['PC{}'.format(i + 1) for i in range(2)])
        df_mt2['Time (ps)'] = df_mt2.index * m2.trajectory.dt
        df_mt2.to_csv('mt2_pc1pc2_data.csv')

        print("third plot")
        g = sns.PairGrid(df_mt2, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df_mt2)))
        g.map(plt.scatter, marker=".")
        plt.savefig("mt2_pc1vspc2.png")

        print("finished")

    def calculate_rmsd(self, wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], mt1_structure="mt1_1GK2-system.pdb", mt1_traj=["mt1_md01.nc","mt1_md02.nc"], mt2_structure="mt2_1GK2-system.pdb", mt2_traj=["mt2_md01.nc","mt2_md02.nc"]):
        wt_u = mda.Universe(wt_structure, wt_traj)
        wt_r = mda.Universe(wt_structure)
        mt1_u = mda.Universe(mt1_structure, mt1_traj)
        mt1_r = mda.Universe(mt1_structure)
        mt2_u = mda.Universe(mt2_structure, mt2_traj)
        mt2_r = mda.Universe(mt2_structure)
        # CALCULATE RMSD
        # Declare the groups to study the RMSD
        # loop_rmsd = 'name CA and around 5 resid 140-146'
        # backbone_rmsd = 'backbone'
        print("start rmsd")
        wt_R_rmsd = rms.RMSD(wt_u,  # universe to align
                          wt_u,  # reference universe or atomgroup
                          select='name CA and around 5 resid 140-146')  # group to superimpose and calculate RMSD
                          # ref_frame=0) # frame index of the reference
                          # groupselections=['name CA and around 5 resid 140-146'],  # groups for RMSD

        wt_R_rmsd.run()
        print("first rmsd run")
        mt1_R_rmsd = rms.RMSD(mt1_u,mt1_u,select='name CA and around 5 resid 140-146')
        mt1_R_rmsd.run()

        mt2_R_rmsd = rms.RMSD(mt2_u,mt2_u,select='name CA and around 5 resid 140-146')
        mt2_R_rmsd.run()

        # Plot the data
        df = pd.DataFrame(wt_R_rmsd.rmsd,columns=['Frame', 'Time (fs)','Loop WT'])

        mt1_df = pd.DataFrame(mt1_R_rmsd.rmsd,columns=['Frame', 'Time (fs)','Loop MT1'])
        mt2_df = pd.DataFrame(mt2_R_rmsd.rmsd,columns=['Frame', 'Time (fs)','Loop MT2'])
        df['Loop MT1'] = mt1_df[['Loop MT1']]
        df['Loop MT2'] = mt2_df[['Loop MT2']]
        df.to_csv('wt_rmsd_data.csv')

        print("picture")
        f = plt.figure(1)
        ax = df.plot(x='Frame', y=['Loop WT', 'Loop MT1', 'Loop MT2'], kind='line')
        ax.set_ylabel(r'RMSD ($\AA$)')
        plt.savefig('wt_rmsd_loop.png')

    def calculate_rmsf(self,wt_structure="1GK2-system.pdb", wt_traj=["md01.nc", "md02.nc"], mt1_structure="mt1_1GK2-system.pdb", mt1_traj=["mt1_md01.nc","mt1_md02.nc"], mt2_structure="mt2_1GK2-system.pdb", mt2_traj=["mt2_md01.nc","mt2_md02.nc"]):
        wt_u = mda.Universe(wt_structure, wt_traj)
        wt_r = mda.Universe(wt_structure)
        # mt1_u = mda.Universe(mt1_structure, mt1_traj)
        # mt1_r = mda.Universe(mt1_structure)
        mt2_u = mda.Universe(mt2_structure, mt2_traj)
        mt2_r = mda.Universe(mt2_structure)
        # CALCUALTE RMSF
        # Align to the first frame, and then average the coordinates
        wt_average = align.AverageStructure(wt_u, wt_u, select='protein and name CA',ref_frame=0).run()
        # mt1_average = align.AverageStructure(mt1_u, mt1_u, select='protein and name CA', ref_frame=0).run()
        mt2_average = align.AverageStructure(mt2_u, mt2_u, select='protein and name CA', ref_frame=0).run()
        print("start align")
        # ref = average.universe
        # Align the trajectory to the average conformation
        # The trajectory can be saved in memory or into a file
        wt_aligner = align.AlignTraj(wt_u, wt_u,
                                  select='protein and name CA',
                                  in_memory=True).run()
        # mt1_aligner = align.AlignTraj(mt1_u, mt1_u,select='protein and name CA',in_memory=True).run()
        mt2_aligner = align.AlignTraj(mt2_u, mt2_u, select='protein and name CA', in_memory=True).run()
        # aligner = align.AlignTraj(u, ref,
        #                           select='protein and name CA',
        #                           filename='aligned_traj.dcd',
        #                           in_memory=False).run()
        # u = mda.Universe(PSF, 'aligned_traj.dcd')
        print("finish align")

        # Now that the trajectory is fitted to the reference, the RMSF can be calculated.
        wt_loop = wt_u.select_atoms('protein and name CA and around 5 resid 140-146')
        # mt1_loop = mt1_u.select_atoms('protein and name CA and around 5 resid 140-146')
        mt2_loop = mt2_u.select_atoms('protein and name CA and around 5 resid 140-146')
        print("start rmsf")
        wt_R_rmsf_loop = rms.RMSF(wt_loop).run()
        # mt1_R_rmsf_loop = rms.RMSF(mt1_loop).run()
        mt2_R_rmsf_loop = rms.RMSF(mt2_loop).run()
        print("finished rmsf")

        # Save the data to a file
        print("start dataframes")
        wt_df_rmsf = pd.DataFrame(wt_R_rmsf_loop.rmsf, columns=['Loop WT'])
        # mt1_df_rmsf = pd.DataFrame(mt1_R_rmsf_loop.rmsf, columns=['Loop MT1'])
        mt2_df_rmsf = pd.DataFrame(mt2_R_rmsf_loop.rmsf, columns=['Loop MT2'])
        # wt_df_rmsf['Loop MT1'] = mt1_df_rmsf[['Loop MT1']]
        wt_df_rmsf['Loop MT2'] = mt2_df_rmsf[['Loop MT2']]
        wt_df_rmsf['Index']=wt_df_rmsf.index
        wt_df_rmsf.to_csv('loop_rmsf_data_wtmt2.csv')

        print("start plot")
        # And plot the RMSF for each selection
        # f = plt.figure(1)
        ax = wt_df_rmsf.plot(x='Index',y=['Loop WT','Loop MT2'])
        ax.set_ylabel(r'RMSF ($\AA$)')
        ax.set_title('RMSF Loop')
        # ax.axvspan(140, 146, zorder=0, alpha=0.2, color='orange', label='Loop')
        ax.set_xlabel('Residue number')
        # plt.xlabel('Residue number')
        # plt.ylabel('RMSF ($\AA$)')
        # plt.title('C_alphas')
        # plt.axvspan(140, 146, zorder=0, alpha=0.2, color='orange', label='Loop')
        # plt.legend();
        plt.savefig('rmsf_loop_wtmt2.png')


if __name__ == "__main__":
    pdb_list = "pdb_list.log"
    analysis = Analysis()
    # analysis.detect_three_residue_tight_loop_serial(pdb_list)
    analysis.detect_three_residue_tight_loop_MPI(pdb_list)
    # analysis.filtering_redundant_results()
    # analysis.replace_pdb_text()
    # analysis.download_pdb_structures()
    # analysis.filter_surface_results()
    # analysis.compare_phi_difference()
    # analysis.pca_analysis_first_part()
    # analysis.pca_analysis_second_part()
    # analysis.calculate_rmsd()
    # analysis.calculate_rmsf()
    # analysis.project_mt_pc_onto_wt()
    # analysis.pc1_vs_pc2_wt()
    # analysis.pc1_vs_pc2_mt1()
    # analysis.pc1_vs_pc2_mt2()