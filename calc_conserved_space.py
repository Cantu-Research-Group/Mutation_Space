#!/usr/bin/env python
import os, sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description='Computes spatially correlated residues based. Spatial alignment of the proteins is done based upon positioning of most conserved distant residues reported in sequence-aligned file, though this can be disabled to use user-aligned structures.')
parser.add_argument('listfile', type=str, help='Name of a file containing the list of PDB (path/)files to be processed, with (path/)names on newlines')
parser.add_argument('msa', type=str, help='Name of the Multiple Sequence Alignment file in ClustalW format')
parser.add_argument('-noalign', action='store_true', help='If flag is given, the program will not align the vectors according to most distant conserved residue positioning')
parser.add_argument('-store', type=str, default='processed_pdbs', help='Name of the directory to generate to store the processed PDBs (default: processed_pdbs)')
args = parser.parse_args()

def retrieve_cons_res(fname): 
    msa = open(fname, 'r').readlines()
    seqlength = 0
    conserved = str()
    aligned_seq = defaultdict(str)

    #Read in the sequences and the conserved residue info
    for line in msa[1:]:
        if line == "\n": continue 
        elif line[0] == " ": 
            conserved += line[-(seqlength+1):-1]
        else:
            seqlength = len(line.split()[1].strip())
            aligned_seq[line.split()[0]] += line.split()[1].strip() 

    #Convert aligned residue indices to sequence indices for each structure
    conserved_indices = [x for x, c in enumerate(conserved) if c == "*"]
    conserved_res_indices = defaultdict(list)
    for seq in aligned_seq.keys():
        for c_index in conserved_indices:
            seq_pos = c_index - aligned_seq[seq][:c_index].count("-")
            conserved_res_indices[seq].append(seq_pos)
    return conserved_res_indices 

def check_files_present(pdbs_list, names_msa):
    for name in pdbs_list:
        if Path(name).is_file() == False: #Confirm PDB files present
            print("WARNING: Could not locate file {0}".format(name))
            pdbs_list.remove(name)
        elif (name not in names_msa) and (name.split("/")[-1] not in names_msa): #Confirm PDBs have corresponding MSA info present
            print("WARNING: Could not identify sequence alignment for file {0}".format(name))
            pdbs_list.remove(name)
    return pdbs_list


def extract_and_save_CAs(file_name, save_dir):
    #Extract coordinate info
    coords_lines = defaultdict(list)
    with open(file_name, 'r') as readpdb:
        for line in readpdb:
            if line.startswith('TER'):    #Terminate at first TER sign
                break
            elif (line.startswith('ATOM') and line[12:16].strip() == "CA") or (line.startswith('HETATM') and line[12:16].strip() == "CA"):
                coords_lines[int(line[22:26])].append(line)    #Save the CA lines according to ResNum
    #Write coords to PDB
    xyz_coords = []
    with open(save_dir+"/"+file_name.rsplit("/",1)[-1], 'w') as savepdb:
        for resnum in coords_lines:
            if len(coords_lines[resnum]) != 1:   #Average alternate conformation atoms based upon Occupancy
                header = coords_lines[resnum][0][:30]
                footer = coords_lines[resnum][0][54:]
                avg_x = 0.0
                avg_y = 0.0
                avg_z = 0.0
                for alt_atom in coords_lines[resnum]:
                    avg_x += float(alt_atom[30:38])*float(alt_atom[54:60])
                    avg_y += float(alt_atom[38:46])*float(alt_atom[54:60]) 
                    avg_z += float(alt_atom[46:54])*float(alt_atom[54:60])
                savepdb.write(header+"{:>8}".format(round(avg_x,3))+"{:>8}".format(round(avg_y,3))+"{:>8}".format(round(avg_z,3))+footer)
                xyz_coords.append([resnum, avg_x, avg_y, avg_z])

            else:
                savepdb.write(coords_lines[resnum][0])
                xyz_coords.append([resnum, float(coords_lines[resnum][0][30:38]), float(coords_lines[resnum][0][38:46]), float(coords_lines[resnum][0][46:54])])

    xyz_coords = pd.DataFrame(xyz_coords, columns=['ResNum','X','Y','Z'])
    return xyz_coords

def calculate_distant_residues(df):
    dist = 0
    res1 = 0
    res2 = 0
    indices = df.index
    for i in range(0, len(indices)):
        for j in range(i, len(indices)):
            if i == j: continue
            x1 = df.iloc[i]['X']
            test_dist = ( (df.iloc[i]['X']-df.iloc[j]['X'])**2 + (df.iloc[i]['Y']-df.iloc[j]['Y'])**2 + (df.iloc[i]['Z']-df.iloc[j]['Z'])**2 )**0.5
            if test_dist > dist:
                dist, res1, res2 = test_dist, i, j
    return df.iloc[res1]['ResNum'], df.iloc[res2]['ResNum']

if __name__ == '__main__':
    
    #Retrieve conserved sequence positions from the MSA
    conserved_res_indices = retrieve_cons_res(args.msa)

    #Check presence of PDB files in files and MSA
    pdb_files = [x.strip() for x in open(args.listfile, 'r').readlines()]
    pdb_files = check_files_present(pdb_files, conserved_res_indices.keys())

    #Generate the PDB storage directory if it is not present:
    if Path(args.store).is_dir() == False:
        os.mkdir(args.store)

    #For each structure
    for pdb_file in pdb_files:
        
        #Extract CA positions data, and save CAs to storage directory
        CA_data = extract_and_save_CAs(pdb_file, args.store)

        #Identify conserved residues
        CA_data_conserved = CA_data[np.isin(CA_data['ResNum'],conserved_res_indices[pdb_file]) == True]

        #Calculate the averaged XYZ center of mass/geometry of all conserved residues
        COM_x = CA_data_conserved["X"].mean()
        COM_y = CA_data_conserved["Y"].mean()
        COM_z = CA_data_conserved["Z"].mean()

        #Calculate the most distant residue indices
        distant_res_i, distant_res_j = calculate_distant_residues(CA_data_conserved)
        
        #If doing conserved alignment
        if args.noalign == False:
            #Translate all coordinates so the center of conserved residues (COM) is at (0,0,0)
            CA_data['X_align'] = CA_data['X']-COM_x
            CA_data['Y_align'] = CA_data['Y']-COM_y
            CA_data['Z_align'] = CA_data['Z']-COM_z

            #Rotate to align vectors in planes so that r1 = <0, 0, d_r1_cr_com>

            #Save the altered coordinates

        #Calculate the vectors for each residue relative to Dist_r1, Dist_r2, and COM_conserved_res
            #Store vector results





