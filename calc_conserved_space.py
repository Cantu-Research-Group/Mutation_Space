#!/usr/bin/env python
import os
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description='Computes spatially correlated residues based. Spatial alignment of the proteins is done based upon positioning of most conserved distant residues reported in sequence-aligned file, though this can be disabled to use user-aligned structures.')
parser.add_argument('listfile', type=str, help='Name of a file containing the list of PDB (path/)files to be processed, with (path/)names on newlines')
parser.add_argument('msa', type=str, help='Name of the Multiple Sequence Alignment file in ClustalW format')
parser.add_argument('-noalign', action='store_true', help='If flag is given, the program will not align the vectors according to most distant conserved residue positioning')
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


if __name__ == '__main__':
    
    #Retrieve conserved sequence positions from the MSA
    conserved_res_indices = retrieve_cons_res(args.msa)

    #Check presence of PDB files in files and MSA
    pdb_files = [x.strip() for x in open(args.listfile, 'r').readlines()]
    pdb_files = check_files_present(pdb_files, conserved_res_indices.keys())

    #For each structure

        #Extract CA positions
            #Average coordinates of alternate positions if needed
        #Identify conserved residues
            #Calculate the ccenter of mass of all conserved residues
            #Calculate the most distant residue indices
        
        #If doing conserved alignment
            #Alter the coordinates
                #Translate to COM_conserved_res as (0,0,0)
                #Rotate to align vectors in planes
                #Save the altered coordinates

        #Calculate the vectors for each residue relative to Dist_r1, Dist_r2, and COM_conserved_res
            #Store vector results





