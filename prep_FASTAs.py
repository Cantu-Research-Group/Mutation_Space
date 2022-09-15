#!/usr/bin/env python
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Extracts the CA atoms for a user-provided list of PDBs and generates a corresponding FASTA formatted file')
parser.add_argument('listfile', type=str, help='Name of a file containing the list of PDB (path/)files to be processed, with (path/)names on newlines. A specific chain for a structure can be chosen by listing the chain letter after the file name, otherwise the default "A" chain will be used.')
parser.add_argument('-o', type=str, default='gen_FASTAs.dat', help='Name of the FASTA outputfile (default gen_FASTAs.dat)')
parser.add_argument('-a', nargs = '+', help='Additional non-canonical amino acid 3-letter to 1-letter code(s). Please provide in the format XXX:X with multiple entries space-separated (e.g. XXX:X YYY:Y ...')
args = parser.parse_args()

aa_three_to_one = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
        'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L',
        'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
        'TRP':'W', 'TYR':'Y', 'VAL':'V'}

def check_files_present(names):
    for name in names:
        check_name = name.split()[0]
        if Path(check_name).is_file() == False:
            print("WARNING: Could not find file {0}".format(check_name))
            names.remove(name)
    return names
        

if __name__ == '__main__':

    #Confirm all user-listed PDBs are present
    fnames = [x.strip() for x in open(args.listfile, 'r').readlines()]
    fnames = check_files_present(fnames)

    #Update the FASTA dictionary if needed:
    if args.a != None:
        for item in args.a:
            aa_three_to_one[item.split(":")[0]] = item.split(":")[1]

    #Extract chain and FASTA info for each PDB
    fastas = {}
    for name in fnames:
        resnum = []
        fasta = str()
        if len(name.split()) == 1:
            chain = 'A'
        else:
            chain = name.split()[-1]
        with open(name.split()[0], 'r') as readpdb:
            for line in readpdb:
                if (line.startswith('ATOM') and line[12:16].strip() == "CA" and line[21] == chain and int(line[22:26]) not in resnum) or \
                        (line.startswith('HETATM') and line[12:16].strip() == "CA" and line[21] == chain and int(line[22:26]) not in resnum):
                    try:
                        letter = aa_three_to_one[line[17:20]] #Store the FASTA sequence
                        fasta += letter
                        resnum.append(int(line[22:26]))
                    except:
                        fasta += "X"
                        resnum.append(int(line[22:26]))
                        print("WARNING: Amino acid {0} in file {1} not found in amino acid key, assigning as X in FASTA sequence".format(line[17:20],name))
        fastas[name] = fasta

    #Write the FASTA file
    with open(args.o, 'w') as savefile:
        for name in fastas:
            pdb_file = name.split()[0].split("/")[-1] #Remove any file location preceeding information
            savefile.write(">"+pdb_file+"\n")
            seq = [fastas[name][i:i+80] for i in range(0,len(fastas[name]), 80)]
            for i in seq:
                savefile.write(i+"\n")


