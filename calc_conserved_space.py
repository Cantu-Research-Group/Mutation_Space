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
parser.add_argument('ref', type=str, help='Name of the PDB file to serve as a reference for the vector analysis')
parser.add_argument('-noalign', action='store_true', help='If flag is given, the program will not align the vectors according to most distant conserved residue positioning')
parser.add_argument('-store', type=str, default='processed_pdbs', help='Name of the directory to generate to store the processed PDBs (default: processed_pdbs)')
args = parser.parse_args()

def retrieve_cons_indices(fname): 
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
    xyz_coords = defaultdict(list)
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
                xyz_coords[len(xyz_coords.keys())] = [resnum, [avg_x, avg_y, avg_z], header, footer]

            else:
                savepdb.write(coords_lines[resnum][0])
                xyz_coords[len(xyz_coords.keys())] = [resnum, [float(coords_lines[resnum][0][30:38]), float(coords_lines[resnum][0][38:46]), float(coords_lines[resnum][0][46:54])], coords_lines[resnum][0][:30], coords_lines[resnum][0][54:]]
    return xyz_coords #Format { Index: [resnum, [X, Y, Z], header, footer], ... }

def calculate_COM_and_distances(df, indices):
    COM_xyz = [0.0, 0.0, 0.0]
    dist = 0.0
    res1 = 0
    res2 = 0
    for i in indices:
        COM_xyz = [x+y for x,y in zip(COM_xyz, df[i][1])]#Sum XYZs
        for j in indices[indices.index(i):]:
            if i == j: continue
            test_dist = ( (df[i][1][0]-df[j][1][0])**2 + (df[i][1][1]-df[j][1][1])**2 + (df[i][1][2]-df[j][1][2])**2)**0.5
            if test_dist > dist:
                dist, res1, res2 = test_dist, i, j #Save most distant residue pairs
    COM_xyz = [x/len(indices) for x in COM_xyz] #Average XYZs
    return COM_xyz, res1, res2

def translate_coords(df, trans):
    for element in df:
        translated = [x-y for x,y in zip(df[element][1], trans)]
        df[element][1] = translated
    return df

def calc_rotation_matrix(axis, df, dist_res_1, dist_res_2):
    def matrix_rot(ax, t):
        if ax == 'X-axis':
            return np.matrix([[1,0,0], [0,np.cos(t),-np.sin(t)], [0,np.sin(t),np.cos(t)]])
        elif ax == 'Y-axis':
            return np.matrix([[np.cos(t),0,np.sin(t)], [0,1,0], [-np.sin(t),0,np.cos(t)]])
        elif ax == 'Z-axis':
            return np.matrix([[np.cos(t),-np.sin(t),0], [np.sin(t),np.cos(t),0], [0,0,1]])

    if axis == 'X-axis':
        p1 = df[dist_res_1][1][1] #Y
        p2 = df[dist_res_1][1][2] #Z
    elif axis == 'Y-axis':
        p1 = df[dist_res_1][1][0] #X
        p2 = df[dist_res_1][1][2] #Z
    elif axis == 'Z-axis':
        p1 = df[dist_res_2][1][1] #Y
        p2 = df[dist_res_2][1][0] #X
    mag_p = (p1**2 + p2**2)**0.5
    theta = np.arccos(p2 / mag_p)
    
    if axis == 'X-axis':
        rot_mat = matrix_rot('X-axis',theta)
        test_mat = np.matrix([[df[dist_res_1][1][0]],[df[dist_res_1][1][1]],[df[dist_res_1][1][2]]])
        test_rot = np.dot(rot_mat, test_mat)
        if -0.01 <= test_rot[1] <= 0.01:
            return rot_mat
        else:
            rot_mat = matrix_rot('X-axis',-theta)
            return rot_mat
    elif axis == 'Y-axis':
        rot_mat = matrix_rot('Y-axis',theta)
        test_mat = np.matrix([[df[dist_res_1][1][0]],[df[dist_res_1][1][1]],[df[dist_res_1][1][2]]])
        test_rot = np.dot(rot_mat, test_mat)
        if -0.01 <= test_rot[0] <= 0.01:
            return rot_mat
        else:
            rot_mat = matrix_rot('Y-axis',-theta)
            return rot_mat
    elif axis == 'Z-axis':
        rot_mat = matrix_rot('Z-axis',theta)
        test_mat = np.matrix([[df[dist_res_2][1][0]],[df[dist_res_2][1][1]],[df[dist_res_2][1][2]]])
        test_rot = np.dot(rot_mat, test_mat)
        if -0.01 <= test_rot[1] <= 0.01:
            return rot_mat
        else:
            rot_mat = matrix_rot('Z-axis',-theta)
            return rot_mat

def rotate_coords(df, mat):
    for res in df:
        start_mat = np.matrix([[df[res][1][0]],[df[res][1][1]],[df[res][1][2]]])
        rotated_mat = np.dot(mat, start_mat).tolist()
        df[res][1] = [rotated_mat[0][0], rotated_mat[1][0], rotated_mat[2][0]] #Save results to XYZ positions
    return df

def calc_vectors_and_cos_sim(df, index_i, index_j, com, ref):

    def vector_magnitude(xyz):
        return ((xyz[0]**2)+(xyz[1]**2)+(xyz[2]**2))**0.5

    def cos_similarity(v1_xyz, v2_xyz, v1_mag, v2_mag): #1 is most similar, 0 is least similar
        if v1_mag == 0 or v2_mag == 0: #NEED TO FIGURE ALTERNATIVE SIMILARITY MEASURE SINCE MAGS FOR SOME VECTORS ARE 0
            return 0.0 
        else:
            return np.dot(v1_xyz, v2_xyz) / (v1_mag*v2_mag)

    res_in_threshold = defaultdict(list)
    for res in df:
        v_1 = [x-y for x,y in zip(df[res][1], df[index_i][1])] #Difference between XYZs of each residue to conserved distant res1
        v_2 = [x-y for x,y in zip(df[res][1], df[index_j][1])] #Difference between XYZs of each residue to conserved distant res2
        v_com = [x-y for x,y in zip(df[res][1], com)] #Difference between XYZs of each residue to COM
        v_1_mag = vector_magnitude(v_1)
        v_2_mag = vector_magnitude(v_2)
        v_com_mag = vector_magnitude(v_com)

        v_1.append(v_1_mag)
        v_2.append(v_2_mag)
        v_com.append(v_com_mag)
        df[res][2] = v_1
        df[res][3] = v_2
        df[res].append(v_com)

        if ref == None: continue #If reference structure not provided, do not proceed

        for ref_res in ref:
            r = ref[ref_res]
            #Calculate average Vector Magnitude Differences to each reference residue
            vec_mag_diff = (((r[2][3]-v_1_mag)**2 + (r[3][3]-v_2_mag)**2 + (r[4][3]-v_com_mag)**2 ) / 3)**0.5
            #Calculate average Vector Cosine Similarity to each reference residue
            v_1_cos_sim = cos_similarity(v_1[:3], r[2][:3], v_1_mag, r[2][3])
            v_2_cos_sim = cos_similarity(v_2[:3], r[3][:3], v_2_mag, r[3][3])
            v_com_cos_sim = cos_similarity(v_com[:3], r[4][:3], v_com_mag, r[4][3])
            avg_cos_sim = (v_1_cos_sim + v_2_cos_sim + v_com_cos_sim)/3

            #If similarities between residues are within threshold, store the value
            if vec_mag_diff < 2.6 * avg_cos_sim - 1.35:
                res_in_threshold[ref_res].append([res, vec_mag_diff])

    if ref == None:
        return df
    else:
        return df, res_in_threshold

def remove_single_pairs(fin, init):
    res_of_singles = defaultdict(list)
    for ref_res in init: #Gather info of subject residues for all single-paired ref residues
        if len(init[ref_res]) == 1:
            res_of_singles[init[ref_res][0][0]].append(ref_res)

    for sub_res in res_of_singles:
        if len(res_of_singles[sub_res]) == 1: #These residues have a 1-ref-res to 1-sub-res correlation and can be directly aligned
            fin[res_of_singles[sub_res][0]] = sub_res
            del init[res_of_singles[sub_res][0]]
        else: #These residues have a multi-ref-res to 1-sub-res correlation. Smallest vector mag will be used to choose alignment
            best_res = res_of_singles[sub_res][0] 
            best_mag = init[best_res][0][1]
            for ref_res in res_of_singles[sub_res]:
                if init[ref_res][0][1] < best_mag:
                    best_res = ref_res
                    best_mag = init[ref_res][0][1]
            fin[best_res] = sub_res
            for ref_res in res_of_singles[sub_res]:
                del init[ref_res]

    for ref_res in init: #Remove used subject residues from available remaining multi-correlated residues
        for pair in init[ref_res]:
            if pair[0] in list(res_of_singles.keys()): 
                init[ref_res].remove(pair)

    return fin, init

def remove_double_pairs(fin, init):
    res_of_doubles = defaultdict(list)
    for ref_res in init: #Gather info of dual subject residues for all double-paired ref residues
        if len(init[ref_res]) == 2:
            res_of_doubles[frozenset([x[0] for x in init[ref_res]])].append(ref_res)

    sub_res_to_remove = []
    for sub_res in res_of_doubles:
        if len(res_of_doubles[sub_res]) == 2: #These are two ref residues pointing to the same two sub residues.
            ref1_res1_id = init[res_of_doubles[sub_res][0]][0][0] #Info on first reference residue's vectors with subjects
            ref1_res1_vec = init[res_of_doubles[sub_res][0]][0][1]
            ref1_res2_id = init[res_of_doubles[sub_res][0]][1][0]
            ref1_res2_vec = init[res_of_doubles[sub_res][0]][1][1]
            ref2_res1_id = init[res_of_doubles[sub_res][1]][0][0]  #Info on second reference residue's vectors with subjects
            ref2_res1_vec = init[res_of_doubles[sub_res][1]][0][1]
            ref2_res2_id = init[res_of_doubles[sub_res][1]][1][0]
            ref2_res2_vec = init[res_of_doubles[sub_res][1]][1][1]

            if ref1_res1_id == ref2_res1_id: #Assign by minimizing net vector magnitudes
                if ref1_res1_vec + ref2_res2_vec < ref1_res2_vec + ref2_res1_vec:
                    fin[res_of_doubles[sub_res][0]] = ref1_res1_id
                    fin[res_of_doubles[sub_res][1]] = ref2_res2_id
                else:
                    fin[res_of_doubles[sub_res][0]] = ref1_res2_id
                    fin[res_of_doubles[sub_res][1]] = ref2_res1_id
            else:
                if ref1_res1_vec + ref2+res1+vec < ref1_res2_vec + ref2_res2_vec:
                    fin[res_of_doubles[sub_res][0]] = ref1_res1_id
                    fin[res_of_doubles[sub_res][1]] = ref2_res1_id
                else:
                    fin[res_of_doubles[sub_res][0]] = ref1_res2_id
                    fin[res_of_doubles[sub_res][1]] = ref2_res2_id
            sub_res_to_remove.extend(list(sub_res)) #Add subject residues to list of residues to remove from alignment options
            del init[res_of_doubles[sub_res][0]] #Remove aligned reference residues from options
            del init[res_of_doubles[sub_res][1]]

    for ref_res in init: #Remove used subject residues from available remaining multi-correlated residues
        for pair in init[ref_res]:
            if pair[0] in sub_res_to_remove:
                init[ref_res].remove(pair)

    return fin, init

def remove_sequential_pairs(fin, init):
    aligned_res = list(fin.keys())
    ref_res = next(iter(init)) 
    i = 1
    a_bool, b_bool = False, False
    a, b = ref_res, ref_res
    while i <= 15: #Find nearest aligned res within up to 15 positions of ref index
        if a_bool == False:
            a = ref_res-i
        if b_bool == False:
            b = ref_res+i
        if a in aligned_res:
            a_bool = True
        if b in aligned_res:
            b_bool = True

        if a_bool == True and b_bool == True:
            break
        else:
            i+=1

    print(a, b, a_bool, b_bool)
    chosen_sub_res = None
    shortest_dist = None
    if a_bool == False and b_bool == True: #Case where only upper bound ref_res was found for comparison
        for sub_res in init[ref_res]: #Identify subject residue (smaller than ref_res subject) and (with shortest sequential dist to ref_res subject)
            if sub_res[0] < fin[b]:
                dist = fin[b]-sub_res[0]
                if chosen_sub_res == None:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
                elif dist < shortest_dist:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
        if chosen_sub_res == None: #If all options are only larger than ref_res subject, identify subject residue with shortest reverse sequential distance to ref_res subject
            for sub_res in init[ref_res]:
                dist = sub_res[0]-fin[b]
                if shortest_dist == None:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
                elif dist < shortest_dist:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
        fin[ref_res] = chosen_sub_res
        del init[ref_res]
        print(fin)
        return fin, init

    elif a_bool == True and b_bool == False: #Case where only the lower bound ref_res was found for comparison
        for sub_res in init[ref_res]: #Identify subject residue (larger than ref_res subject) and (with shortest sequential dist to ref_res subject)
            if sub_res[0] > fin[a]:
                dist = sub_res[0]-fin[a]
                if chosen_sub_res == None:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
                elif dist < shortest_dist:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
        if chosen_sub_res == None: #If all options are only smaller than ref_res subject, identify subject residue with shortest reverse sequential distance to ref_res subject
            for sub_res in init[ref_res]:
                dist = fin[a]-sub_res[0]
                if shortest_dist == None:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
                elif dist < shortest_dist:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
        fin[ref_res] = chosen_sub_res
        del init[ref_res]
        print(fin)
        return fin, init

    elif a_bool == True and b_bool == True:
        for sub_res in init[ref_res]: #Identify subject residue between both ref_res a and ref_res b subjects
            if fin[a] < sub_res[0] < fin[b] or fin[b] < sub_res[0] < fin[a]:
                dist = min([abs(sub_res[0]-fin[a]), abs(sub_res[0]-fin[b])])
                if chosen_sub_res == None:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
                elif dist < shortest_dist:
                    shortest_dist = dist
                    chosen_sub_res = sub_res[0]
        if chosen_sub_res == None: #If all options are not between either ref_res a or ref_res b subjects, choose the one with shortest distance
            for sub_res in init[ref_res]:
                dist = min([abs(sub_res[0]-fin[a]), abs(sub_res[0]-fin[b])]):


    else:
        sys.exit()


    return fin, init




    
    sys.exit()

if __name__ == '__main__':
    
    #Retrieve conserved sequence positions from the MSA
    conserved_res_indices = retrieve_cons_indices(args.msa)
    
    #Check presence of PDB files in files and MSA
    pdb_files = [x.strip() for x in open(args.listfile, 'r').readlines()]
    pdb_files = check_files_present(pdb_files, conserved_res_indices.keys())

    #Generate the PDB storage directory if it is not present:
    if Path(args.store).is_dir() == False:
        os.mkdir(args.store)

    #Extract reference/target structure information
    ref_data = extract_and_save_CAs(args.ref, args.store)
    
    #Calculate the reference/target averaged XYZ center of mass/geometry of all conserved residues & Identify distant conserved res
    ref_COM, ref_distant_res_i, ref_distant_res_j = calculate_COM_and_distances(ref_data, conserved_res_indices[args.ref])

    #If doing conserved alignment
    if args.noalign == False:
        #Translate all coordinates so the center of conserved residues (COM) is at (0,0,0)
        ref_data = translate_coords(ref_data, ref_COM)
        ref_COM = [x-y for x,y in zip(ref_COM, ref_COM)] #Reassign COM to (0,0,0)
        
        #Rotate coords to align vectors in planes
        rotation_matrix_X = calc_rotation_matrix('X-axis', ref_data, ref_distant_res_i, ref_distant_res_j)
        ref_data = rotate_coords(ref_data, rotation_matrix_X)
        rotation_matrix_Y = calc_rotation_matrix('Y-axis', ref_data, ref_distant_res_i, ref_distant_res_j)
        ref_data = rotate_coords(ref_data, rotation_matrix_Y)
        rotation_matrix_Z = calc_rotation_matrix('Z-axis', ref_data, ref_distant_res_i, ref_distant_res_j)
        ref_data = rotate_coords(ref_data, rotation_matrix_Z)
        
        #Save the altered coordinates
        with open(args.store+"/aligned-"+args.ref.rsplit("/",1)[-1], 'w') as savefile:
            for res in ref_data:
                #Save header + rounded X + rounded Y + rounded Z + footer
                savefile.write(ref_data[res][2]+"{:>8}".format(round(ref_data[res][1][0],3))+"{:>8}".format(round(ref_data[res][1][1],3))+"{:>8}".format(round(ref_data[res][1][2],3))+ref_data[res][3])
    
    #Calculate reference vectors and vector properties - replace header&footer info
    #New dict format: { Index: [resnum, [X, Y, Z], [v1_x, v1_y, v1_z, v1_mag], [v2_x, v2_y, v2_z, v2_mag], [com_x, com_y, com_z, com_mag] ] ...}
    ref_data = calc_vectors_and_cos_sim(ref_data, ref_distant_res_i, ref_distant_res_j, ref_COM, None)


    #Extract info and calculate properties for each non-reference subject structure
    for pdb_file in [x for x in pdb_files if x != args.ref]:
        
        #Extract subject CA positions data, and save CAs to storage directory
        sub_data = extract_and_save_CAs(pdb_file, args.store)
       
        #Extract subject averaged XYZ center of mass/geometry of all conserved residues & Identify distant conserved res
        sub_COM, sub_distant_res_i, sub_distant_res_j = calculate_COM_and_distances(sub_data, conserved_res_indices[pdb_file])

        #If doing conserved alignment
        if args.noalign == False:
            #Translate all subject coordinates so the center of conserved residues (COM) is at (0,0,0)
            sub_data = translate_coords(sub_data, sub_COM)
            sub_COM = [x-y for x,y in zip(sub_COM, sub_COM)]

            #Rotate to align vectors in planes
            rotation_matrix_X = calc_rotation_matrix('X-axis', sub_data, sub_distant_res_i, sub_distant_res_j)
            sub_data = rotate_coords(sub_data, rotation_matrix_X)
            rotation_matrix_Y = calc_rotation_matrix('Y-axis', sub_data, sub_distant_res_i, sub_distant_res_j)
            sub_data = rotate_coords(sub_data, rotation_matrix_Y)
            rotation_matrix_Z = calc_rotation_matrix('Z-axis', sub_data, sub_distant_res_i, sub_distant_res_j)
            sub_data = rotate_coords(sub_data, rotation_matrix_Z)

            #Save the altered coordinates
            with open(args.store+"/aligned-"+args.ref.rsplit("/",1)[-1], 'w') as savefile:
                for res in sub_data:
                    savefile.write(sub_data[res][2]+"{:>8}".format(round(sub_data[res][1][0],3))+"{:>8}".format(round(sub_data[res][1][1],3))+"{:>8}".format(round(sub_data[res][1][2],3))+sub_data[res][3])
        
        #Calculate subject vectors and vector properties - replace header&footer info
        #Also calculate vector magnitude differences and cosine similarities to reference structure
        #Store residues within similarity threshold. Format: {ref_res_index: [ [sub_res_index, vector_mag_diff] ...] }
        sub_data, conserved_res = calc_vectors_and_cos_sim(sub_data, sub_distant_res_i, sub_distant_res_j, sub_COM, ref_data)

        final_alignment = defaultdict(dict)
        for item in conserved_res:
            print(item, conserved_res[item])
            print("*****")
        while True: #Remove Single and Double-correlation pairs
            unchanged_single_double = len(final_alignment)
            
            #Align single-correlation pairs
            while True:
                unchanged_single = len(final_alignment)
                final_alignment, conserved_res = remove_single_pairs(final_alignment, conserved_res)
                for item in conserved_res:
                    print(item, conserved_res[item])
                print("~~~~~~~~~~~~~~~~~~~~~")
                if len(final_alignment) == unchanged_single:
                    break
            
            #Align double-correlation pairs
            final_alignment, conserved_res = remove_double_pairs(final_alignment, conserved_res)
            for item in conserved_res:
                print(item, conserved_res[item])
            print("*********************************************************************************")
            if len(final_alignment) == unchanged_single_double:
                break
        print(pdb_file)
        for item in conserved_res:
            print(item, conserved_res[item])
        print(final_alignment)
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        #All Single and Double-correlation pairs are exhausted. Align by sequence relative to aligned residues
        while True:
            unchanged = len(conserved_res)
            for item in conserved_res:
                print(item, conserved_res[item], "\t", ref_data[item][0], sub_data[conserved_res[item][0][0]][0], sub_data[conserved_res[item][1][0]][0])
            final_alignment, conserved_res = remove_sequential_pairs(final_alignment, conserved_res)

            if len(conserved_res) == unchanged:
                print(conserved_res)
                break


