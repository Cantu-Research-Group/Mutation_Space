#!/usr/bin/env python
import os, sys
import argparse, math, scoring_matrix
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description='Identifies spatially correlated residues and calculates their MSSC scores. Spatial alignment of the proteins is done based upon positioning of most distant conserved residues reported in sequence-aligned file, though this can be disabled to use user-aligned structures.')
parser.add_argument('listfile', type=str, help='Name of a file containing the list of PDB (path/)files to be processed, with (path/)names on newlines')
parser.add_argument('msa', type=str, help='Name of the Multiple Sequence Alignment file in ClustalW format')
parser.add_argument('ref', nargs = '+', help='(Path/)Name of the PDB file to serve as a reference/target for the vector analysis. If using non-"A" protein chain, indicate chain letter after the file name (e.g. A123.pdb B)')
parser.add_argument('-save', type=str, default='results.xlsx', help='Name of .xlsx results file to generate, default results.xlsx')
parser.add_argument('-noalign', action='store_true', help='If flag is given, the program will not align the proteins according to most distant conserved residue positioning')
parser.add_argument('-store', type=str, default='processed_pdbs', help='Name of the directory to generate and store the processed PDBs (default: processed_pdbs)')
parser.add_argument('-criteria', type=float, default=0.1, help='Threshold for TS-SS similarity, default=0.1')
parser.add_argument('-tightcriteria', type=float, default=0.03, help='Threshold for TS-SS similarity tight convergence where atoms with a TS-SS score less than the value are automatically considered to be correlated, default=0.03')
parser.add_argument('-pydist', action='store_const', const='pymol_distances.log', help='If flag is given, the program will output a file to store pymol commands useful for visualizing interconnections between correlated atoms. Alternate name for output file can be specified after flag (default: pymol_distances.log)')
parser.add_argument('-pyheat', action='store_const', const='pymol_heatmap.log', help='If flag is given, the program will output a file to store pymol commands useful for visualizing MSSC score heatmap distributions. Alternate name for output file can be specified after flag (default: pymol_heatmap.log)')
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
        check_name = name.split()[0]
        if Path(check_name).is_file() == False: #Confirm PDB files present
            print("WARNING: Could not locate file {0}".format(check_name))
            pdbs_list.remove(name)
        elif (check_name not in names_msa) and (check_name.split("/")[-1] not in names_msa): #Confirm PDBs have corresponding MSA info present
            print("WARNING: Could not identify sequence alignment for file {0}".format(check_name))
            pdbs_list.remove(name)
    return pdbs_list


def extract_and_save_CAs(file_name, save_dir):
    #Extract coordinate info
    coords_lines = defaultdict(list)
    if len(file_name) == 1:
        chain = 'A'
    else:
        chain = file_name[1]
    with open(file_name[0], 'r') as readpdb:
        for line in readpdb:
            if (line.startswith('ATOM') and line[12:16].strip() == "CA" and line[21] == chain) or (line.startswith('HETATM') and line[12:16].strip() == "CA" and line[21] == chain):
                coords_lines[int(line[22:26])].append(line)    #Save the CA lines according to ResNum
    #Write coords to PDB
    xyz_coords = defaultdict(list)
    with open(save_dir+"/"+file_name[0].rsplit("/",1)[-1], 'w') as savepdb:
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

def calc_vectors_sim(df, index_i, index_j, com, ref, threshold):

    def vector_magnitude(xyz):
        return math.sqrt(sum(math.pow(i,2) for i in xyz))

    def calc_TS_SS(v1_xyz, v2_xyz): #TS-SS algorithm referenced from DOI: 10.1109/BigDataService.2016.14
        def calc_vec_cosine(v1, v2) :
            return sum(i*j for i,j in zip(v1,v2)) / (vector_magnitude(v1) * vector_magnitude(v2))

        def calc_theta(v1, v2) :
            return math.acos(round(calc_vec_cosine(v1,v2),10)) + math.radians(10)

        def calc_triangle(v1, v2) :
            angle = math.radians(calc_theta(v1,v2))
            return (vector_magnitude(v1) * vector_magnitude(v2) * math.sin(angle)) / 2

        def calc_sector(v1, v2) :
            euclid_dist = math.sqrt(sum(math.pow((i-j),2) for i,j in zip(v1, v2)))
            mag_diff = abs(vector_magnitude(v1)-vector_magnitude(v2))
            angle = calc_theta(v1, v2)
            return math.pi * math.pow((euclid_dist+mag_diff),2) * angle/360

        return calc_triangle(v1_xyz, v2_xyz) * calc_sector(v1_xyz, v2_xyz)

    res_in_threshold = defaultdict(list)
    for res in df:
        resname = df[res][2][17:20]

        v_1 = [(x-y+0.001) for x,y in zip(df[res][1], df[index_i][1])] #Difference between XYZs of each residue to conserved distant res1
        v_2 = [(x-y+0.001) for x,y in zip(df[res][1], df[index_j][1])] #Difference between XYZs of each residue to conserved distant res2
        v_com = [(x-y+0.001) for x,y in zip(df[res][1], com)] #Difference between XYZs of each residue to COM
        v_1_mag = vector_magnitude(v_1)
        v_2_mag = vector_magnitude(v_2)
        v_com_mag = vector_magnitude(v_com)

        v_1.append(v_1_mag)
        v_2.append(v_2_mag)
        v_com.append(v_com_mag)
        df[res][2] = v_1
        df[res][3] = v_2
        df[res].append(v_com)
        df[res].append(resname)

        if ref == None: continue #If reference structure not provided, do not proceed

        for ref_res in ref:
            r = ref[ref_res]

            #If TS-SS similarities between residues are within TS-SS threshold, store the value
            compare_v1 = calc_TS_SS(v_1[:3], r[2][:3])
            compare_v2 = calc_TS_SS(v_2[:3], r[3][:3])
            compare_com = calc_TS_SS(v_com[:3], r[4][:3])
            if (compare_v1+compare_v2+compare_com)/3 < threshold:
                res_in_threshold[ref_res].append([res, (compare_v1+compare_v2+compare_com)/3, resname])

    if ref == None:
        return df
    else:
        return df, res_in_threshold

def remove_tight_pairs(fin, init, threshold):
    res_in_threshold = defaultdict(list)
    for ref_res in init: #Gather info of subject residues with tightest correlations within threshold for each ref_res
        tightest_res = None
        tightest_val = None
        for sub_res in init[ref_res]:
            if sub_res[1] <= threshold:
                if tightest_res == None:
                    tightest_res = sub_res[0]
                    tightest_val = sub_res[1]
                elif sub_res[1] <= threshold and sub_res[1] <= tightest_val:
                    tightest_res = sub_res[0]
                    tightest_val = sub_res[1]
        if tightest_res != None:
            res_in_threshold[tightest_res].append(ref_res)

    for sub_res in list(res_in_threshold):
        if len(res_in_threshold[sub_res])==1: #If subject residue has within-threshold correlation with only 1 ref residue, add to alignment
            fin[res_in_threshold[sub_res][0]] = sub_res
            del init[res_in_threshold[sub_res][0]]
        else:
            del res_in_threshold[sub_res]
    for ref_res in init: #Remove used subject residues from available remaining multi-correlated residues
        init[ref_res] = [x for x in init[ref_res] if x[0] not in list(res_in_threshold)]

    return fin, init


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
    for ref_res in list(init): #Remove used subject residues from available remaining multi-correlated residues
        if not init[ref_res]: #Remove any empty elements
            del init[ref_res]
            continue
        init[ref_res] = [x for x in init[ref_res] if x[0] not in list(res_of_singles)]

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
        if not init[ref_res]: #Remove any empty elements
            del init[ref_res]
            continue
        init[ref_res] = [x for x in init[ref_res] if x[0] not in sub_res_to_remove]

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
        if chosen_sub_res == None: #If all options are only larger than ref_res subject, identify subject residue with shortest vector factor
            for sub_res in init[ref_res]:
                if shortest_dist == None:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]
                elif sub_res[1] < shortest_dist:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]

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
        if chosen_sub_res == None: #If all options are only smaller than ref_res subject, choose subject residue with shortest vector factor
            for sub_res in init[ref_res]:
                if shortest_dist == None:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]
                elif sub_res[1] < shortest_dist:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]

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
        if chosen_sub_res == None: #If all options are not between either ref_res a or ref_res b subjects, choose the one with smallest vector factor
            for sub_res in init[ref_res]:
                if chosen_sub_res == None:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]
                elif sub_res[1] < shortest_dist:
                    shortest_dist = sub_res[1]
                    chosen_sub_res = sub_res[0]
        
    else: #If there are no other aligned res within 15 of the residue, choose the smallest vector mag diff
        for sub_res in init[ref_res]:
            if chosen_sub_res == None:
                shortest_dist = sub_res[1]
                chosen_sub_res = sub_res[0]
            elif sub_res[1] < shortest_dist:
                shortest_dist = sub_res[1]
                chosen_sub_res = sub_res[0]

    fin[ref_res] = chosen_sub_res
    del init[ref_res]
    for res in init:
        init[res] = [x for x in init[res] if x[0] != chosen_sub_res]
        if not init[res]:
            del init[res]

    return fin, init

def calc_score(row, dist):
    orig = row[0]
    score = 0
    for col in row[1:]:
        if pd.isnull(col):
            continue #"Adding" 0 to score for missing residues
        else:
            score += dist[orig][col]
    return score

if __name__ == '__main__':
    
    #Retrieve conserved sequence positions from the MSA
    conserved_res_indices = retrieve_cons_indices(args.msa)
    
    #Check presence of PDB files in files and MSA
    pdb_files = [x.strip() for x in open(args.listfile, 'r').readlines()]
    pdb_files = check_files_present(pdb_files, conserved_res_indices.keys())

    #Generate the PDB storage directory if it is not present:
    if Path(args.store).is_dir() == False:
        os.mkdir(args.store)
    #Remove previous pymol logfiles if going to generate them
    if args.pydist:
        try: os.remove(args.pydist)
        except OSError: pass

    #Extract reference/target structure information
    ref_data = extract_and_save_CAs(args.ref, args.store)
    
    #Calculate the reference/target averaged XYZ center of mass/geometry of all conserved residues & Identify distant conserved res
    ref_COM, ref_distant_res_i, ref_distant_res_j = calculate_COM_and_distances(ref_data, conserved_res_indices[args.ref[0].split("/")[-1]])

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
        with open(args.store+"/aligned-"+args.ref[0].rsplit("/",1)[-1], 'w') as savefile:
            for res in ref_data:
                #Save header + rounded X + rounded Y + rounded Z + footer
                savefile.write(ref_data[res][2]+"{:>8}".format(round(ref_data[res][1][0],3))+"{:>8}".format(round(ref_data[res][1][1],3))+"{:>8}".format(round(ref_data[res][1][2],3))+ref_data[res][3])
    
    #Calculate reference vectors and vector properties - replace header&footer info
    #New dict format: { Index: [resnum, [X, Y, Z], [v1_x, v1_y, v1_z, v1_mag], [v2_x, v2_y, v2_z, v2_mag], [com_x, com_y, com_z, com_mag], ResName] ...}
    ref_data = calc_vectors_sim(ref_data, ref_distant_res_i, ref_distant_res_j, ref_COM, None, args.criteria)

    #Initialize alignment storage dictionary
    results_dataframe = defaultdict(dict)
    results_dataframe_resnum = defaultdict(dict)
    for pos in ref_data:
        results_dataframe[args.ref[0]][ref_data[pos][0]] = ref_data[pos][5]
        results_dataframe_resnum[args.ref[0]][ref_data[pos][0]] = ref_data[pos][0]

    #Extract info and calculate properties for each non-reference subject structure
    for pdb_file in [x.split() for x in pdb_files]:
        if pdb_file[0] == args.ref[0]: continue
        print("Processing structure: {}".format(pdb_file[0].split("/")[-1]))
        
        #Extract subject CA positions data, and save CAs to storage directory
        sub_data = extract_and_save_CAs(pdb_file, args.store)
        
        #Extract subject averaged XYZ center of mass/geometry of all conserved residues & Identify distant conserved res
        sub_COM, sub_distant_res_i, sub_distant_res_j = calculate_COM_and_distances(sub_data, conserved_res_indices[pdb_file[0].split("/")[-1]])

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
            with open(args.store+"/aligned-"+pdb_file[0].rsplit("/",1)[-1], 'w') as savefile:
                for res in sub_data:
                    savefile.write(sub_data[res][2]+"{:>8}".format(round(sub_data[res][1][0],3))+"{:>8}".format(round(sub_data[res][1][1],3))+"{:>8}".format(round(sub_data[res][1][2],3))+sub_data[res][3])
        
        #Calculate subject vectors and vector properties - replace header&footer info
        #Also calculate vector TS-SS similarities to reference structure
        #Store residues within similarity threshold. Format: {ref_res_index: [ [sub_res_index, vector_mag_diff] ...] }
        sub_data, conserved_res = calc_vectors_sim(sub_data, sub_distant_res_i, sub_distant_res_j, sub_COM, ref_data, args.criteria)
        
        final_alignment = defaultdict(list)

        #Remove tight-correlation pairs
        final_alignment, conserved_res = remove_tight_pairs(final_alignment, conserved_res, args.tightcriteria)

        while len(conserved_res) > 0: #Remove Single and Double-correlation pairs
            unchanged_single_double = len(final_alignment)
            
            #Align single-correlation pairs
            while len(conserved_res) > 0:
                unchanged_single = len(final_alignment)
                final_alignment, conserved_res = remove_single_pairs(final_alignment, conserved_res)
                if len(final_alignment) == unchanged_single:
                    break
            
            #Align double-correlation pairs
            final_alignment, conserved_res = remove_double_pairs(final_alignment, conserved_res)
            if len(final_alignment) == unchanged_single_double:
                break
        
        #All Single and Double-correlation pairs are exhausted. Align by sequence relative to aligned residues
        while len(conserved_res) > 0:
            unchanged = len(conserved_res)
            final_alignment, conserved_res = remove_sequential_pairs(final_alignment, conserved_res)
            if len(conserved_res) == unchanged or len(conserved_res) == 0:
                break

        #Store the results
        if args.pydist: #Geerate optional interconnected residue pymol commands
            with open(args.pydist, 'a') as logfile:
                for item in final_alignment:
                    results_dataframe[pdb_file[0].split("/")[-1]][ref_data[item][0]] = sub_data[final_alignment[item]][5]
                    results_dataframe_resnum[pdb_file[0].split("/")[-1]][ref_data[item][0]] = sub_data[final_alignment[item]][0]
                    logfile.write("distance (model aligned-{} and res {}), (model aligned-{} and res {})\n".format(args.ref[0].split("/")[-1].split(".")[0], ref_data[item][0], pdb_file[0].split("/")[-1].split(".")[0], sub_data[final_alignment[item]][0]))
        else:
            for item in final_alignment:
                results_dataframe[pdb_file[0]][ref_data[item][0]] = sub_data[final_alignment[item]][5]
                results_dataframe_resnum[pdb_file[0]][ref_data[item][0]] = sub_data[final_alignment[item]][0]
        
    #Generate the dataframes of the aligned results
    results_dataframe = pd.DataFrame(results_dataframe)
    results_dataframe['Score'] = results_dataframe.apply(lambda x: calc_score(x, scoring_matrix.granthams), axis=1)
    results_dataframe['MSCC'] = results_dataframe['Score'] * ((results_dataframe.drop('Score',axis=1).nunique(axis=1)-1) / (results_dataframe.drop('Score',axis=1).count(axis=1)-1))
    results_dataframe['SOR'] = (results_dataframe.drop(columns=['Score','MSCC'],axis=1).count(axis=1)-1)/(len(results_dataframe.keys())-2)

    results_dataframe_resnum = pd.DataFrame(results_dataframe_resnum)
    results_dataframe_resnum['MSCC'] = results_dataframe['MSCC']
    results_dataframe_resnum['MSCC_norm'] = ( (results_dataframe_resnum['MSCC']-results_dataframe_resnum['MSCC'].min()) / (results_dataframe_resnum['MSCC'].max()-results_dataframe_resnum['MSCC'].min()))*100
    
    #Save results to xlsx file
    with pd.ExcelWriter(args.save) as savefile:
        results_dataframe.to_excel(savefile, sheet_name='ResID', index=False)
        results_dataframe_resnum.to_excel(savefile, sheet_name='ResNum', index=False)

    #Generate optional MSSC heatmap pymol commands
    if args.pyheat:
        def write_heatmap(row, names, savefile):
            for i in range(0,len(row)-2):
                if pd.isnull(row[i]): continue
                else:
                    savefile.write("alter (model {} and res {}), b= {}\n".format(names[i].split("/")[-1].split(".")[0], int(row[i]), round(row[-1],2)))

        with open(args.pyheat, 'w') as logfile:
            results_dataframe_resnum.apply(lambda x: write_heatmap(x, x.keys(), logfile), axis=1)
