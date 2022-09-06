#!pyton3

"""
Author: Benjamin Caswell
Created: 18Nov2021
Purpose: This will load Ca coords from pdb files, read MSA data, and produce a
table of triangulation data for each Ca using two conserved residues and the
center of mass.
"""

import os, sys
import re, argparse
import numpy

parser = argparse.ArgumentParser(description='Computes spatially correlated residues based. Spatial alignment of the proteins is done based upon positioning of most conserved distant residues reported in sequence-aligned file, though this can be disabled to use user-aligned structures.')
parser.add_argument('PDBs', type=str, help='Name of directory containing A-chain CA PDBs')
parser.add_argument('msa', type=str, help='Name of the Multiple Sequence Alignment file in ClustalW format')
parser.add_argument('ref', type=str, help='Name of the PDB file to serve as a reference for the vector analysis')
parser.add_argument('-noalign', action='store_true', help='If flag is given, the program will not align the vectors according to most distant conserved residue positioning')
args = parser.parse_args()


alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
source = args.PDBs
output = "./matches/"
msa = args.msa
pivot = args.ref

if not os.path.exists('./matches'): os.mkdir('./matches')
if not os.path.exists('./pruned_pdbs'): os.mkdir('./pruned_pdbs')

def dist_euclidian(v1, v2):
	"""return euclidian distance between two points extracted from vectors"""
	return numpy.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)


def dist_manhattan(v1, v2):
	"""return manhattan distance between two points"""
	return numpy.abs(v1[0]-v2[0]) + numpy.abs(v1[1]-v2[1])+ numpy.abs(v1[2]-v2[2])


def cos_similarity(v1, v2):
	"""return cos of the angle between two vectors. 1 = perfect, 0=way off"""
	v1_mag = numpy.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
	v2_mag = numpy.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
	return numpy.dot(v1, v2)/(v1_mag*v2_mag)


def vector_magnitude(v1):
	"""returnt the magnitude of some vector"""
	return numpy.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


def get_res_num(file, index):
	"""return the residue number from the data dictionary"""
	return data[file]['structure'][index]['Ca_number']


def load_multiprot_lines(file):
	numbers = '123456789'
	"""this takes and aligned multiprot Ca only file and loads the lines to
	memory, averaging any partial occupancy lines for consistency"""
	with open(file, 'r') as fopen:
		flines = fopen.readlines()
	#flines = flines[1:]
	output = []
	r = 0
	while r < len(flines) - 1:
		if flines[r][16] == ' ':
			output.append(flines[r])
			r += 1
		elif flines[r][16] in alphabet:
			alphabet_index = alphabet.index(flines[r][16])
			n = 0
			x_avg = 0
			y_avg = 0
			z_avg = 0
			for s in range(r + 1, len(flines)):
				if flines[s][16] != alphabet[s - r + alphabet_index]:
					n = s - r
					for c in range(n):
						x_avg += float(flines[r + c][30:38].strip()) * float(flines[r + c][56:60].strip())
						y_avg += float(flines[r + c][38:46].strip()) * float(flines[r + c][56:60].strip())
						z_avg += float(flines[r + c][46:54].strip()) * float(flines[r + c][56:60].strip())
					break
			x_out = (str(x_avg) + '00000000')[:8]
			y_out = (str(y_avg) + '00000000')[:8]
			z_out = (str(z_avg) + '00000000')[:8]

			outline = flines[r]
			outline = outline[:30] + x_out + y_out + z_out + outline[54:]
			output.append(outline)
			r += n
		elif flines[r][16] in numbers:
			number_index = numbers.index(flines[r][16])
			n = 0
			x_avg = 0
			y_avg = 0
			z_avg = 0
			for s in range(r + 1, len(flines)):
				if flines[s][16] != numbers[s - r + number_index]:
					n = s - r
					for c in range(n):
						x_avg += float(flines[r + c][30:38].strip()) * float(flines[r + c][56:60].strip())
						y_avg += float(flines[r + c][38:46].strip()) * float(flines[r + c][56:60].strip())
						z_avg += float(flines[r + c][46:54].strip()) * float(flines[r + c][56:60].strip())
					break
			x_out = (str(x_avg) + '00000000')[:8]
			y_out = (str(y_avg) + '00000000')[:8]
			z_out = (str(z_avg) + '00000000')[:8]

			outline = flines[r]
			outline = outline[:30] + x_out + y_out + z_out + outline[54:]
			output.append(outline)
			r += n
	return output


def rindex(lst, value):
	lst = lst[::-1]
	i = lst.index(value)
	lst = lst[::-1]
	return len(lst) - i - 1


def print_trig_anal():
	for combo, dat in triangulation_analysis.items():
		f1 = combo.split(':')[0]
		f2 = combo.split(':')[1]
		print(f1, f2, 'data', sep='\t')
		for ca1, dat1 in dat.items():
			print(get_res_num(f1, ca1), end='\t')
			for ca2, dist in dat1.items():
				print(get_res_num(f2, ca2), end='\t')
			print('')
	exit()


def print_final_align():
	for combo, matches in final_align.items():
		f1 = combo.split(':')[0]
		f2 = combo.split(':')[1]
		print(f1, f2, sep='\t')
		for i1, i2 in matches.items():
			print(i1, i2, sep='\t')
		print('')
	exit()


def save_final_align(folder):
	for combo, matches in final_align.items():
		f1 = combo.split(':')[0]
		f2 = combo.split(':')[1]
		with open(folder + f1 + '_' + f2 + '_correlated_residues.dat', 'w') as fopen:
			fopen.write(combo + '\n')
			for i1, i2 in matches.items():
				fopen.write(str(i1) + '\t' + str(i2) + '\n')
			fopen.write('\n')


def get_avg_mag_dif(f1, i1, f2, i2):
	i1r1 = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_1'])
	i1r2 = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_2'])
	i1rcom = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_com'])
	i2r1 = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_1'])
	i2r2 = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_2'])
	i2rcom = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_com'])
	return (numpy.abs(i1r1 - i2r1) + numpy.abs(i1r2 - i2r2) + numpy.abs(i1rcom - i2rcom))/3


def get_mag_dif_factor(f1, i1, f2, i2):
	i1r1 = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_1'])
	i1r2 = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_2'])
	i1rcom = vector_magnitude(triangulation_tables[f1]['cons_res'][i1]['v_com'])
	i2r1 = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_1'])
	i2r2 = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_2'])
	i2rcom = vector_magnitude(triangulation_tables[f2]['cons_res'][i2]['v_com'])
	r1_mag_dif = numpy.abs(i1r1 - i2r1)
	r2_mag_dif = numpy.abs(i1r2 - i2r2)
	com_mag_dif = numpy.abs(i1rcom - i2rcom)
	return numpy.sqrt((r1_mag_dif ** 2 + r2_mag_dif ** 2 + com_mag_dif ** 2) / 3)


def get_cos_sim(f1, i1, f2, i2):
	v11 = triangulation_tables[f1]['cons_res'][i1]['v_1']
	v12 = triangulation_tables[f1]['cons_res'][i1]['v_2']
	v1com = triangulation_tables[f1]['cons_res'][i1]['v_com']
	v21 = triangulation_tables[f2]['cons_res'][i2]['v_1']
	v22 = triangulation_tables[f2]['cons_res'][i2]['v_2']
	v2com = triangulation_tables[f2]['cons_res'][i2]['v_com']
	v_1_csim = cos_similarity(v11, v21)
	v_2_csim = cos_similarity(v12, v22)
	v_com_csim = cos_similarity(v1com, v2com)
	return numpy.average([v_1_csim, v_2_csim, v_com_csim])


def print_most_dist_residues():
	for file, inds in most_distant_residues.items():
		print(file, data[file]['structure'][inds[0]]['Ca_number'], data[file]['structure'][inds[1]]['Ca_number'])
	exit()


def save_aligned_pdbs(pdb_A_folder, output_folder):
	for file in os.listdir(pdb_A_folder):
		lines = load_multiprot_lines(pdb_A_folder + file)
		with open(output_folder + file[:-6] + '_pruned.pdb', 'a') as wopen:
			#wopen.write(lines[0])
			for c, line in enumerate(lines[:-1]):
				x_new = str(round(data[file]['structure'][c]['x'], 3))
				y_new = str(round(data[file]['structure'][c]['y'], 3))
				z_new = str(round(data[file]['structure'][c]['z'], 3))
				x_new = (12 - len(x_new)) * ' ' + str(x_new)
				y_new = (8 - len(y_new)) * ' ' + str(y_new)
				z_new = (8 - len(z_new)) * ' ' + str(z_new)
				newline = line[:26] + x_new + y_new + z_new + line[54:]
				wopen.write(newline)
	exit()


def print_conserved_res_num(file):
	for r in range(len(conserved_indices)):
		print(get_res_num(file, conserved_indices[r][file]))
	#exit()


# load data from pdb files
print('loading PDB files')
data = {}  # {file:{'structure':{Ca_index:{'Ca_number':#, 'res':'', 'x':#, 'y':#, 'z':#}},
#					'cutoff': #
#					'com':{'x':#, 'y':#, 'z':#},
#					'c_term':{'x':#, 'y':#, 'z':#},
#  					'n_term':{'x':#, 'y':#, 'z':#},
#  					'cons_res_1':{'x':#, 'y':#, 'z':#},
#  					'cons_res_2':{'x':#, 'y':#, 'z':#}
#  					}
for file in os.listdir(source):
	data[file] = {'structure':{}, 'cutoff':0., 'com':{}, 'cr_com':{}, 'c_term':{},
				  'n_term':{}, 'cons_res_1':{}, 'cons_res_2':{}}
	Ca_positions = load_multiprot_lines(source+file)
	num_Cas = len(Ca_positions)
	for index in range(num_Cas):
		data[file]['structure'][index] = {'Ca_number':0, 'res':'', 'x':0., 'y':0., 'z':0.}
		data[file]['structure'][index]['Ca_number'] = int(Ca_positions[index][22:26].strip())
		data[file]['structure'][index]['res'] = Ca_positions[index][17:20]
		data[file]['structure'][index]['x'] = float(Ca_positions[index][30:38].strip())
		data[file]['structure'][index]['y'] = float(Ca_positions[index][38:46].strip())
		data[file]['structure'][index]['z'] = float(Ca_positions[index][46:54].strip())

	# get cutoff distance
	adjacent_dists = []
	for index in range(len(data[file])-1):
		x1 = data[file]['structure'][index]['x']
		x2 = data[file]['structure'][index+1]['x']
		y1 = data[file]['structure'][index]['y']
		y2 = data[file]['structure'][index + 1]['y']
		z1 = data[file]['structure'][index]['z']
		z2 = data[file]['structure'][index + 1]['z']
		adjacent_dist = numpy.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
		adjacent_dists.append(adjacent_dist)
	data[file]['cutoff'] = numpy.average(adjacent_dists)

	# Set C and N coordinates
	data[file]['c_term']['x'] = float(Ca_positions[0][30:38].strip())
	data[file]['c_term']['y'] = float(Ca_positions[0][38:46].strip())
	data[file]['c_term']['z'] = float(Ca_positions[0][46:54].strip())
	data[file]['n_term']['x'] = float(Ca_positions[-1][30:38].strip())
	data[file]['n_term']['y'] = float(Ca_positions[-1][38:46].strip())
	data[file]['n_term']['z'] = float(Ca_positions[-1][46:54].strip())

	# Set center of mass (com) position
	x, y, z = [], [], []
	adjacent_dists = []
	for index in range(len(data[file]['structure'])):
		x.append(data[file]['structure'][index]['x'])
		y.append(data[file]['structure'][index]['y'])
		z.append(data[file]['structure'][index]['z'])
	data[file]['com']['x'] = numpy.average(x)
	data[file]['com']['y'] = numpy.average(y)
	data[file]['com']['z'] = numpy.average(z)

# read MSA data to find good conserved residues
print('reading MSA data')
msa_data = {'conserved':'', 'files':{}}  # {file:'ASPQLIY...', 'conserved':' *::: **.  :'}
with open(msa) as msa_open:
	msa_lines = msa_open.readlines()
for r in range(1, len(msa_lines)):
	if re.match('\S', msa_lines[r]):  # fasta line
		file = msa_lines[r][:msa_lines[r].index(' ')]
		if file not in msa_data['files'].keys():  # add first line to sequence
			msa_data['files'][file] = msa_lines[r][rindex(msa_lines[r], ' '):].strip()
		elif file in msa_data['files'].keys():  # add later lines to sequence
			msa_data['files'][file] += msa_lines[r][rindex(msa_lines[r], ' '):].strip()

		if msa_lines[r+1].startswith(' '):  # conserved residue line
			msa_data['conserved'] += msa_lines[r+1][16:-1]  # add conservation data

# determine indices of conserved residues for each file
conserved_indices = []
for r in range(len(msa_data['conserved'])):
	if msa_data['conserved'][r] == '*':
		to_add = {}
		for file in msa_data['files']:
			gaps = msa_data['files'][file][:r].count('-')
			file_index = r - gaps
			to_add[file] = file_index
		conserved_indices.append(to_add)


def print_all_conserved_res():
	for file in sorted(msa_data['files']):
		print(file, end='\t')
	print('')
	for r in range(len(conserved_indices)):
		for file in sorted(msa_data['files']):
			print(get_res_num(file=file, index=conserved_indices[r][file]), end='\t')
		print('')
	exit()
# print_all_conserved_res()

def pymol_cons_res_select():
	print('select ', end='')
	for c, file in enumerate(sorted(msa_data['files'])):
		print(file[:-6] + '_pruned///', end='')
		for r in range(len(conserved_indices)):
			print(get_res_num(file=file, index=conserved_indices[r][file]), end='')
			if r + 1 < len(conserved_indices):
				print('+', end='')
			else:
				pass
		if c + 1 < len(msa_data['files']):
			print('/ or ', end='')
		else:
			print('/')
	exit()
#pymol_cons_res_select()

# Find the two conserved residues that are furthest apart
# Also determine conserved residue com - cr_com
most_distant_residues = {}
for file in data.keys():
	cr_com_calc = {'x': [], 'y': [], 'z': []}
	most_distant_residues[file] = [0,0]
	largest_dist = 0
	for r in range(len(conserved_indices)):
		car_index = conserved_indices[r][file]
		x2 = data[file]['structure'][car_index]['x']
		y2 = data[file]['structure'][car_index]['y']
		z2 = data[file]['structure'][car_index]['z']
		cr_com_calc['x'].append(x2)  # add values to cr_com_calc to avg later
		cr_com_calc['y'].append(y2)
		cr_com_calc['z'].append(z2)
		for c in range(r+1):
			cac_index = conserved_indices[c][file]
			x1 = data[file]['structure'][cac_index]['x']
			y1 = data[file]['structure'][cac_index]['y']
			z1 = data[file]['structure'][cac_index]['z']
			dist = numpy.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
			if dist > largest_dist:
				largest_dist = dist
				most_distant_residues[file][1] = car_index
				most_distant_residues[file][0] = cac_index
	# avg com data for this file and add to data
	data[file]['cr_com']['x'] = numpy.average(cr_com_calc['x'])
	data[file]['cr_com']['y'] = numpy.average(cr_com_calc['y'])
	data[file]['cr_com']['z'] = numpy.average(cr_com_calc['z'])
# print_most_dist_residues()

# Update cons_res_1 and cons_res_2 data with furthest conserved residues
for file in most_distant_residues.keys():
	dub_1 = most_distant_residues[file][0]
	dub_2 = most_distant_residues[file][1]
	data[file]['cons_res_1']['x'] = data[file]['structure'][dub_1]['x']
	data[file]['cons_res_2']['x'] = data[file]['structure'][dub_2]['x']
	data[file]['cons_res_1']['y'] = data[file]['structure'][dub_1]['y']
	data[file]['cons_res_2']['y'] = data[file]['structure'][dub_2]['y']
	data[file]['cons_res_1']['z'] = data[file]['structure'][dub_1]['z']
	data[file]['cons_res_2']['z'] = data[file]['structure'][dub_2]['z']

if args.noalign == False:
	# Shift coords so that cons_res center of mass = (0, 0, 0)
	for file, dat in data.items():
			print('translating ' + file)
			x_trans = dat['cr_com']['x']
			y_trans = dat['cr_com']['y']
			z_trans = dat['cr_com']['z']
			# shift structure points
			for index, coords in dat['structure'].items():
				coords['x'] -= x_trans
				coords['y'] -= y_trans
				coords['z'] -= z_trans
			# shift special points
			for res_type in ['com', 'cr_com', 'cons_res_1', 'cons_res_2', 'n_term', 'c_term']:
				dat[res_type]['x'] -= x_trans
				dat[res_type]['y'] -= y_trans
				dat[res_type]['z'] -= z_trans

# Rotate points so that r1 = <0, 0, d_r1_cr_com>
	for file, dat in data.items():
		print('rotating ' + file)
	# find angle for x_axis rotation
		cry = 0  # looking down x axis,
		crz = 1  # set rotational destination to <x, 0, 1>
		py = dat['cons_res_1']['y']
		pz = dat['cons_res_1']['z']  # coords of point to be rotated
		mag_p = numpy.sqrt(py**2 + pz**2)  # mag_cr = 1 so denominator = mag_p
		dot_prod = crz * pz  # cry = 0 so z is all that affects dot prot
		th_x = numpy.arccos(dot_prod / mag_p)
		x_rot_matrix = numpy.matrix([[1, 0, 0],
									 [0, numpy.cos(th_x), -numpy.sin(th_x)],
									 [0, numpy.sin(th_x), numpy.cos(th_x)]])
		# rotate about x axis
		cr1_matrix = numpy.matrix([[data[file]['cons_res_1']['x']],
								   [data[file]['cons_res_1']['y']],
								   [data[file]['cons_res_1']['z']]])
		x_rotated_cr1_matrix = numpy.dot(x_rot_matrix, cr1_matrix)
		q=0
		if -0.01 <= x_rotated_cr1_matrix[1] <= 0.01:
			data[file]['cons_res_1']['x'] = float(x_rotated_cr1_matrix[0])
			data[file]['cons_res_1']['y'] = float(x_rotated_cr1_matrix[1])
			data[file]['cons_res_1']['z'] = float(x_rotated_cr1_matrix[2])
			for c in ['cons_res_2', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				x_rotated_pt_matrix = numpy.dot(x_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(x_rotated_pt_matrix[0])
				data[file][c]['y'] = float(x_rotated_pt_matrix[1])
				data[file][c]['z'] = float(x_rotated_pt_matrix[2])
		else:
			th_x = -th_x
			x_rot_matrix = numpy.matrix([[1, 0, 0],
										 [0, numpy.cos(th_x), -numpy.sin(th_x)],
										 [0, numpy.sin(th_x), numpy.cos(th_x)]])
			for c in ['cons_res_1', 'cons_res_2', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				x_rotated_pt_matrix = numpy.dot(x_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(x_rotated_pt_matrix[0])
				data[file][c]['y'] = float(x_rotated_pt_matrix[1])
				data[file][c]['z'] = float(x_rotated_pt_matrix[2])
		for index, coords in dat['structure'].items():
			pt_matrix = numpy.matrix([[coords['x']],
									  [coords['y']],
									  [coords['z']]])
			x_rotated_pt_matrix = numpy.dot(x_rot_matrix, pt_matrix)
			coords['x'] = float(x_rotated_pt_matrix[0])
			coords['y'] = float(x_rotated_pt_matrix[1])
			coords['z'] = float(x_rotated_pt_matrix[2])

	# find angle for y_axis rotation
		crx = 0
		crz = 1
		px = dat['cons_res_1']['x']
		pz = dat['cons_res_1']['z']
		mag_p = numpy.sqrt(px ** 2 + pz ** 2)
		dot_prod = crz * pz
		th_y = numpy.arccos(dot_prod / mag_p)
		y_rot_matrix = numpy.matrix([[numpy.cos(th_y), 0, numpy.sin(th_y)],
									 [0, 1, 0],
									 [-numpy.sin(th_y), 0, numpy.cos(th_y)]])
		# rotate about y axis
		cr1_matrix = numpy.matrix([[data[file]['cons_res_1']['x']],
								   [data[file]['cons_res_1']['y']],
								   [data[file]['cons_res_1']['z']]])
		y_rotated_cr1_matrix = numpy.dot(y_rot_matrix, cr1_matrix)
		q=0
		if -0.01 <= y_rotated_cr1_matrix[0] <= 0.01:
			data[file]['cons_res_1']['x'] = float(y_rotated_cr1_matrix[0])
			data[file]['cons_res_1']['y'] = float(y_rotated_cr1_matrix[1])
			data[file]['cons_res_1']['z'] = float(y_rotated_cr1_matrix[2])
			for c in ['cons_res_2', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				y_rotated_pt_matrix = numpy.dot(y_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(y_rotated_pt_matrix[0])
				data[file][c]['y'] = float(y_rotated_pt_matrix[1])
				data[file][c]['z'] = float(y_rotated_pt_matrix[2])
		else:
			th_y = -th_y
			y_rot_matrix = numpy.matrix([[numpy.cos(th_y), 0, numpy.sin(th_y)],
										 [0, 1, 0],
										 [-numpy.sin(th_y), 0, numpy.cos(th_y)]])
			for c in ['cons_res_1', 'cons_res_2', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				y_rotated_pt_matrix = numpy.dot(y_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(y_rotated_pt_matrix[0])
				data[file][c]['y'] = float(y_rotated_pt_matrix[1])
				data[file][c]['z'] = float(y_rotated_pt_matrix[2])
		for index, coords in dat['structure'].items():
			pt_matrix = numpy.matrix([[coords['x']],
									  [coords['y']],
									  [coords['z']]])
			y_rotated_pt_matrix = numpy.dot(y_rot_matrix, pt_matrix)
			coords['x'] = float(y_rotated_pt_matrix[0])
			coords['y'] = float(y_rotated_pt_matrix[1])
			coords['z'] = float(y_rotated_pt_matrix[2])
		# find angle for z_axis rotation
		crx = 1
		cry = 0
		px = dat['cons_res_2']['x']
		py = dat['cons_res_2']['y']
		mag_p = numpy.sqrt(px ** 2 + py ** 2)
		dot_prod = crx * px
		th_z = numpy.arccos(dot_prod / mag_p)
		z_rot_matrix = numpy.matrix([[numpy.cos(th_z), -numpy.sin(th_z), 0],
									 [numpy.sin(th_z), numpy.cos(th_z), 0],
									 [0, 0, 1]])
		# rotate about z axis
		cr2_matrix = numpy.matrix([[data[file]['cons_res_2']['x']],
								   [data[file]['cons_res_2']['y']],
								   [data[file]['cons_res_2']['z']]])
		z_rotated_cr2_matrix = numpy.dot(z_rot_matrix, cr2_matrix)
		q=0
		if -0.01 <= z_rotated_cr2_matrix[1] <= 0.01:
			data[file]['cons_res_2']['x'] = float(z_rotated_cr2_matrix[0])
			data[file]['cons_res_2']['y'] = float(z_rotated_cr2_matrix[1])
			data[file]['cons_res_2']['z'] = float(z_rotated_cr2_matrix[2])
			for c in ['cons_res_1', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				z_rotated_pt_matrix = numpy.dot(z_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(z_rotated_pt_matrix[0])
				data[file][c]['y'] = float(z_rotated_pt_matrix[1])
				data[file][c]['z'] = float(z_rotated_pt_matrix[2])
		else:
			th_z = -th_z
			z_rot_matrix = numpy.matrix([[numpy.cos(th_z), -numpy.sin(th_z), 0],
										 [numpy.sin(th_z), numpy.cos(th_z), 0],
										 [0, 0, 1]])
			for c in ['cons_res_1', 'cons_res_2', 'n_term', 'c_term', 'com']:
				pt_matrix = numpy.matrix([[data[file][c]['x']],
										  [data[file][c]['y']],
										  [data[file][c]['z']]])
				z_rotated_pt_matrix = numpy.dot(z_rot_matrix, pt_matrix)
				data[file][c]['x'] = float(z_rotated_pt_matrix[0])
				data[file][c]['y'] = float(z_rotated_pt_matrix[1])
				data[file][c]['z'] = float(z_rotated_pt_matrix[2])
		for index, coords in dat['structure'].items():
			pt_matrix = numpy.matrix([[coords['x']],
									  [coords['y']],
									  [coords['z']]])
			z_rotated_pt_matrix = numpy.dot(z_rot_matrix, pt_matrix)
			coords['x'] = float(z_rotated_pt_matrix[0])
			coords['y'] = float(z_rotated_pt_matrix[1])
			coords['z'] = float(z_rotated_pt_matrix[2])
# calculate triangulation tables
triangulation_tables = {}  # {file:{'c_n_com':{Ca_index:{'d_c':#, 'd_n':#, 'd_com':#, 'Ca_number':#}}
#							'cons_res':{Ca_index:{'d1':#, 'd2':#, 'd_com':#, 'Ca_number':#}}
#							}
for file in data.keys():
	print('calculating ' + file + ' triangulation tables')
	triangulation_tables[file] = {'cons_res': {}}
	# calc cons_res tables
	for index in range(len(data[file]['structure'])):
		triangulation_tables[file]['cons_res'][index] = {'v_1': (), 'v_2': (), 'v_com': ()}

		index_x = data[file]['structure'][index]['x']
		index_y = data[file]['structure'][index]['y']
		index_z = data[file]['structure'][index]['z']

		cr1_x = data[file]['cons_res_1']['x']
		cr1_y = data[file]['cons_res_1']['y']
		cr1_z = data[file]['cons_res_1']['z']

		cr2_x = data[file]['cons_res_2']['x']
		cr2_y = data[file]['cons_res_2']['y']
		cr2_z = data[file]['cons_res_2']['z']

		cr_com_x = data[file]['cr_com']['x']
		cr_com_y = data[file]['cr_com']['y']
		cr_com_z = data[file]['cr_com']['z']
		"""
		d1 = numpy.sqrt((cr1_x - index_x) ** 2 + (cr1_y - index_y) ** 2 + (cr1_z - index_z) ** 2)
		d2 = numpy.sqrt((cr2_x - index_x) ** 2 + (cr2_y - index_y) ** 2 + (cr2_z - index_z) ** 2)
		dcom = numpy.sqrt((com_x - index_x) ** 2 + (com_y - index_y) ** 2 + (com_z - index_z) ** 2)
		"""
		cr1_to_pt_vect = numpy.array([index_x - cr1_x, index_y - cr1_y, index_z - cr1_z])
		cr2_to_pt_vect = numpy.array([index_x - cr2_x, index_y - cr2_y, index_z - cr2_z])
		com_to_pt_vect = numpy.array([index_x - cr_com_x, index_y - cr_com_y, index_z - cr_com_z])

		triangulation_tables[file]['cons_res'][index]['v_1'] = cr1_to_pt_vect
		triangulation_tables[file]['cons_res'][index]['v_2'] = cr2_to_pt_vect
		triangulation_tables[file]['cons_res'][index]['v_com'] = com_to_pt_vect
# Establish file-pairs for comparison
triangulation_analysis = {}
files = list(data.keys())
files.remove(pivot)
combos = []
for file in files:
	combo = [pivot, file]
	combos.append(combo)


def singles_cycle():
	print('singles')
	while True:
		changed = False
		todel = set()  # list of ca2 indices that had a 1:1 ca1:ca2 match
		for i1, dat1 in dat.items():
			res1 = get_res_num(f1, i1)
			if len(dat1) == 1:  # if there is a 1:1 ca1:ca2 match
				i2 = list(dat1.keys())[0]
				res2 = get_res_num(f2, i2)
				temp_align[res1] = res2
				todel.add(i2)
		for i1, dat1 in dat.items():
			if dat1.keys():
				overlap = set(dat1.keys()).intersection(todel)
				if overlap:
					for index in overlap:
						dat1.pop(index)
					changed = True
					main_changed = True
		if not changed:
			for i1, i2 in temp_align.items():
				final_align[combo][i1] = i2
			break


def doubles_cycle():
	print('doubles')
	while True:
		changed = False
		todel = set()
		for i1 in range(len(dat.keys())):
			if len(dat[i1]) == 2:
				res1 = get_res_num(f1, i1)
				match1 = sorted(list(dat[i1].keys()))
				for i2 in range(i1 + 1, len(dat.keys())):
					res2 = get_res_num(f1, i2)
					match2 = sorted(list(dat[i2].keys()))
					if match2 == match1:
						dub_1 = match1[0]
						dub_2 = match1[1]
						if dat[i1][dub_1]['mag_dif_factor'] + dat[i2][dub_2]['mag_dif_factor'] < \
								dat[i1][dub_2]['mag_dif_factor'] + dat[i2][dub_1]['mag_dif_factor']:
							todel.add(dub_1)
							todel.add(dub_2)
							temp_align[res1] = get_res_num(f2, dub_1)
							temp_align[res2] = get_res_num(f2, dub_2)
						else:
							todel.add(dub_1)
							todel.add(dub_2)
							temp_align[res1] = get_res_num(f2, dub_2)
							temp_align[res2] = get_res_num(f2, dub_1)

		for i1, dat1 in dat.items():
			if dat1.keys():
				overlap = set(dat1.keys()).intersection(todel)
				if overlap:
					for index in overlap:
						dat1.pop(index)
					changed = True
					main_changed = True
		if not changed:
			for i1, i2 in temp_align.items():
				final_align[combo][i1] = i2
			break


# Add possible matches to triangulation analysis
for combo in combos:  # for each possible file combination
	f1 = combo[0]
	f2 = combo[1]

	triangulation_analysis[f1 + ':' + f2] = {}

	for ca1_index, data1 in triangulation_tables[f1]['cons_res'].items():
		triangulation_analysis[f1+':'+f2][ca1_index] = {}
		v11 = data1['v_1']
		v12 = data1['v_2']
		v1com = data1['v_com']
		v11_mag = vector_magnitude(v11)
		v12_mag	= vector_magnitude(v12)
		v1com_mag = vector_magnitude(v1com)

		for ca2_index, data2 in triangulation_tables[f2]['cons_res'].items():
			v21 = data2['v_1']
			v22 = data2['v_2']
			v2com = data2['v_com']
			v21_mag = vector_magnitude(v21)
			v22_mag = vector_magnitude(v22)
			v2com_mag = vector_magnitude(v2com)

			v_1_csim = cos_similarity(v11, v21)
			v_2_csim = cos_similarity(v12, v22)
			v_com_csim = cos_similarity(v1com, v2com)

			v_1_ed = dist_euclidian(v11, v21)
			v_2_ed = dist_euclidian(v12, v22)
			v_com_ed = dist_euclidian(v1com, v2com)

			v_1_mn = dist_manhattan(v11, v21)
			v_2_mn = dist_manhattan(v12, v22)
			v_com_mn = dist_manhattan(v1com, v2com)

			r1_mag_dif = numpy.abs(v11_mag-v21_mag)
			r2_mag_dif = numpy.abs(v12_mag-v22_mag)
			com_mag_dif = numpy.abs(v1com_mag-v2com_mag)

			avg_mag_dif = (r1_mag_dif + r2_mag_dif + com_mag_dif) / 3
			avg_cos_sim = numpy.average([v_1_csim, v_2_csim, v_com_csim])
			mag_dif_factor = numpy.sqrt((r1_mag_dif**2 + r2_mag_dif**2 + com_mag_dif**2)/3)
			# if avg_cos_sim > 0.9 and mag_dif_factor < 1:
			if mag_dif_factor < 2.6 * avg_cos_sim - 1.35:
				#print(get_res_num(f1, ca1_index), get_res_num(f2, ca2_index), v11_mag, v12_mag, v1com_mag, v21_mag, v22_mag, v2com_mag, sep='\t')
				triangulation_analysis[f1+':'+f2][ca1_index][ca2_index] = \
					{'avg_mag_dif':avg_mag_dif, 'mag_dif_factor':mag_dif_factor, 'avg_cos_sim':avg_cos_sim}

	# add most distant conserved residues to triangulation_analysis
	f1_cons_res1 = most_distant_residues[f1][0]
	f1_cons_res2 = most_distant_residues[f1][1]
	f2_cons_res1 = most_distant_residues[f2][0]
	f2_cons_res2 = most_distant_residues[f2][1]
	triangulation_analysis[f1 + ':' + f2][f1_cons_res1] = \
		{f2_cons_res1: {'avg_mag_dif': 0, 'mag_dif_factor': 0, 'avg_cos_sim': 0}}
	triangulation_analysis[f1 + ':' + f2][f1_cons_res2] = \
		{f2_cons_res2: {'avg_mag_dif': 0, 'mag_dif_factor': 0, 'avg_cos_sim': 0}}
final_align = {}
for combo, dat in triangulation_analysis.items():
	print('assigning ' + combo)
	final_align[combo] = {}
	temp_align = {}
	f1 = combo.split(':')[0]
	f2 = combo.split(':')[1]
	# main assignment cycle, singles-doubles
	while True:
		main_changed = False
		print('enter main')
		# find 1:1 matches, add to temp align, remove from triganal
		singles_cycle()

		temp_align = {}
		# find matching 1:2 matches, assign optimally by avg_dif_factor
		doubles_cycle()

		if not main_changed:
			break

	temp_align = {}
	# assign oddballs by sequence
	while True:
		print('cleanup')
		changed = False
		todel = set()  # f2 indices to remove from all positions
		for i1, dat1 in dat.items():
			toclear = False
			res1 = get_res_num(f1, i1)
			if dat1:
				matches = sorted(list(dat1.keys()))
				for ind in todel:
					if ind in matches:
						matches.remove(ind)
				match_cycle = 0
				expand_cycle = 1
				while match_cycle < len(matches):
					match_ca_number = get_res_num(f2, matches[match_cycle])
					try:
						if final_align[combo][res1 - expand_cycle] < match_ca_number < \
								final_align[combo][res1 + expand_cycle]:
							final_align[combo][res1] = match_ca_number
							todel.add(matches[match_cycle])
							toclear = True
							break
						elif (final_align[combo][res1 - expand_cycle] > match_ca_number) or \
								(final_align[combo][res1 + expand_cycle] < match_ca_number):
							match_cycle += 1
					except KeyError:
						expand_cycle += 1
						if expand_cycle > 15:  # look through 1/10th of list
							match_cycle += 1
							expand_cycle = 1
			if toclear:  # i1 was matched, clear unused 'possible' matches
				dat[i1] = {}

		for i1, dat1 in dat.items():
			if dat1.keys():
				overlap = set(dat1.keys()).intersection(todel)
				if overlap:
					for index in overlap:
						dat1.pop(index)
					changed = True

		if not changed:
			for i1, i2 in temp_align.items():
				final_align[combo][i1] = i2
			break

	# add non-matches to final_align
	for i1, dat1 in dat.items():
		if get_res_num(f1, i1) not in final_align[combo].keys():
			final_align[combo][get_res_num(f1, i1)] = {}

save_final_align(output)
save_aligned_pdbs(source, './pruned_pdbs/')
