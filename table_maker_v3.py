#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser(description='Generates table of grantham distance for residues in matches.')
parser.add_argument('pdb', type=str, help='Path location directory for where the pdbs are stored ')
parser.add_argument('matches', type=str, help='Path location for where the matches are stored')
parser.add_argument('ref', type=str, help='Name of the pdb file to treat as a reference')
parser.add_argument('-o', type=str, default='distances.csv', help='Name of the output .csv file (default: distances.csv)')
args = parser.parse_args()


def load_multiprot_lines(file):
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
	return output

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
granthams = {'SER': {'ARG':110, 'LEU':145, 'PRO':74, 'THR':58, 'ALA':99, 'VAL':124,
					 'GLY':56, 'ILE':142, 'PHE':155, 'TYR':144, 'CYS':112, 'HIS':89,
					 'GLN':68, 'ASN':46, 'LYS':121, 'ASP':65, 'GLU':80, 'MET':135,
					 'TRP':177, 'SER':0, '':0, 'MSE':135},
			 'ARG':{'SER':110, 'LEU':102, 'PRO':103, 'THR':71, 'ALA':112, 'VAL':96,
					 'GLY':125, 'ILE':97, 'PHE':97, 'TYR':77, 'CYS':180, 'HIS':29,
					 'GLN':43, 'ASN':86, 'LYS':26, 'ASP':96, 'GLU':54, 'MET':91,
					 'TRP':101, 'ARG':0, '':0, 'MSE':91},
			 'LEU':{'SER':145, 'ARG':102, 'PRO':98, 'THR':92, 'ALA':96, 'VAL':32,
					 'GLY':138, 'ILE':5, 'PHE':22, 'TYR':36, 'CYS':198, 'HIS':99,
					 'GLN':113, 'ASN':153, 'LYS':107, 'ASP':172, 'GLU':138, 'MET':15,
					 'TRP':61, 'LEU':0, '':0, 'MSE':15},
			 'PRO':{'SER':74,'ARG':103, 'LEU':98, 'THR':38, 'ALA':27, 'VAL':68,
					 'GLY':42, 'ILE':95, 'PHE':114, 'TYR':110, 'CYS':169, 'HIS':77,
					 'GLN':76, 'ASN':91, 'LYS':103, 'ASP':108, 'GLU':93, 'MET':87,
					 'TRP':147, 'PRO':0, '':0, 'MSE':87},
			 'THR':{'SER':58, 'ARG':71, 'LEU':92, 'PRO':38, 'ALA':58, 'VAL':69,
					 'GLY':59, 'ILE':89, 'PHE':103, 'TYR':92, 'CYS':149, 'HIS':47,
					 'GLN':42, 'ASN':65, 'LYS':78, 'ASP':85, 'GLU':65, 'MET':81,
					 'TRP':128, 'THR':0, '':0, 'MSE':81},
			 'ALA':{'SER':99, 'ARG':112, 'LEU':96, 'PRO':27, 'THR':58, 'VAL':64,
					 'GLY':60, 'ILE':94, 'PHE':113, 'TYR':112, 'CYS':195, 'HIS':86,
					 'GLN':91, 'ASN':111, 'LYS':106, 'ASP':126, 'GLU':107, 'MET':84,
					 'TRP':148, 'ALA':0, '':0 ,'MSE':84},
			 'VAL':{'SER':124, 'ARG':96, 'LEU':32, 'PRO':68, 'THR':69, 'ALA':64,
					 'GLY':109, 'ILE':29, 'PHE':50, 'TYR':55, 'CYS':192, 'HIS':84,
					 'GLN':96, 'ASN':133, 'LYS':97, 'ASP':152, 'GLU':121, 'MET':21,
					 'TRP':88, 'VAL':0, '':0, 'MSE':21},
			 'GLY':{'SER':56, 'ARG':125, 'LEU':138, 'PRO':42, 'THR':59, 'ALA':60,
					'VAL':109, 'ILE':135, 'PHE':153, 'TYR':147, 'CYS':159, 'HIS':98,
					'GLN':87, 'ASN':80, 'LYS':127, 'ASP':94, 'GLU':98, 'MET':127,
					'TRP':184, 'GLY':0, '':0, 'MSE':127},
			 'ILE':{'SER':142, 'ARG':97, 'LEU':5, 'PRO':95, 'THR':89, 'ALA':94,
					'VAL':29, 'GLY':135, 'PHE':21, 'TYR':33, 'CYS':198, 'HIS':94,
					'GLN':109, 'ASN':149, 'LYS':102, 'ASP':168, 'GLU':134, 'MET':10,
					'TRP':61, 'ILE':0, '':0, 'MSE':10},
			 'PHE':{'SER':155, 'ARG':97, 'LEU':22, 'PRO':114, 'THR':103, 'ALA':113,
					'VAL':50, 'GLY':153, 'ILE':21, 'TYR':22, 'CYS':205, 'HIS':100,
					'GLN':116, 'ASN':158, 'LYS':102, 'ASP':177, 'GLU':140, 'MET':28,
					'TRP':40, 'PHE':0, '':0, 'MSE':28},
			 'TYR':{'SER':144, 'ARG':77, 'LEU':36, 'PRO':110, 'THR':92, 'ALA':112,
					'VAL':55, 'GLY':147, 'ILE':33, 'PHE':22, 'CYS':194, 'HIS':83,
					'GLN':99, 'ASN':143, 'LYS':85, 'ASP':160, 'GLU':122, 'MET':36,
					'TRP':37, 'TYR':0, '':0, 'MSE':36},
			 'CYS':{'SER':112, 'ARG':180, 'LEU':198, 'PRO':169, 'THR':149, 'ALA':195,
					'VAL':192, 'GLY':159, 'ILE':198, 'PHE':205, 'TYR':194, 'HIS':174,
					'GLN':154, 'ASN':139, 'LYS':202, 'ASP':154, 'GLU':170, 'MET':196,
					'TRP':215, 'CYS':0, '':0, 'MSE':196},
			 'HIS':{'SER':89, 'ARG':29, 'LEU':99, 'PRO':77, 'THR':47, 'ALA':86,
					'VAL':84, 'GLY':98, 'ILE':94, 'PHE':100, 'TYR':83, 'CYS':174,
					'GLN':24, 'ASN':68, 'LYS':32, 'ASP':81, 'GLU':40, 'MET':87,
					'TRP':115, 'HIS':0, '':0, 'MSE':87},
			 'GLN':{'SER':68, 'ARG':43, 'LEU':113, 'PRO':76, 'THR':42, 'ALA':91,
					'VAL':96, 'GLY':87, 'ILE':109, 'PHE':116, 'TYR':99, 'CYS':154,
					'HIS':24, 'ASN':46, 'LYS':53, 'ASP':61, 'GLU':29, 'MET':101,
					'TRP':130, 'GLN':0, '':0, 'MSE':101},
			 'ASN':{'SER':46, 'ARG':86, 'LEU':153, 'PRO':91, 'THR':65, 'ALA':111,
					'VAL':133, 'GLY':80, 'ILE':149, 'PHE':158, 'TYR':143, 'CYS':139,
					'HIS':68, 'GLN':46, 'LYS':94, 'ASP':23, 'GLU':42, 'MET':142,
					'TRP':174, 'ASN':0, '':0, 'MSE':142},
			 'LYS':{'SER':121, 'ARG':26, 'LEU':107, 'PRO':103, 'THR':78, 'ALA':106,
					'VAL':97, 'GLY':127, 'ILE':102, 'PHE':102, 'TYR':85, 'CYS':202,
					'HIS':32, 'GLN':53, 'ASN':94, 'ASP':101, 'GLU':56, 'MET':95,
					'TRP':110, 'LYS':0, '':0, 'MSE':95},
			 'ASP':{'SER':65, 'ARG':96, 'LEU':172, 'PRO':108, 'THR':85, 'ALA':126,
					'VAL':152, 'GLY':94, 'ILE':168, 'PHE':177, 'TYR':160, 'CYS':154,
					'HIS':81, 'GLN':61, 'ASN':23, 'LYS':101, 'GLU':45, 'MET':160,
					'TRP':181, 'ASP':0, '':0, 'MSE':160},
			 'GLU':{'SER':80, 'ARG':54, 'LEU':138, 'PRO':93, 'THR':65, 'ALA':107,
					'VAL':121, 'GLY':98, 'ILE':134, 'PHE':140, 'TYR':122, 'CYS':170,
					'HIS':40, 'GLN':29, 'ASN':42, 'LYS':56, 'ASP':45, 'MET':126,
					'TRP':152, 'GLU':0, '':0, 'MSE':126},
			 'MET':{'SER':135, 'ARG':91, 'LEU':15, 'PRO':87, 'THR':81, 'ALA':84,
					'VAL':21, 'GLY':127, 'ILE':10, 'PHE':28, 'TYR':36, 'CYS':196,
					'HIS':87, 'GLN':101, 'ASN':142, 'LYS':95, 'ASP':160, 'GLU':126,
					'TRP':67, 'MET':0, '':0, 'MSE':0},
			 'MSE':{'SER':135, 'ARG':91, 'LEU':15, 'PRO':87, 'THR':81, 'ALA':84,
					'VAL':21, 'GLY':127, 'ILE':10, 'PHE':28, 'TYR':36, 'CYS':196,
					'HIS':87, 'GLN':101, 'ASN':142, 'LYS':95, 'ASP':160, 'GLU':126,
					'TRP':67, 'MET':0, '':0, 'MSE':0},
			 'TRP':{'SER':177, 'ARG':101, 'LEU':61, 'PRO':147, 'THR':128, 'ALA':148,
					'VAL':88, 'GLY':184, 'ILE':61, 'PHE':40, 'TYR':37, 'CYS':215,
					'HIS':115, 'GLN':130, 'ASN':174, 'LYS':110, 'ASP':181, 'GLU':152,
					'MET':67, 'TRP':0, '':0, 'MSE':67}}
matches = args.matches + "/"
pdbs = args.pdb + "/"
pivot = args.ref
files = os.listdir(pdbs)
table = {'files':[pivot]}

with open(pdbs + pivot) as popen:
	plines = popen.readlines()
for line in plines:
	table[line[22:26].strip()] = [[line[22:26].strip(), line[17:20].strip()]]

for file in sorted(files):
	if file not in table['files']:
		table['files'].append(file)

for file in table['files'][1:]:
	with open(matches + pivot + '_' + file + '_correlated_residues.dat') as mopen:
		mlines = mopen.readlines()[1:-1]
	for line in mlines:
		pivot_res = line.strip().split('\t')[0]
		subject_res = line.strip().split('\t')[1]
		if subject_res == '{}':
			table[pivot_res].append(['', ''])
		else:
			flines = load_multiprot_lines(pdbs + file)
			for fline in flines:
				aa = fline[17:20].strip()
				if fline[22:26].strip() == subject_res:
					table[pivot_res].append([subject_res, aa])

for piv_res_num, matches in table.items():
	if piv_res_num != 'files':
		piv_aa = matches[0][1]
		matches[0].append(0)
		for r in range(1, len(matches)):
			sub_aa = matches[r][1]
			matches[r].append(granthams[piv_aa][sub_aa])
		score = 0
		count = 0
		for r in range(len(matches)):
			if matches[r][1]:
				count += 1
				score += matches[r][2]
		matches.append(['p_ave', float(count/len(matches)), ''])
		matches.append(['score', score, ''])

with open(args.o, 'w') as savefile:
	for i in table['file']:
		savefile.write(str(i)+",")
	savefile.write("\n")
	for k, v in table.items():
		if k != 'files':
			savefile.write(str(k)+",")
			for i in v:
				savefile.write(str(i[1])+",")
			savefile.write("\n")

#print('files', end='\t')
#for i in table['files']:
#	print(i, end='\t')
#print('')
#for k, v in table.items():
#	if k != 'files':
#		print(k, end='\t')
#		for i in v:
#			print(i[1], end='\t')
#		print('')
