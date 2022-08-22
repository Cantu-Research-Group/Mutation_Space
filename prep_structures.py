import os, argparse



def find_nth(haystack, needle, n):
	start = haystack.find(needle)
	while start >= 0 and n > 1:
		start = haystack.find(needle, start+len(needle))
		n -= 1
	return start


def remove_DS(some_list):
	try:
		some_list.remove('.DS_Store')
	except ValueError:
		pass


amino_acids = ['ALA', 'THR', 'LEU', 'ILE', 'GLU', 'PRO', 'CYS', 'TRP', 'ASN',
			   'VAL', 'ARG', 'GLY', 'GLN', 'SER', 'ASP', 'LYS', 'PHE', 'HIS',
			   'TYR', 'MET', 'MSE']

fam = 'KS4'

source = 'C:/Users/k1l3r/Desktop/Mutation_space/'+fam+'/'+fam+'_pdbs/'
destination = 'C:/Users/k1l3r/Desktop/Mutation_space/'+fam+'/'+fam+'_A_pdbs/'


files = os.listdir(source)
for file in files:
	start_index = 0
	end_index = 0
	with open(source+file) as fopen:
		flines = fopen.readlines()
	for c, line in enumerate(flines):
		if start_index == 0 and (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20] in amino_acids:
			start_index = c
		elif end_index == 0 and line.startswith('TER') and line[17:20] in amino_acids:
			end_index = c+1
	flines = flines[start_index:end_index]

	temp = []
	with open(destination+file[:-4]+'_A.pdb', 'a') as output:
		for line in flines:
			if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('TER'):
				temp.append(line)

		for line in temp:
			if line.startswith('MODEL'):
				output.write(line)
			elif line[13:15] == 'CA' or line.startswith('TER'):
				output.write(line)
			elif line.startswith('ENDMDL'):
				output.close()
