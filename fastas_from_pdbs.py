import os


def load_multiprot_lines(file):
	"""this takes and aligned multiprot Ca only file and loads the lines to
	memory, averaging any partial occupancy lines for consistency"""
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	numbers = '123456789'
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

family = 'KS4'
source = 'C:/Users/k1l3r/Desktop/Mutation_space/'+family+'/'+family+'_A_pdbs/'
output = 'C:/Users/k1l3r/Desktop/Mutation_space/'+family+'/'+family+'_fastas_from_pdbs.txt'
three_letter = ['ALA', 'THR', 'LEU', 'ILE', 'GLU', 'PRO', 'CYS', 'TRP', 'ASN',
			   'VAL', 'ARG', 'GLY', 'GLN', 'SER', 'ASP', 'LYS', 'PHE', 'HIS',
			   'TYR', 'MET', 'MSE', 'SCY', 'OCS', 'CSO', 'CSJ', 'CSD']
two_letter = ['A', 'T', 'L', 'I', 'E', 'P', 'C', 'W', 'N', 'V', 'R', 'G', 'Q', 'S',
		  'D', 'K', 'F', 'H', 'Y', 'M', 'M', 'C', 'C', 'C', 'C', 'C']

abbreviations = {}
for r in range(len(three_letter)):
	abbreviations[three_letter[r]] = two_letter[r]

if '.DS_Store' in os.listdir(source):
	os.remove(source+'.DS_Store')

fastas = {}
for file in os.listdir(source):
	fastas[file] = ''
	lines = load_multiprot_lines(source+file)
	for line in lines:
		res = line[17:20]
		fastas[file] += abbreviations[res]

with open(output, 'w') as oopen:
	for k, v in fastas.items():
		oopen.write('>' + k + '\n')
		oopen.write(v + '\n')
