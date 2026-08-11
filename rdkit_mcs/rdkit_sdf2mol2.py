import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import glob, os, sys
import Mol2writer


def rename_mol2(inp):
	filename = f'{inp}.mol2'
	fp=open(filename,'r')
	gp=open('tmp.mol2','w')
	c_count = 1 
	h_count = 1
	o_count = 1
	f_count = 1
	n_count = 1
	s_count = 1
	cl_count = 1
	line=fp.readline()
	while line:
		if line[0:13] == '@<TRIPOS>ATOM':
			gp.write(line)
			line=fp.readline()
			while line[0:13] != '@<TRIPOS>BOND':
				tmp=line.split()
				if tmp[1] == 'C':
					tmp[1] = tmp[1]+str(c_count)
					c_count += 1
				if tmp[1] == 'O':
					tmp[1] = tmp[1]+str(o_count)
					o_count += 1
				if tmp[1] == 'N':
					tmp[1] = tmp[1]+str(n_count)
					n_count += 1
				if tmp[1] == 'S':
					tmp[1] = tmp[1]+str(s_count)
					s_count += 1
				if tmp[1] == 'F':
					tmp[1] = tmp[1]+str(f_count)
					f_count += 1
				if tmp[1] == 'H':
					tmp[1] = tmp[1]+str(h_count)
					h_count += 1
				if tmp[1] == 'Cl':
					tmp[1] = tmp[1]+str(cl_count)
					cl_count += 1
				#print(len(tmp),tmp)
				gp.write("%7s %-6s  %10s%10s%10s %-6s%3s  %-4s     %9s\n" % (tmp[0],tmp[1],
						tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8]))
				line=fp.readline()

		else:
			gp.write(line)
			line=fp.readline()

	fp.close()
	gp.close()
	os.system(f"mv tmp.mol2 {inp}.mol2")
	#os.system(f"rm tmp.mol2")

def unique_atoms(mol):
	atom_list = []
	with open(f'{mol}.mol2') as f2:
		reading_atoms = False
		for line in f2:
			if line.startswith("@<TRIPOS>ATOM"):
				reading_atoms = True
			if line.startswith("@<TRIPOS>BOND"):
				reading_atoms = False
			if reading_atoms:
				tmp = line.split()
				if tmp:
					#print(tmp)
					try:
						atom_list.append(tmp[1])
					except IndexError: None
						
	if len(set(atom_list)) != len(atom_list):
		print(f"Duplicate atom names found for {mol}.")
		quit()

#=============================================================================================#
#=============================================================================================#
if os.path.exists('Mol2writer.py'):
	None
else:
	print("We need the function for writing Mol2 File.")
	quit()


with open('mol_list.txt', 'r') as f0:
	for lines in f0.readlines():
		mol_name = str(lines.split('\n')[0]).split(' ')[0]
		lig = mol_name
		
		print(f'processing ligand named: {lig}')
		sdf = Chem.SDMolSupplier(f'{lig}.sdf')

		Mol2writer.MolToMol2File(sdf[0], f'{lig}.mol2', confId=-1)
		rename_mol2(lig)
		unique_atoms(lig)

quit()
