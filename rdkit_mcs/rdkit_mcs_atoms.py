import rdkit
print('\n-------------------------------')
print('\n rdkit Version : ', rdkit.__version__)
print('\n-------------------------------')
from rdkit import Chem
import numpy as np
from rdkit.Chem import rdFMCS
#import joblib 
import time
import fnmatch

mol_list = open('mol_list.txt', 'r')
use_cats = False
use_hybrid = True
user_mcs = False
user_mcs_pattern = 'C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C(=O)[O-]'


'''
Things to add:
1. take care of hydrogens/F/Cl when writing selections for CATS - done
2. user needs to specify the resid/segid for each mol for CATS
3. most of the time (unless it casues ring breaking) have at least one heavy atom as part of the substituent 

Possible failures:
1. Degenrate MCS
2.  
'''
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
mols = []
names = []
for pdb in mol_list:
	mols.append(Chem.MolFromPDBFile(f'{pdb.split()[0]}.pdb', removeHs=False))
	names.append(f'{pdb.split()[0]}')

nmols = len(mols)

def atom_by_index(pdb, index):
    at = pdb.GetAtomWithIdx(index)
    name = at.GetPDBResidueInfo().GetName()
    return name
    
def get_atom_name(pdb):
    idx = 0
    atoms = []
    for atom in pdb.GetAtoms():
        at = pdb.GetAtomWithIdx(idx)
        ri = at.GetPDBResidueInfo()
        #print(ri.GetName())
        #print(ri.GetResidueName())
        idx += 1
        atoms.append(ri.GetName())
    return atoms

#===============================================================================#
#-------------------------- MCS Rules ------------------------------------------#
#===============================================================================#
params = rdFMCS.MCSParameters()
params.AtomTyper = rdFMCS.AtomCompare.CompareElements
params.BondTyper = rdFMCS.BondCompare.CompareOrder
#params.AtomTyper = rdFMCS.AtomCompare.CompareAnyHeavyAtom
params.BondCompareParameters.RingMatchesRingOnly = True
params.BondCompareParameters.CompleteRingsOnly = True
params.BondCompareParameters.MatchFusedRings = True
#===============================================================================#

if user_mcs != True:
	start = time.time()
	res = rdFMCS.FindMCS(mols, params, timeout=3600)
	elapsed_time = time.time() - start
	print(f"Time taken for MCS search: {elapsed_time:.4f} seconds")
	print('\n')

highlightAtoms = []
atoms_list = []
site_atoms = []
anchor_atoms = []

def get_neighbour(mol, idx):
	atoms_to_look = []
	potential_neighbours = mol.GetAtomWithIdx(idx).GetNeighbors()
	#print(f'number of neighbours = {len(potential_neighbours)}')
	for neigh in potential_neighbours:
		#print(neigh.GetSymbol(), neigh.GetIdx())
		if mol.GetAtomWithIdx(neigh.GetIdx()).GetMass() > 2:
			#print(mol.GetAtomWithIdx(neigh.GetIdx()).GetMass())
			atom_to_look = atom_by_index(mol, neigh.GetIdx())
			atoms_to_look.append(atom_to_look)
	#print(atoms_to_look)
	return atoms_to_look

#===============================================================================#
# get common core atoms
for mol in mols:
	if user_mcs == True:
		query_mol = Chem.MolFromSmiles(user_mcs_pattern, sanitize=False)
		hit = mol.GetSubstructMatch(query_mol)
		highlightAtoms.append(hit)
	else:
		hit = mol.GetSubstructMatch(res.queryMol)
		highlightAtoms.append(hit)
	#atom_list = [atom_by_index(mol, atom) if mol.GetAtomWithIdx(atom).GetAtomicNum() != 1 else None for atom in hit]
	atom_list = [atom_by_index(mol, atom) for atom in hit]

	# get site n substituents 
	not_hit = [i for i in range(mol.GetNumAtoms()) if i not in hit]
	site_atom = [atom_by_index(mol, atom) for atom in not_hit]

	# take the heavy atom out of common core if there are only hydrogens left in the substituents
		
	
	# get anchor atom
	for atom in not_hit:
		atoms_to_look = get_neighbour(mol, atom)
		for anchor in atoms_to_look:
			if anchor in atom_list:
				#print(anchor)
				anchor_atom = anchor

	atoms_list.append(atom_list)
	site_atoms.append(site_atom)
	anchor_atoms.append(anchor_atom)


# make sure all the mols have the same number of atoms in the core
natoms_core = []
for i in range(nmols):
	natoms_core.append(int(len(atoms_list[i])))
if all(x == natoms_core[0] for x in natoms_core) != True:
	print('The common core region does not contain same number of atoms for all the molecules.')
	print(f'{natoms_core}')
	print('Exiting without writing any files.')
	quit()
#===============================================================================#

#===============================================================================#
#------------------printing info for MSLD --------------------------------------#
#===============================================================================#	
if use_cats == True:
	for i in range(len(atoms_list[0])):
    	#print(list_mol1[i], list_mol2[i], list_mol3[i], list_mol4[i]) #, list_mol5[i])
		#if str(atoms_list[0][i]) != '*H*':
		parts = []
		if fnmatch.fnmatch(str(atoms_list[0][i]), '*H*') == False:
			if fnmatch.fnmatch(str(atoms_list[0][i]), '*F*') == False:
				if fnmatch.fnmatch(str(atoms_list[0][i]), '*Cl*') == False:
			#print(f'   cats sele atom LIG1 1 {atoms_list[0][i]} .or. atom LIG2 1 {atoms_list[1][i]} .or. atom LIG3 1 {atoms_list[2][i]} .or. LIG4 1 {atoms_list[3][i]} .or. atom LIG5 1 {atoms_list[4][i]}')
					for j in range(nmols):
						parts.append(f'atom LIG{j+1} 1 {atoms_list[j][i]}')
		if parts:
			print('   cats sele ' + ' .or. '.join(parts) + ' end')


print('\n')

if use_hybrid == True:
	print(f'# file written using rdkit version {rdkit.__version__} mcs')
	print(f'# {nmols} molecules processed')
	print('')
	print(f'NSUBS {nmols}')	
	print('')
	print(f'REFLIG {names[0]}')
	print('')
	print('CORE')
	for i in range(nmols):
		print(f"{names[i]} {' '.join(map(str, atoms_list[i]))}") #{x for x in atoms_list[i]}')
	print('')
	print('ANCHOR ATOMS')
	for i in range(nmols):
		print(f'{names[i]} {anchor_atoms[i]}')
	print('')
	print('SITE 1 FRAGMENTS')
	for i in range(nmols):
		print(f"{names[i]} {' '.join(map(str, site_atoms[i]))}")
	print('')
	print('END')

quit()

'''
References:
https://greglandrum.github.io/rdkit-blog/posts/2022-06-23-3d-mcs.html
https://greglandrum.github.io/rdkit-blog/posts/2023-10-27-mcswhatsnew.html
https://github.com/rdkit/rdkit/discussions/7277
'''
#==================================================================================#
