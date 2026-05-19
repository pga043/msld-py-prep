import rdkit
print('\n-------------------------------')
print('\n rdkit Version : ', rdkit.__version__)
print('\n-------------------------------')
from rdkit import Chem
import numpy as np
from rdkit.Chem import rdFMCS


mol1 = Chem.MolFromPDBFile('EOS14090.pdb', removeHs=False)
mol2 = Chem.MolFromPDBFile('EOS22516.pdb', removeHs=False)
mol3 = Chem.MolFromPDBFile('EOS15144.pdb', removeHs=False)
mol4 = Chem.MolFromPDBFile('EOS14090P.pdb', removeHs=False)

mols = [mol1, mol2, mol3, mol4,] # mol5] 

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


params = rdFMCS.MCSParameters()
params.AtomTyper = rdFMCS.AtomCompare.CompareElements
params.BondTyper = rdFMCS.BondCompare.CompareOrder
params.BondCompareParameters.RingMatchesRingOnly = True
params.BondCompareParameters.CompleteRingsOnly = True

res = rdFMCS.FindMCS(mols, params)

highlightAtoms_mol1 = mol1.GetSubstructMatch(res.queryMol)
highlightAtoms_mol2 = mol2.GetSubstructMatch(res.queryMol)
highlightAtoms_mol3 = mol3.GetSubstructMatch(res.queryMol)
highlightAtoms_mol4 = mol4.GetSubstructMatch(res.queryMol)


#print('\n\n')
#print('res.queryMol : ' , res.queryMol)
#print('res.queryMol.GetAtoms : ' , [ i.GetIdx() for i in res.queryMol.GetAtoms()])
#print('\n\n')
#print('highlightAtoms_mol_1 : ' , highlightAtoms_mol1)

#print('\n\n')
#print('highlightAtoms_mol_2 : ' , highlightAtoms_mol2)
#print('\n\n')

list_mol1 = [atom_by_index(mol1, atom) for atom in highlightAtoms_mol1]
list_mol2 = [atom_by_index(mol2, atom) for atom in highlightAtoms_mol2]
list_mol3 = [atom_by_index(mol3, atom) for atom in highlightAtoms_mol3]
list_mol4 = [atom_by_index(mol4, atom) for atom in highlightAtoms_mol4]

for i in range(len(list_mol1)):
    #print(list_mol1[i], list_mol2[i], list_mol3[i], list_mol4[i]) #, list_mol5[i])
	print(f'   cats sele atom LIG 1 {list_mol1[i]} .or. atom L1 1 {list_mol2[i]} .or. atom L2 1 {list_mol3[i]} .or. atom L3 1 {list_mol4[i]} end')

#for i in range(len(list_mol1)):
#    #print(list_mol1[i], list_mol2[i], list_mol3[i], list_mol4[i]) #, list_mol5[i])
#    print(f'   cats sele atom LIG 1 {list_mol1[i]} .or. atom L2 1 {list_mol3[i]} end')

#for i in range(len(list_mol1)):
#    #print(list_mol1[i], list_mol2[i], list_mol3[i], list_mol4[i]) #, list_mol5[i])
#    print(f'   cats sele atom LIG 1 {list_mol1[i]} .or. atom L3 1 {list_mol4[i]} end')


#list_of_tuples = list(zip(list_mol1, list_mol2, list_mol3, list_mol4)) #, list_mol5))
#print(list_of_tuples)


quit()
#==================================================================================#

from rdkit import Chem
from rdkit.Chem.Draw import  MolDraw2DCairo
from io import BytesIO
from PIL import Image
import cv2
import numpy as np



mols = [Chem.MolFromSmiles(smi) for smi in [
    "Oc3cc(C1CCCCC1)cc(C2CCCC2)c3",
    "NC4CCC(c3cc(O)cc(C2CCCc1ccccc12)c3)C4"]
]


mol_1 = mols[0]

mol_1.SetProp("_Name","mol_1")

mol_2 = mols[1]

mol_2.SetProp("_Name","mol_2")



def show_mol(d2d,mol,legend='',highlightAtoms=[]):
    d2d.DrawMolecule(mol,legend=legend, highlightAtoms=highlightAtoms)
    d2d.FinishDrawing()
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)
def show_images(imgs,buffer=5):
    height = 0
    width = 0
    for img in imgs:
        height = max(height,img.height)
        width += img.width
    width += buffer*(len(imgs)-1)
    res = Image.new("RGBA",(width,height))
    x = 0
    for img in imgs:
        res.paste(img,(x,0))
        x += img.width + buffer
    return res

imgs = []

d2d = MolDraw2DCairo(350,300)
d2d.drawOptions().addAtomIndices = True

dopts = d2d.drawOptions()
dopts.setBackgroundColour((0,.9,.9,.3))

imgs.append(show_mol(d2d,mol_1,legend = mol_1.GetProp('_Name')))



d2d = MolDraw2DCairo(350,300)
d2d.drawOptions().addAtomIndices = True

dopts = d2d.drawOptions()
dopts.setBackgroundColour((0,.9,.9,.3))

imgs.append(show_mol(d2d,mol_2,legend = mol_2.GetProp('_Name')))




cv2.namedWindow("Input")
cv2.imshow("Input", cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))
cv2.waitKey(0)


cv2.imwrite('input.png' , cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))


from rdkit.Chem import rdFMCS


params = rdFMCS.MCSParameters()
params.AtomTyper = rdFMCS.AtomCompare.CompareElements
params.BondTyper = rdFMCS.BondCompare.CompareOrder
params.BondCompareParameters.RingMatchesRingOnly = True
params.BondCompareParameters.CompleteRingsOnly = True

res = rdFMCS.FindMCS(mols, params)


imgs = []

d2d = MolDraw2DCairo(350,300)
d2d.drawOptions().addAtomIndices = True

dopts = d2d.drawOptions()
dopts.setBackgroundColour((0,.9,.9,.3))

imgs.append(show_mol(d2d,res.queryMol,legend = 'res'))


cv2.namedWindow("MCS")
cv2.imshow("MCS", cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))
cv2.waitKey(0)

cv2.imwrite('mcs.png' , cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))



imgs = []

d2d = MolDraw2DCairo(350,300)
d2d.drawOptions().addAtomIndices = True

dopts = d2d.drawOptions()
dopts.setBackgroundColour((0,.9,.9,.3))


highlightAtoms_mol_1 = mol_1.GetSubstructMatch(res.queryMol)

imgs.append(show_mol(d2d,mol_1,legend = mol_1.GetProp('_Name'), highlightAtoms = highlightAtoms_mol_1 ))



d2d = MolDraw2DCairo(350,300)
d2d.drawOptions().addAtomIndices = True

dopts = d2d.drawOptions()
dopts.setBackgroundColour((0,.9,.9,.3))

highlightAtoms_mol_2 = mol_2.GetSubstructMatch(res.queryMol)

imgs.append(show_mol(d2d,mol_2,legend = mol_2.GetProp('_Name'), highlightAtoms = highlightAtoms_mol_2 ))



cv2.namedWindow("Output")
cv2.imshow("Output", cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))
cv2.waitKey(0)

cv2.imwrite('output.png' , cv2.cvtColor(np.array(show_images(imgs)), cv2.COLOR_BGR2RGB))


print('\n\n')
print('res.queryMol : ' , res.queryMol)
print('res.queryMol.GetAtoms : ' , [ i.GetIdx() for i in res.queryMol.GetAtoms()])

print('\n\n')
print('highlightAtoms_mol_1 : ' , highlightAtoms_mol_1)

print('\n\n')
print('highlightAtoms_mol_2 : ' , highlightAtoms_mol_2)

print('\n\n')
list_of_tuples = list(zip(highlightAtoms_mol_1 ,highlightAtoms_mol_2)) 
print(list_of_tuples)

