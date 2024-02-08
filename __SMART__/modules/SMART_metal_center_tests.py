import sys, os
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')
pt = Chem.GetPeriodicTable()

structure = 'R-TCPTTL_1.sdf'
metal_id = 1
binding_atoms = [2, 22, 51, 71, 102]
geometry = 'octahedral'

supp = Chem.SDMolSupplier(structure, removeHs = False)
for mol in supp:
    print(mol.GetNumAtoms())
    for atom in mol.GetAtoms():
        if atom.GetIdx() == metal_id-1:
            print(atom.GetSymbol(), atom.GetIdx())
            for a in binding_atoms:

                for b in binding_atoms:
                    if a == b:
                        continue
                    angle = rdMolTransforms.GetAngleDeg(mol.GetConformer(), a-1, metal_id-1, b-1)
                    
                    print(a, metal_id, b, angle)
            #for n in atom.GetNeighbors():
            #    print(n.GetSymbol(), n.GetIdx())
