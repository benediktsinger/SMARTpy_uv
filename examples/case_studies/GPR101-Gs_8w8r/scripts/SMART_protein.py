import os, sys
import pandas as pd
import numpy as np

sys.path.append('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/')
import SMART as smart
from SMART import conf_search as search
from SMART import descriptors as desc
import rdkit
from rdkit import Chem

props = pd.DataFrame()
residue_atoms = {'P30':
                    {'tip':[22,21,20,19,18], 'ref':[19,18], 'dist':2.0},
                'W441':
                    {'tip':[1934], 'ref':[1935,1937], 'dist':2.0},
                'T437':
                    {'tip':[1902], 'ref':[1901], 'dist':2.0}}
                #'T111':
                #    {'tip':[646,647], 'ref':[645], 'dist':2.0}}
                #{'L181':
                #    {'tip':[1205,1204], 'ref':[1203], 'dist':2.0}}
atom_pos = {}

file = '8w8r_trunc.pdb'
path = os.path.join('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/examples/case_studies/GPR101-Gs_8w8r/', file)
mol = Chem.MolFromPDBFile(path, removeHs=False)
for res in residue_atoms.keys():
    atom_pos[res] = {}
    tip, ref = [], []
    for atom in mol.GetAtoms():
        if atom.GetIdx() in residue_atoms[res]['tip']:
            tip.append(np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx())))

        if atom.GetIdx() in residue_atoms[res]['ref']:
            ref.append(np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx())))
    if res != 'P30':
        tip_ = np.array(tip).mean(axis=0)
        ref_ = np.array(ref).mean(axis=0)
        atom_pos[res]['tip'] = tip_
        atom_pos[res]['ref'] = ref_
    elif res == 'P30':
        tip_ = np.array(tip).mean(axis=0)
        atom_pos[res]['tip'] = tip_
        atom_pos[res]['ref'] = ref

probe_name = 'O_SiH2_6_linear'
if search.TEMPLATE.is_template(probe_name):
    # conformer search
    print('get template')
    search.TEMPLATE.GetTemplate(probe_name)
else:
    # make template + conformer search
    print('gen template')
    search.TEMPLATE.GenerateTemplate(probe_name)

structure = smart.ReadFile(path)
for res in atom_pos.keys():
    structure = smart.ReadFile(path)
    if res == 'P30':
        structure.reference_angle(tip_pos=atom_pos[res]['tip'], refs_pos=[atom_pos[res]['ref'][1],atom_pos[res]['ref'][0]], dist=residue_atoms[res]['dist'])
        structure.export_alignment(out=res+'_align', dummy='Mg')
    else:
        structure.reference_vector(tip_pos=atom_pos[res]['tip'], ref_pos=atom_pos[res]['ref'], dist=residue_atoms[res]['dist'])
        structure.export_alignment(out=res+'_align', dummy='Mg')

    # run template conformer search
    search.PARAMS.read_parameters({'NSTEP':100})
    prbs_ = search.TEMPLATE_SEARCH(structure)
    if not prbs_:
        continue
    smart.ExportStructure(prbs_, file.split('.')[0]+'_'+res+'_'+probe_name)
    continue

    # gather parameters
    properties = desc.get_all_properties(prbs_, structure.MOL, id=0, prox_radius=5.0, alpha=0)
    properties['Name'] = name
    properties['Residue'] = res
    if props.empty:
        props = properties
    else:
        props = pd.concat([props, properties], join='inner', ignore_index=True)

props.to_excel('8w8r.xlsx')
print('DONE')
