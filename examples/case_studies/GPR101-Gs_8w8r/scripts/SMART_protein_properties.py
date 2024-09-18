import os, sys
import pandas as pd
import numpy as np

sys.path.append('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/')
import SMART as smart
from SMART import conf_search as search
from SMART import descriptors as desc
import rdkit
from rdkit import Chem

file = '8w8r_trunc.pdb'
structure = smart.ReadFile(file)

props = pd.DataFrame()
residue_atoms = {'T111':
                    {'tip':[646,647], 'ref':645, 'dist':2.0},
                'L181':
                    {'tip':[1205,1204], 'ref':1203, 'dist':2.0}}
for res in residue_atoms.keys():
    prbs_ = desc.read_CAVITY_file(file.split('.')[0]+'_'+res+'_O_SiH2_6_linear.sdf')
    # gather parameters
    #properties = desc.PyVista_Properties(structure.MOL, prbs_, id=0, prox_radius=3.5, a=0)
    properties = desc.Morfeus_Properties(structure.MOL, prbs_, id=residue_atoms[res]['ref']+1, prox_radius=3.5, vol=True, sasa=False, sterimol=True)
    properties['Res'] = res

    if props.empty:
        props = pd.DataFrame(properties, index=[0])
    else:
        props = pd.concat([props, pd.DataFrame(properties, index=[0])], join='inner', ignore_index=True)

props.to_excel('GPR101_.xlsx')
print('DONE')
