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
files = {'8dsg_trunc.pdb':[3671,3183], '5ucw_trunc.pdb':[3680,3198]}
for file in files.keys():
    structure = smart.ReadFile(file)
    prbs_ = desc.read_CAVITY_file(file.split('.')[0]+'_O_SiH2_6_linear.sdf')
    # gather parameters
    #properties = desc.PyVista_Properties(structure.MOL, prbs_, id=0, prox_radius=3.5, a=0)
    properties = desc.Morfeus_Properties(structure.MOL, prbs_, id=files[file][0], prox_radius=3.5, vol=True, sasa=False, sterimol=True)
    #properties = desc.get_all_properties(prbs_, structure.MOL, id=files[file][0], prox_radius=4.0, alpha=0)
    properties['Enzyme'] = file.split('.')[0]

    if props.empty:
        props = pd.DataFrame(properties, index=[0])
    else:
        props = pd.concat([props, pd.DataFrame(properties, index=[0])], join='inner', ignore_index=True)

props.to_excel('enzymes.xlsx')
print('DONE')
