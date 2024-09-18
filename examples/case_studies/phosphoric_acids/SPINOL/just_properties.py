
# BINOL/SPINOL-CPA case studies (adapted from disclosed CPA dataset: https://pubs.acs.org/doi/10.1021/acscatal.9b03535)
import os, sys
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

sys.path.append('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/')
import SMART as smart
from SMART import descriptors as desc

props = pd.DataFrame({})
for file in os.listdir(os.getcwd()):
    if file.endswith('xtbopt.xyz'):
    #if file == 'SPINOL_c.xtbopt.xyz' or file == 'SPINOL_i.xtbopt.xyz':
        name = file.split('.')[0]
        if not os.path.isfile(name+'_cut.sdf'):
            continue
        print(name)
        structure = smart.ReadFile(name+'_cut.sdf')
        if os.path.isfile(name+'_S_SiH2_10_cyclic.sdf'):
            prbs_ = desc.read_CAVITY_file(name+'_S_SiH2_10_cyclic.sdf')
            #properties = desc.get_all_properties(prbs_, structure.MOL, id=0, prox_radius=5.0, alpha=0)
            properties = desc.Morfeus_QO_Properties(structure.MOL, prbs_, id=0, z_axis=[0], xz_plane=[13,14], quadrant=True, octant=True)
            properties['Name'] = name
            properties['Probe'] = 'S_SiH2_10_cyclic'
            properties = pd.DataFrame(properties, index=[0])
            if props.empty:
                props = properties
            else:
                props = pd.concat([props, properties], join='inner', ignore_index=True)

props.to_excel('SPINOL_OCT.xlsx')
print('DONE')
