import os, sys
import pandas as pd

sys.path.append('/Users/rebecca/Documents/SMART-molecular-descriptors/')
import SMART as smart
from SMART import descriptors as desc

pyv_props = pd.DataFrame()
mf_props = pd.DataFrame()
oct_props = pd.DataFrame()
for file in os.listdir(os.getcwd()):
    if file.endswith('.mol2'):
        name = file.split('.')[0]
        cav_name = name+'_S_SiH2_12_cyclic.sdf'
        cav_file = os.path.join(os.getcwd(), cav_name)
        if not os.path.isfile(cav_file):
            continue
        #print('reading', cav_name, file)
        probe = desc.read_CAVITY_file(cav_file)
        struc = smart.ReadFile(file)
        #print(probe, struc)

        xy = []
        for atom in struc.MOL.GetAtoms():
            if atom.GetIdx() != 1:
                continue
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'Rh':
                    continue
                elif n.GetSymbol() == 'O':
                    xy.append(n.GetIdx())

        pyv_properties = desc.get_Cloud_Properties(struc, probe, id=1, prox_radius=5.0, alpha=0)
        pyv_properties['Name'] = name
        #mf_properties = desc.get_BuriedVolume_Properties(struc, probe, id=1, prox_radius=5.0, vol=True, sasa=False, sterimol=False)
        #mf_properties['Name'] = name
        #oct_properties = get_Octant_Properties(struc, probe, id=1, z_axis=[1], xz_plane=xy, prox_radius=5.0, octant=True, quadrant=False)
        #oct_properties['Name'] = name
        if pyv_props.empty:
            pyv_props = pd.DataFrame(pyv_properties, index=[0])
        else:
            pyv_props = pd.concat([pyv_props, pd.DataFrame(pyv_properties, index=[0])], join='inner', ignore_index=True)

        #if mf_props.empty:
        #    mf_props = mf_properties
        #else:
        #    mf_props = pd.concat([mf_props, mf_properties], join='inner', ignore_index=True)

        #if oct_props.empty:
        #    oct_props = oct_properties
        #else:
        #    oct_props = pd.concat([oct_props, oct_properties], join='inner', ignore_index=True)

pyv_props.to_excel('SMART_pyvista.xlsx')
#mf_props.to_excel('SMART_morfeus.xlsx')
#oct_props.to_excel('SMART_morfeus_octant.xlsx')
