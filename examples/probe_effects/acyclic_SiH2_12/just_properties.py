import os, sys
import pandas as pd
MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
sys.path.append(MODS)
import SMART_descriptors as SMART_des

props = pd.DataFrame()
for file in os.listdir(os.getcwd()):
    if file.endswith('.mol2'):
        print(file)
        for i in ['0', '1']:
            try:
                name = file.split('.')[0]
                cav_name = name+'_S_SiH2_12_acyclic_Rh'+i+'_cavity.sdf'

                print('reading', cav_name, file)
                probe = SMART_des.read_CAVITY_file(cav_name)
                struc = SMART_des.read_STRUCTURE_file(file)

                z = 0
                xy = []
                for atom in struc.GetAtoms():
                    if atom.GetIdx() != int(i):
                        continue
                    for n in atom.GetNeighbors():
                        if n.GetSymbol() == 'Rh':
                            z = n.GetIdx()
                        elif n.GetSymbol() == 'O':
                            xy.append(n.GetIdx())
                print('getting props')
                properties = SMART_des.get_all_properties(struc, probe, int(i), prox_radius=5.0, alpha=0)
                properties['Name'] = name
                p = SMART_des.Morfeus_Octant_Properties(struc, probe, int(i), z_axis=[z], xz_plane=xy, prox_radius=5.0)
                properties = pd.DataFrame(p, index=[0]).merge(pd.DataFrame(properties, index=[0]), how='inner', left_index=True, right_index=True)
                if props.empty:
                    props = properties
                else:
                    props = pd.concat([props, properties], join='inner', ignore_index=True)
            except:
                continue
props.to_excel('S_SiH2_12_cyclic_OUT.xlsx')
