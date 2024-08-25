import os, sys
import pandas as pd
#MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_descriptors as SMART_des
#from rdkit import Chem

probe = SMART_des.read_CAVITY_file('8dsg_O_SiH2_6_linear_cavity.sdf')
struc = SMART_des.read_STRUCTURE_file('8dsg_trunc.pdb')

properties = SMART_des.get_all_properties(struc, probe, 0, prox_radius=5.0, alpha=0)
#props = SMART_des.Morfeus_Octant_Properties(struc, probe, id, z_axis=[3225], xz_plane=[3717,3719,3720], prox_radius=5.0)
#properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
#print(properties)
properties.to_excel('protein.xlsx')
