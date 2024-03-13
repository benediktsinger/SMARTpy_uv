import os, sys
import pandas as pd
MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
sys.path.append(MODS)
import SMART_descriptors as SMART_des

probe, struc, cmplx = SMART_des.read_files('R-TCPTTL_1_S_SiF2_8_cyclic')
properties = SMART_des.get_all_properties(struc, probe, cmplx, 0, prox_radius=5.0, alpha=0)
print(properties)
properties.to_excel('OUT.xlsx')
