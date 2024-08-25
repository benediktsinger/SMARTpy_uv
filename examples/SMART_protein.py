import os, sys
import pandas as pd
#MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
from rdkit import Chem

file = '8w8q_O_SiH2_6_linear_docked.pdb'

docked = SMART_conf.read_structure(file)
print(docked)
SMART_conf.PARAMS.read_parameters({'NSTEP':20})
try:
    # Conf search (custom template method)
    cav, struc, cmplx = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked)
except:
    sys.exit()
SMART_conf.save_out(file.split('.')[0]+'_O_SiH2_6_linear')
# Compute properties
properties = SMART_des.get_all_properties(struc, cav, cmplx, 0, prox_radius=5.0, alpha=0)
print(properties)
properties.to_excel('8dsg.xlsx')
