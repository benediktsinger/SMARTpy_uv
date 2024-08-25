import os, sys
import pandas as pd
#MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
from rdkit import Chem

file = '8w8q.pdb'
print(Chem.MolFromMol2File(file, sanitize=False))
binding_atom = 675-2
ref_atom = 673-2

# Prepare structure
structure = SMART.ReadFile(file)
structure.reference_vector(binding_atom, ref_atom, dist=1.0) #tip=0, tail =1

#structure.reference_geometry(binding_atom, 'trigonal', ref_atoms, dist=2.0)

# Prepare probe
probe = SMART.ReadProbe('O_SiH2_6_linear.mol2') #FIX PROBE ATOM NUMS

# Adjust clash tolerance (default: 0.0)
#SMART.PARAMS.CLASHTOL = 0.05 # 0.05 Å further
#SMART.PARAMS.ClashTol = -0.05 # 0.05 Å closer

# Dock probe to structure
print('\nadding probe O_SiH2_6_linear.mol2 to '+file)
docked = SMART.add_probe(structure, probe)
SMART.ExportStructure(docked, outname=file.split('.')[0]+'_O_SiH2_6_linear_docked')

SMART_conf.PARAMS.read_parameters({'NSTEP':20})
try:
    # Conf search (custom template method)
    cav, struc, cmplx = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked)
except:
    sys.exit()
SMART_conf.save_out(file.split('.')[0]+'_O_SiH2_6_linear')
# Compute properties
properties = SMART_des.get_all_properties(structure.MOL, cav, binding_atom, prox_radius=3.5, alpha=0)
print(properties)
properties.to_excel('8w8q.xlsx')
