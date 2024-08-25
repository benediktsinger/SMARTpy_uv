import os, sys
import pandas as pd
#MODS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'modules')
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
from rdkit import Chem

file = '8dsg.pdb'
binding_atom = 0
ref_atom = 0

# Prepare structure
structure = SMART.ReadFile(file)

for atom in structure.MOL.GetAtoms():
    if atom.GetSymbol() == 'Fe':
        binding_atom = atom.GetIdx()
        for n in atom.GetNeighbors():
            if n.GetSymbol() == 'O':
                ref_atom = n.GetIdx()
                break
        break

# Set reference vector
structure.reference_vector(binding_atom, ref_atom, dist=1.0) #tip=0, tail =1

# Prepare probe
probe = SMART.ReadProbe('S_SiH2_8_cyclic.mol2')
# Adjust clash tolerance (default: 0.0)
#SMART.PARAMS.CLASHTOL = 0.05 # 0.05 Å further
#SMART.PARAMS.ClashTol = -0.05 # 0.05 Å closer

# Dock probe to structure
print('\nadding probe S_SiH2_8_cyclic.mol2 to '+file)
docked = SMART.add_probe(structure, probe)
SMART.ExportStructure(docked, outname=file.split('.')[0]+'_S_SiH2_8_cyclic_docked')

# Set template steps to 20
SMART_conf.PARAMS.read_parameters({'NSTEP':20})
try:
    # Conf search (custom template method)
    cav, struc, cmplx = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked)
except:
    sys.exit()
SMART_conf.save_out(file.split('.')[0]+'_S_SiH2_8_cyclic')
# Compute properties
properties = SMART_des.get_all_properties(structure.MOL, cav, binding_atom, prox_radius=3.5, alpha=0)
print(properties)
properties.to_excel('8dsg.xlsx')
