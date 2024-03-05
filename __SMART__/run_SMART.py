import os, sys
import pandas as pd
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des

file = 'R-TCPTTL_1.mol2'

# Prepare structure
structure = SMART.ReadFile(file)
#print(structure.NumAtoms)
structure.reference_vector(0, 1, dist=2.0)

# Prepare probe
probe = SMART.ReadProbe('S_SiF2_12_cyclic.mol2')
#probe_count = probe.NumAtoms

# Dock probe to structure
print('\nadding probe',probe)
docked = SMART.add_probe(structure, probe)

# Conf search (custom template method)
cav, struc, cmplx = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked)
#SMART.ExportStructure(docked, file.split('.')[0]+'_S_SiF2_8_cyclic')

# Compute properties
properties = SMART_des.get_all_properties(struc, cav, id, prox_radius=3.5, alpha=0)
print(properties)
print('DONE')
