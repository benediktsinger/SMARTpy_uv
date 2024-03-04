import os, sys
import pandas as pd
MODS = os.path.join(os.getcwd(), 'modules')
sys.path.append(MODS)
import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
#PROBES = os.path.join(os.getcwd(), 'Probes')

#pt = Chem.GetPeriodicTable()

file = 'R-TCPTTL_1.mol2'
# Prepare structure
structure = SMART.ReadFile(file)
print(structure.NumAtoms)
structure.Reference_Vector(tip=0, tail=1, dist=2.0)
print(structure.Vector)
sys.exit()
# Prepare probe
probe = SMART.Probe('S_SiF2_8_cyclic')
probe_count = probe.NATOMS
# Dock probe to structure
print('\nadding probe',probe)
docked = SMART.addProbe(structure, probe, dist=2.0)
SMART.ExportStructure(docked, file.split('.')[0]+'_S_SiF2_8_cyclic')

print('DONE')
