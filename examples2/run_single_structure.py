# Code template for computing SMART descriptors for a single file passed as a cmd line argument
# Run -
#   python3 run_single_structure.py [structure.ext] [probe.mol2]
#

import os, sys
import pandas as pd

import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
PROBES = os.path.join(os.getcwd(), 'Probes')

structure_file = argv[1] # structure.mol2
probe_file = argv[2] # probe.mol2
binding_atom_idx = 0
ref_atom_idx = 1

# Initialize structure
structure = SMART.ReadFile(structure_file.split('.')[1], structure_file.split('.')[0], tip=binding_atom_idx, tail=ref_atom_idx)

# Initialize probe
probe = SMART.Probe(probe_file.split('.')[0])

# Dock probe to structure
print('\nadding probe',probe)
docked = SMART.addProbe(structure, probe, dist=2.0)

# Export docked structure (Optional)
SMART.ExportStructure(docked, file.split('.')[0]+'_'+probefile.split('.')[0])

# Initialize conf search parameters
print('\nconformer search')
PAR = {'FIXAT':[-probe.NATOMS, docked.GetNumAtoms(), 1], 'NSTEP':500, 'MAXROTATION':330,'MINROTATION':30, 'SEED':42}
SMART_conf.PARAMS.read_parameters(PAR)

confs = None
try:
    # Run conformer search
    confs = SMART_conf.start(docked)
except RecursionError:
    print('simulation failed for', file.split('.')[0],probefile.split('.')[0])
    sys.exit()

# Export conformer ensemble (Optional)
SMART_conf.save_out('SMART_confs_'+file.split('.')[0]+'_'+probefile.split('.')[0])

# Molecular descriptors
print('\ncomputing SASA')
rdkit_properties = SMART_des.RDKit_Properties(confs, binding_atom_idx)
print('\ncomputing Sterimol/VBur')
sterimol_properties = SMART_des.DBSTEP_Properties(confs, binding_atom_idx)
print('\ncomputing Delauney triangulation')
pyvista_properties = SMART_des.PyVista_Properties(confs, binding_atom_idx)

# Export descriptors
df = pd.concat([pd.DataFrame(rdkit_properties), pd.DataFrame(sterimol_properties), pd.DataFrame(pyvista_properties)])
with pd.ExcelWriter('SMART_'+file+'_'+probefile.split('.')[0]+'.xlsx', engine='openpyxl') as writer:
    df.to_excel(writer)

print('DONE')
