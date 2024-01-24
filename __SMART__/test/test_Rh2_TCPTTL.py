# Test case 2:
#  - read multiple conformations per input file
#  - pass RDKit.Mol object to SMART codes
#  - 2 atoms to define reference plane (tail)
#
import os, sys
import pandas as pd
import numpy as np
import morfeus as mf

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import ForceField
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
PROBES = os.path.join(os.getcwd(), 'Probes')

def Mol2MolSupplier(file=None, sanitize=True):
    mols=[]
    with open(file, 'r') as f:
        line =f.readline()
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol = []
                mol.append(line)
                line = f.readline()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                mol[-1] = mol[-1].rstrip() # removes blank line at file end
                block = ",".join(mol).replace(',','')
                m=Chem.MolFromMol2Block(block,sanitize=sanitize,removeHs=False)
            mols.append(m)
    return(mols)

##########################################################################

dfs = {}
file = 'R-TCPTTL_1'

structure = SMART.ReadFile('mol', file, tip=0, tail=1)
probe = SMART.Probe('acyclic_6SiH2')

print('\nadding probe')
docked = SMART.addProbe(structure, probe, dist=2.0)
SMART.ExportStructure(docked, file+'_probe')

print('\nconformer search')
PAR = {'FIXAT':[-22, docked.GetNumAtoms(), 1], 'EWIN':41.9, 'NSTEP':500, 'MAXROTATION':330,'MINROTATION':30, 'SEED':42}
SMART_conf.PARAMS.read_parameters(PAR)
confs = None
try:
    confs = SMART_conf.start(docked)
except RecursionError:
    print('simulation failed for', file)
    sys.exit()

if confs == None:
    print('simulation failed for',file)
    sys.exit()

SMART_conf.save_out('SMART_confs_'+file)

print('\ncomputing SASA')
properties = SMART_des.RDKit_Properties(confs, 0)

dfs['SASAcavity'] = properties['SASAcavity']
dfs['SASAcat'] = properties['SASAcat']
dfs['SASAcomplx'] = properties['SASAcomplx']
dfs['CSA'] = properties['CSA']
dfs['ESA'] = properties['ESA']
dfs['Asphericity'] = properties['Asphericity']
dfs['Eccentricity'] = properties['Eccentricity']
dfs['InertialShapeFactor'] = properties['InertialShapeFactor']
dfs['SpherocityIndex'] = properties['SpherocityIndex']
dfs['NormalizedInertiaRatio1/3'] = properties['NormalizedInertiaRatio1/3']

print('\ncomputing Sterimol/VBur')
sterimol = SMART_des.DBSTEP_Properties(confs, 0)
dfs['Sterimol-L'] = sterimol['Sterimol-L']
dfs['Sterimol-B1'] = sterimol['Sterimol-B1']
dfs['Sterimol-B5'] = sterimol['Sterimol-B5']
dfs['VBurPROXcavity'] = sterimol['VBurPROXcavity']
dfs['%VBurPROXcavity'] = sterimol['%VBurPROXcavity']
dfs['VBurFULLcavity'] = sterimol['VBurFULLcavity']

df = pd.DataFrame(dfs, index=[0])
with pd.ExcelWriter('SMART_'+file+'.xlsx', engine='openpyxl') as writer:
    df.to_excel(writer)

print()
print('DONE')
