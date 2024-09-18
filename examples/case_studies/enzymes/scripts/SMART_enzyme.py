import os, sys
import pandas as pd
import numpy as np

sys.path.append('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/')
import SMART as smart
from SMART import conf_search as search
from SMART import descriptors as desc
import rdkit
from rdkit import Chem

props = pd.DataFrame()
#files = {'5ucw_trunc.pdb':[3680,3198], '1fag_trunc.pdb':[3712,3225]}
files = {'8dsg_trunc.pdb':[3671,3183]}
for file in files.keys():
    path = os.path.join('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/examples/case_studies/enzymes/', file)
    mol = Chem.MolFromPDBFile(path, sanitize=False, proximityBonding=True)
    binding_atom, ref_atom = files[file][0], files[file][1]

    # Prepare structure
    structure = smart.ReadFile(path)
    # Set reference vector
    structure.reference_vector(tip_id=binding_atom, ref_id=ref_atom, dist=2.0)
    structure.export_alignment(out=file.split('.')[0]+'_align', dummy='Mg')

    # Prepare probe
    probe_name = 'O_SiH2_6_linear'#O_SiH2_6_linear
    if search.TEMPLATE.is_template(probe_name):
        # conformer search
        print('get template')
        search.TEMPLATE.GetTemplate(probe_name)
    else:
        # make template + conformer search
        print('gen template')
        search.TEMPLATE.GenerateTemplate(probe_name)

    search.PARAMS.read_parameters({'NSTEP':50, 'VERBOSE':True})
    prbs_ = search.TEMPLATE_SEARCH(structure)
    if not prbs_:
        sys.exit()
    smart.ExportStructure(prbs_, file.split('.')[0]+'_'+probe_name)
    continue

    # gather parameters
    properties = desc.get_all_properties(prbs_, structure.MOL, id=0, prox_radius=4.0, alpha=0)
    properties['Name'] = name
    properties['Enzyme'] = file.split('.')[0]

    if props.empty:
        props = properties
    else:
        props = pd.concat([props, properties], join='inner', ignore_index=True)

#props.to_excel('enzymes.xlsx')
print('DONE')
