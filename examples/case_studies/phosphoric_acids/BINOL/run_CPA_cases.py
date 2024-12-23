# BINOL/SPINOL-CPA case studies (adapted from disclosed CPA dataset: https://pubs.acs.org/doi/10.1021/acscatal.9b03535)
import os, sys
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

sys.path.append('/Users/rebecca/Desktop/Coding_Projects/SMART-molecular-descriptors/')
import SMART as smart
from SMART import conf_search as search
from SMART import descriptors as desc

props = pd.DataFrame()
for file in os.listdir(os.getcwd()):
    #if file == 'SPINOL_c.xtbopt.xyz' or file == 'SPINOL_i.xtbopt.xyz':
    if file.endswith('xtbopt.xyz'):
        name = file.split('.')[0]
        if os.path.isfile(name+'_S_SiH2_10_cyclic.sdf'):
            continue
        print(name)
        # Prepare input file in RDKit
        mol = Chem.MolFromXYZFile(file)
        rdDetermineBonds.DetermineConnectivity(mol)
        #rdDetermineBonds.DetermineBondOrders(mol, charge=-1)
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol)

        # save xyz of P atom for probe position
        # get xyz of backbone O atoms for vector reference position
        # remove P, P-OH, P=O
        oids, p_pos = [], None
        ids_rmv = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'P':
                p_pos = np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx()))
                for i in atom.GetNeighbors():
                    ids_rmv.append(i.GetIdx())
                    if (len(i.GetNeighbors()) == 1 and ['P'] == [j.GetSymbol() for j in i.GetNeighbors()]) or (len(i.GetNeighbors()) == 2 and 'H' in [j.GetSymbol() for j in i.GetNeighbors()]):
                        ids_rmv.append(i.GetIdx())
                        for j in i.GetNeighbors():
                            if j.GetSymbol() == 'H':
                                ids_rmv.append(j.GetIdx())
                    else:
                        oids.append(np.array(mol.GetConformer().GetAtomPosition(i.GetIdx())))
                ids_rmv.append(atom.GetIdx())
                break

        ids_rmv.sort(reverse=True)
        mol_ = Chem.RWMol(mol)
        for i in ids_rmv:
            mol_.RemoveAtom(i)

        # define structure
        structure = smart.ReadMol(mol_)
        # set probe tether at P xyz coordinates
        # set backbone O atoms as reference atoms
        structure.reference_geometry(tip_pos=p_pos, refs_pos=oids, geom='trigonal', dist=0.1)
        smart.ExportStructure(structure.MOL, outname=name+'_cut')
        #continue
        probe_name = 'S_SiH2_10_cyclic'
        search.PARAMS.read_parameters({'NSTEP':50,'VERBOSE':True})#, 'CLASHTOL':0.05
        #try:
        if search.TEMPLATE.is_template(probe_name):
            # conformer search
            print('get template')
            search.TEMPLATE.GetTemplate(probe_name)
        else:
            # make template + conformer search
            print('generate template')
            search.TEMPLATE.GenerateTemplate(probe_name)
        # run template conformer search
        prbs_ = search.TEMPLATE_SEARCH(structure)
        if not prbs_:
            continue
        smart.ExportStructure(prbs_, outname=name+'_'+probe_name)
        continue
    	# gather all parameters
        properties = desc.get_all_properties(prbs_, structure.MOL, id=0, prox_radius=5.0, alpha=0)
        #p = desc.Morfeus_QO_Properties(structure.MOL, prbs_, id=0, z_axis=bkbone, xz_plane=bkbone, prox_radius=5.0)
        properties['Name'] = name
        properties['Probe'] = probe_name
        if props.empty:
            props = properties
        else:
            props = pd.concat([props, properties], join='inner', ignore_index=True)
        #except Exception as e:
        #    print(e)
        #    continue

#props.to_excel('BINOL.xlsx')
print('DONE')
