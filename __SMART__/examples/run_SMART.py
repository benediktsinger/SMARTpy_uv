import os, sys
import pandas as pd

import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
import SMART_descriptors as SMART_des
PROBES = os.path.join(os.getcwd(), 'Probes')
#pt = Chem.GetPeriodicTable()

dfs = {}
for file in os.listdir(os.getcwd()):
    if file.endswith('.sdf'):
        print('reading',file)
        # Prepare structure
        structure = SMART.ReadFile('mol2', file.split('.')[0], tip=0, tail=1)
        # Prepare probe
        for probefile in os.listdir(PROBES):
            if 'acyclic' not in probefile:
                continue

            if os.path.isfile(os.path.join(os.getcwd(),'SMART_'+file+'_'+probefile.split('.')[0]+'.xlsx')):
                continue

            probe = SMART.Probe(probefile.split('.')[0])
            probe_count = probe.NATOMS
            # Dock probe to structure
            print('\nadding probe',probe)
            docked = SMART.addProbe(structure, probe, dist=2.0)
            SMART.ExportStructure(docked, file.split('.')[0]+'_'+probefile.split('.')[0])

            print('\nconformer search')
            PAR = {'FIXAT':[-probe_count, docked.GetNumAtoms(), 1], 'NSTEP':500, 'MAXROTATION':330,'MINROTATION':30, 'SEED':42} #'EWIN':41.9
            SMART_conf.PARAMS.read_parameters(PAR)
            confs = None
            try:
                confs = SMART_conf.start(docked)
            except RecursionError:
                print('simulation failed for', file.split('.')[0],probefile.split('.')[0])
                continue
            SMART_conf.save_out('SMART_confs_'+file.split('.')[0]+'_'+probefile.split('.')[0])
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
            with pd.ExcelWriter('SMART_'+file+'_'+probefile.split('.')[0]+'.xlsx', engine='openpyxl') as writer:
                df.to_excel(writer)

            print()
            print('DONE')



print('DONE')
