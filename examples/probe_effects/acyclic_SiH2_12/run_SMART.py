import os, sys
import pandas as pd

import SMART_probe_utils as SMART
import SMART_conf_search as SMART_conf
#import SMART_descriptors as SMART_des
PROBES = os.path.join(os.getcwd(), 'Probes')
templates = os.path.join(os.getcwd(), 'templates')
#pt = Chem.GetPeriodicTable()

dfs = {}
for file in os.listdir(os.getcwd()):
    if file.endswith('.mol2'):
        print('reading',file)
        for i in [[0,1],[1,0]]:
            # Prepare structure
            fname = file.split('.')[0]
            structure = SMART.ReadFile(fname+'.mol2')
            structure.reference_vector(i[0], i[1], dist=2.0)
            # Prepare probe
            pname= 'S_SiH2_12_cyclic'
            probe = SMART.ReadProbe(pname+'.mol2', path=os.path.join(os.getcwd(), 'Probes'))
            #probe_count = probe.NATOMS
            # Dock probe to structure
            print('\nadding probe',pname,'to',fname)
            docked = SMART.add_probe(structure, probe)
            #SMART.ExportStructure(docked, file.split('.')[0]+'_'+probefile.split('.')[0])
            print('\nconformer search')
            #PAR = {'FIXAT':[-probe_count, docked.GetNumAtoms(), 1], 'NSTEP':500, 'MAXROTATION':330,'MINROTATION':30, 'SEED':42} #'EWIN':41.9
            #SMART_conf.PARAMS.read_parameters(PAR)
            #SMART_conf.TEMPLATE.get_template('S_SiF2_8_cyclic.sdf')
            confs = None
            try:
                probes, mol, cplx  = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked)
                #SMART_conf.save_out(fname+'_'+pname+'_Rh'+str(i[0]))
            except RecursionError:
                print('simulation failed for', fname,pname)
                continue

            if probes == None:
                continue
            SMART_conf.save_out(fname+'_'+pname+'_Rh'+str(i[0]))
            #SMART_conf.save_out(fname+'_'+pname+'_Rh'+str(i[0]))

print('DONE')
