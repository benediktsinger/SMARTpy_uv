import os, sys
import pandas as pd

sys.path.append('/uufs/chpc.utah.edu/common/home/u1319506/SMART_manuscript/')
import SMART as smart
from SMART import conf_search as search
from SMART import descriptors as desc

PROBES = os.path.join(os.getcwd(), 'Probes')
templates = os.path.join(os.getcwd(), 'templates')
#pt = Chem.GetPeriodicTable()

dfs = {}
for file in os.listdir(os.getcwd()):
    if file.endswith('.mol2'):
        print('reading',file)
        for i in [[0,1],[1,0]]:
            # Prepare structure
            name = file.split('.')[0]
            structure = smart.ReadFile(name+'.mol2')
            structure.reference_vector(tip_id=i[0], ref_id=i[1], dist=2.0)
            # Prepare probe
            probe_name= 'S_SiH2_12_cyclic'
            if search.TEMPLATE.is_template(probe_name):
                print('get template')
                search.TEMPLATE.GetTemplate(probe_name)
            else:
                # make template + conformer search
                print('generate template')
                search.TEMPLATE.GenerateTemplate(probe_name)

            print('\nconformer search')
            search.PARAMS.read_parameters({'NSTEP':50})

            prbs_ = search.TEMPLATE_SEARCH(structure)
            if not prbs_:
                continue
            smart.ExportStructure(prbs_, outname=name+'_'+probe_name)
            #continue


print('DONE')
