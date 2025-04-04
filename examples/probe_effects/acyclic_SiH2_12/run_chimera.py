import os
import SMART_chimera_descriptors as SMART

for file in os.listdir(os.getcwd()):
    if file.endswith('.mol2'):
        for i in ['0', '1']:
            try:
                cavity_file = file.split('.')[0]+'_S_SiH2_12_acyclic_Rh'+i+'_cavity.sdf'
                SMART.get_chimera_descriptors(file, cavity_file, path=os.getcwd(), radius=5.0, lvl_mol=0.07, lvl_cav=4)
            except:
                continue
