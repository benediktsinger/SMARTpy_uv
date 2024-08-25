# to run:
# /path/to/chimera --nogui --nostatus --script "SMART_chimera_descriptors.py <lvl_cav> <lvl_mol> <rad>" > SMART.txt
#
import os, sys, argparse
from chimera import runCommand as rc

def get_chimera_descriptors(filename, path=os.getcwd(), radius=3.5, lvl_mol=0.07, lvl_cav=4):
    cavity_file = filename + '_cavity.sdf'
    mol_file = filename + '_mol.sdf'

    cavity_path = os.path.join(path, cavity_file)
    mol_path = os.path.join(path, mol_file)

    try:
        rc("open "+cavity_path) #0
        rc("shape sphere center @S1 radius "+str(rad)) #1
        rc("select #0@S1") # remove tether atom, (S1)
        rc("delete sel")
        rc("molmap #0 2.8") #2
        rc("volume #2 level "+str(lvl_cav))

        rc("open "+mol_path) #3
        rc("molmap #3 2.8") #3.1
        rc("volume #3.1 level "+str(lvl_mol))

        rc("measure volume #2")
        rc("measure area #2")
        rc("measure contactArea #2 #3.1 1.1 color black offset 0")

        rc("mask #2 #1") #4
        rc("measure volume #4")
        rc("measure area #4")
        rc("measure contactArea #4 #3.1 1.1 color black offset 0")

        rc("close all")
        rc("close session")
    except:
        print('failed for '+filename)
        rc("close all")
        rc("close session")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe Conformer-Generation Utility Package.',epilog='Probe conformer search by torsional/template search')
    parser.add_argument('-f', required=True) #filename (filename_cavity.sdf, filename_mol.sdf, filename_cplx.sdf)
    parser.add_argument('-r', default=3.5, required=False) #proximal radius (Å)
    parser.add_argument('-path', default=os.getcwd(), required=False)
    parser.add_argument('-lvl_cav', default=4, required=True) #volume level cavity
    parser.add_argument('-lvl_mol', default=0.07, required=True) #volume level mol

    args = parser.parse_args()

    get_chimera_descriptors(args.f, args.r, args.path, args.lvl_mol, args.lvl_cav)
