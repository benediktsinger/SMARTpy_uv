#!/usr/bin/python
# SMART probe conformational search utility
#   version: b2.1 (released 03-14-2024)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#
# Dev Notes:
#  Features of FullMonte inspired development of this code.
#  - FullMonte is an archived public GitHub repository. (https://github.com/patonlab/FullMonte)
#  - Robert Paton (https://github.com/patonlab/)
#
############# ------- DEPENDENCIES
import sys, os, random, math, argparse, time, mathutils
sys.setrecursionlimit(10000)
import numpy as np
#from multiprocessing import Pool
from scipy.spatial.transform import Rotation as R
#from scipy.stats import ortho_group
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit import ForceField
from rdkit import RDLogger
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__),'templates')
PROBE_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Probes')
############# ------- UTILITY CLASSES
class PARAMS:
    def __init__(self):
        pass
    def read_parameters(PAR):
        # read parameter specification
        ERRTAGS = []
        for i in PAR.keys():
            if i not in dir(PARAMS):
                ERRTAGS.append(i)
                continue
            else:
                if i == 'SEED':
                    #PARAMS.SEED = int(PAR[i])
                    print('SEED not yet supported!')
                elif i == 'NSTEP':
                    PARAMS.NSTEP = int(PAR[i])
                elif i == 'MAXROTATION':
                    PARAMS.MAXROTATION = float(PAR[i])
                elif i == 'MINROTATION':
                    PARAMS.MINROTATION = float(PAR[i])
                elif i == 'BINDING_POS':
                    PARAMS.BINDING_POS = np.array(PAR[i])
                elif i == 'VERBOSE':
                    PARAMS.VERBOSE = PAR[i]
                elif i == 'NPROCS':
                    PARAMS.NPROCS = int(PAR[i])
        return ERRTAGS
    def read_parameter_file(Pin):
        # open file
        PAR = {}
        with open(Pin, 'r') as param:
            for line in param:
                if len(line.split()) > 2:
                    PAR[line.split()[0]] = line.split()[1:]
                else:
                    PAR[line.split()[0]] = line.split()[1]
        PARAMS.read_parameters(PAR)
PARAMS.MAXROTATION = 330
PARAMS.MINROTATION = 30
PARAMS.NSTEP = 10
PARAMS.SEED = None
PARAMS.CLASHTOL = 0.0
PARAMS.BINDING_POS = None
PARAMS.VERBOSE = True
PARAMS.NPROCS = 4

class TEMPLATE():
    def __init__(self):
        pass
    def is_template(template):
        global TEMPLATE_PATH
        if os.path.isfile(os.path.join(TEMPLATE_PATH, template+'.sdf')):
            return True
        else:
            return False
    def GenerateTemplate(template, numConfs=200, pruneRmsThresh=0.1):
        global TEMPLATE_PATH, PROBE_PATH
        try:
            starttime = time.time()
            # get probe MOL
            TEMPLATE.NAME = template
            TEMPLATE.MOLS = read_structure(os.path.join(PROBE_PATH, template+'.mol2'))
            # generate template conformers
            AllChem.EmbedMultipleConfs(TEMPLATE.MOLS, numConfs=numConfs, pruneRmsThresh=pruneRmsThresh, useRandomCoords=True, useMacrocycleTorsions=True, useSmallRingTorsions=True)
            TEMPLATE.NCONFS = TEMPLATE.MOLS.GetNumConformers()
            # save template to sdf file
            writer = Chem.SDWriter(os.path.join(TEMPLATE_PATH, template+'.sdf'))
            print('writing template to file: ', template+'.sdf')
            for cid in [conf.GetId() for conf in TEMPLATE.MOLS.GetConformers()]:
                writer.write(TEMPLATE.MOLS, confId=cid)
            print('***************** Done *****************')
            print("--- Time Elapsed: %s seconds ---" % (time.time() - starttime))
        except Exception as e:
            print('Invalid probe structure')
            print(e)
            sys.exit()
    def GetTemplate(template):
        global TEMPLATE_PATH
        try:
            # get template ENS
            TEMPLATE.NAME = template
            TEMPLATE.MOLS = None
            for mol in Chem.SDMolSupplier(os.path.join(TEMPLATE_PATH, template+'.sdf'), removeHs=False):
                if not TEMPLATE.MOLS:
                    TEMPLATE.MOLS = mol
                else:
                    TEMPLATE.MOLS.AddConformer(mol.GetConformer(), assignId=True)
            TEMPLATE.NCONFS = TEMPLATE.MOLS.GetNumConformers()
        except Exception as e:
            print('Invalid template structure')
            print(e)
            sys.exit()
TEMPLATE.MOLS = None
TEMPLATE.NAME = None
TEMPLATE.NCONFS = None

class MOL_INIT:
    pass
MOL_INIT.FNAME = None
MOL_INIT.TYPE = None
MOL_INIT.INIT_MOL = None
MOL_INIT.SAVE_CONFS = None

############# ------- MAIN FNXS
def save_out(name, mol=None):
    writer = Chem.SDWriter(name+'.sdf')
    print('writing conformers to file',name+'.sdf')
    if mol:
        m = mol
    else:
        m = MOL_INIT.SAVE_CONFS

    for cid in [conf.GetId() for conf in m.GetConformers()]:#range(1, CONFS.PROBES.GetNumConformers()):
        writer.write(m, confId=cid)

def ConfToMol(mol, conf):
    id = conf.GetId()
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(mol.GetConformer(id))
    return new_mol

def read_structure(Fin):
    # open file and read structure
    MOL_INIT.FNAME = Fin.split('.')[0]
    MOL_INIT.TYPE = Fin.split('.')[1]
    if MOL_INIT.TYPE == 'mol2':
        MOL_INIT.MOL = Chem.MolFromMol2File(Fin, removeHs=False)
    elif MOL_INIT.TYPE == 'mol':
        MOL_INIT.MOL = Chem.MolFromMolFile(Fin, removeHs=False)
    elif MOL_INIT.TYPE == 'xyz':
        MOL_INIT.MOL = Chem.MolFromXYZFile(Fin, removeHs=False)
    elif MOL_INIT.TYPE == 'pdb':
        MOL_INIT.MOL = Chem.MolFromPDBFile(Fin, removeHs=False)
    elif MOL_INIT.TYPE == 'sdf':
        MOL_INIT.MOL = Chem.SDMolSupplier(Fin, removeHs=False)[0]
    return MOL_INIT.MOL

def clash_check(mol, probe, cid):

    for p in range(1, probe.GetNumAtoms()):
        for s in range(mol.GetNumAtoms()):
                i_ = probe.GetConformer(cid).GetAtomPosition(p)
                i_r = Chem.GetPeriodicTable().GetRvdw(probe.GetAtomWithIdx(p).GetAtomicNum())
                j_ = mol.GetConformer().GetAtomPosition(s)
                j_r = Chem.GetPeriodicTable().GetRvdw(mol.GetAtomWithIdx(s).GetAtomicNum())
                d = (i_-j_).LengthSq()
                r_d = i_r + j_r + PARAMS.CLASHTOL
                if (d - r_d) <= 0:
                    return True
    return False

def vectorize(frag, tip, atoms, rotvec, theta):
    # rotate 3D about vector
    vecs = []
    for i in atoms:
        a = np.array(frag.GetConformer().GetAtomPosition(i.GetIdx()))
        dotproduct = rotvec[0]*(float(a[0]) - float(tip[0])) + rotvec[1]*(float(a[1]) - float(tip[1])) + rotvec[2]*(float(a[2]) - float(tip[2]))
        centre = [float(tip[0]) + dotproduct*rotvec[0], float(tip[1]) + dotproduct*rotvec[1], float(tip[2]) + dotproduct*rotvec[2]]
        v = [float(a[0]) - centre[0], float(a[1]) - centre[1], float(a[2]) - centre[2]]
        d = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        px = v[0]*math.cos(theta) + v[1]*math.sin(theta)*rotvec[2] - v[2]*math.sin(theta)*rotvec[1]
        py = v[1]*math.cos(theta) + v[2]*math.sin(theta)*rotvec[0] - v[0]*math.sin(theta)*rotvec[2]
        pz = v[2]*math.cos(theta) + v[0]*math.sin(theta)*rotvec[1] - v[1]*math.sin(theta)*rotvec[0]
        newv = [px + centre[0], py + centre[1], pz + centre[2]]

        vecs.append(newv)
    return vecs

def rotation(theta):
    # rotate 3D about vector
    cid_ = TEMPLATE.MOLS.GetConformers()[0].GetId()
    rotvec = TEMPLATE.MOLS.GetConformer(cid_).GetAtomPosition(0) - TEMPLATE.MOLS.GetConformer(cid_).GetAtomPosition(1) / np.linalg.norm(TEMPLATE.MOLS.GetConformer(cid_).GetAtomPosition(0) - TEMPLATE.MOLS.GetConformer(cid_).GetAtomPosition(1))
    quat = mathutils.Quaternion(rotvec, theta)
    r = R.from_quat(quat)
    r.as_quat()

    for cid in [conf.GetId() for conf in TEMPLATE.MOLS.GetConformers()]:
        newv = r.apply([np.array(TEMPLATE.MOLS.GetConformer(cid).GetAtomPosition(i)) for i in range(TEMPLATE.MOLS.GetNumAtoms())])
        for i in range(TEMPLATE.MOLS.GetNumAtoms()):
            TEMPLATE.MOLS.GetConformer(cid).SetAtomPosition(i, newv[i])
    #save_out(str(math.degrees(theta)), TEMPLATE.MOLS)

def template_search(step=1):
    if PARAMS.VERBOSE:
        print('step ',step,' / ',PARAMS.NSTEP)
    RANDOM_ROT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    for conf in TEMPLATE.MOLS.GetConformers():
        translation_v =  PARAMS.BINDING_POS - np.array(conf.GetAtomPosition(0))
        for i in range(TEMPLATE.MOLS.GetNumAtoms()):
            pt = np.array(conf.GetAtomPosition(i))
            translate = pt+translation_v
            conf.SetAtomPosition(i, translate)
        # check for probe-structure clashes
        if not clash_check(MOL_INIT.MOL, TEMPLATE.MOLS, conf.GetId()):
            if not MOL_INIT.SAVE_CONFS:
                MOL_INIT.SAVE_CONFS = ConfToMol(TEMPLATE.MOLS, conf)
                if PARAMS.VERBOSE:
                    print('-saving conformer-', 1)
            else:
                MOL_INIT.SAVE_CONFS.AddConformer(conf, assignId=True)
                if PARAMS.VERBOSE:
                    print('-saving conformer-', MOL_INIT.SAVE_CONFS.GetNumConformers())

    # rotate template for next iteration
    rotation(math.radians(RANDOM_ROT))
    # end conditions
    if step < PARAMS.NSTEP:
        step +=1
        return template_search(step)
    return

############# ------- HANDLER FNXS
def TEMPLATE_SEARCH(mol=None):
    '''initiate conformational search by ensemble template'''
    if not mol:
        print('No docked MOL specified')
        sys.exit()

    if mol:
        MOL_INIT.MOL = mol
        MOL_INIT.SAVE_CONFS = None
    else:
        print('no input structure!')
        sys.exit()

    # print search parameters
    if PARAMS.VERBOSE:
        print('!! TEMPLATE search algorithm !!')
        print('***************** PARAMETERS *****************')
        print('NSTEP:', PARAMS.NSTEP)
        print('SEED: ', PARAMS.SEED)
        print('MINROTATION', PARAMS.MINROTATION)
        print('MAXROTATION', PARAMS.MAXROTATION)
        print('***************** RUNNING... *****************')
    starttime = time.time()

    template_search()

    if PARAMS.VERBOSE:
        print('***************** Finished *****************')
        print("--- Time Elapsed: %s seconds ---" % (time.time() - starttime))
    if MOL_INIT.SAVE_CONFS != None and MOL_INIT.SAVE_CONFS.GetNumConformers() != 0:
        if PARAMS.VERBOSE:
            print('---', MOL_INIT.SAVE_CONFS.GetNumConformers(), 'conformers saved ---')
    else:
        print('SIMULATION FAILED')
        sys.exit()

    return MOL_INIT.SAVE_CONFS

def main(Fin, Pin, name):
    dock = read_structure(Fin)
    errtags = PARAMS.read_parameter_file(Pin)
    if errtags:
        print('Could not parse parameters in file',errtags)
        sys.exit()

    save_out(name)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe Conformer-Generation Utility Package.',epilog='Probe conformer search by torsional/template search')
    parser.add_argument('-f', required=True) #structure file
    parser.add_argument('-p', required=True) #parameter file
    parser.add_argument('-o', required=False, default='SMART_ENS_') #output filename

    args = parser.parse_args()
    main(args.f, args.p, args.o)
