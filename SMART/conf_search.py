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
#from multiprocessing import Pool # TODO
from scipy.spatial.transform import Rotation as R
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit import ForceField
from rdkit import RDLogger
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')

import SMART as smart
from SMART import ReadProbe
from SMART.utils import conf_to_mol_, clash_check_, vectorize_, rotation_

TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)),'templates')
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
                elif i == 'REF_POS':
                    PARAMS.REF_POS = np.array(PAR[i])
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
PARAMS.REF_POS = None
PARAMS.VERBOSE = True
PARAMS.NPROCS = 4

class TEMPLATE():
    def __init__(self):
        super().__init__()
        pass
    def is_template(template):
        global TEMPLATE_PATH
        if os.path.isfile(os.path.join(TEMPLATE_PATH, template+'.sdf')):
            return True
        else:
            return False
    def GenerateTemplate(template, numConfs=200, pruneRmsThresh=0.1, overwrite=True):
        global TEMPLATE_PATH, PROBE_PATH
        try:
            starttime = time.time()
            # get probe MOL
            TEMPLATE.NAME = template
            TEMPLATE.MOLS = ReadProbe(os.path.join(PROBE_PATH, template+'.mol2')).MOL
            # generate template conformers
            AllChem.EmbedMultipleConfs(TEMPLATE.MOLS, numConfs=numConfs, pruneRmsThresh=pruneRmsThresh, useRandomCoords=True, useBasicKnowledge=True, useMacrocycleTorsions=True, useSmallRingTorsions=True)
            TEMPLATE.NCONFS = TEMPLATE.MOLS.GetNumConformers()
            # save template to sdf file
            if overwrite:
                if os.path.isfile(os.path.join(TEMPLATE_PATH, template+'.sdf')):
                    os.remove(os.path.join(TEMPLATE_PATH, template+'.sdf'))
            else:
                if os.path.isfile(os.path.join(TEMPLATE_PATH, template+'.sdf')):
                    print('already exits')
                    sys.exit()
            writer = Chem.SDWriter(os.path.join(TEMPLATE_PATH, template+'.sdf'))
            print('writing template to file: ', template+'.sdf')
            for cid in [conf.GetId() for conf in TEMPLATE.MOLS.GetConformers()]:
                writer.write(TEMPLATE.MOLS, confId=cid)
            print('***************** Done *****************')
            print("--- Template Construction Time: %s seconds ---" % (time.time() - starttime))
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
def template_search(step=1):
    if PARAMS.VERBOSE:
        print('step ',step,' / ',PARAMS.NSTEP)
    RANDOM_ROT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    ref_rotvec = (PARAMS.BINDING_POS - PARAMS.REF_POS) / np.linalg.norm(PARAMS.BINDING_POS - PARAMS.REF_POS)
    for conf in TEMPLATE.MOLS.GetConformers():
        rot_v = np.array(np.array(conf.GetAtomPosition(0)-conf.GetAtomPosition(TEMPLATE.MOLS.GetNumAtoms()-1)))/np.linalg.norm(np.array(np.array(conf.GetAtomPosition(0)-conf.GetAtomPosition(TEMPLATE.MOLS.GetNumAtoms()-1))))
        r = R.align_vectors([ref_rotvec], [rot_v])
        r[0].as_rotvec()
        # calculate vector-aligned atom position
        newv = r[0].apply([np.array(conf.GetAtomPosition(i)) for i in range(TEMPLATE.MOLS.GetNumAtoms())])
        for i in range(TEMPLATE.MOLS.GetNumAtoms()):
            # set vector-aligned atom position
            conf.SetAtomPosition(i, newv[i])
        #smart.ExportStructure(conf_to_mol_(TEMPLATE.MOLS, conf), str(conf.GetId())+'_rotate')
        # calculate translation vector
        translation_v =  PARAMS.BINDING_POS - np.array(conf.GetAtomPosition(0))
        for i in range(TEMPLATE.MOLS.GetNumAtoms()):
            pt = np.array(conf.GetAtomPosition(i))
            # calculate translated atom position
            translate = pt+translation_v
            # set to translated atom position
            conf.SetAtomPosition(i, translate)
        #smart.ExportStructure(conf_to_mol_(TEMPLATE.MOLS, conf), str(conf.GetId())+'_transl')
        # check for probe-structure clashes
        if not clash_check_(MOL_INIT.MOL, TEMPLATE.MOLS, conf.GetId(), PARAMS.CLASHTOL):
            if not MOL_INIT.SAVE_CONFS:
                MOL_INIT.SAVE_CONFS = conf_to_mol_(TEMPLATE.MOLS, conf)
                #MOL_INIT.SAVE_CONFS = conf_
                if PARAMS.VERBOSE:
                    print('-saving conformer-', 1)
            else:
                MOL_INIT.SAVE_CONFS.AddConformer(conf, assignId=True)
                if PARAMS.VERBOSE:
                    print('-saving conformer-', MOL_INIT.SAVE_CONFS.GetNumConformers())
        #break
    # rotate template for next iteration
    TEMPLATE.MOLS = rotation_(TEMPLATE.MOLS, math.radians(RANDOM_ROT), ref_rotvec)
    # end conditions
    if step < PARAMS.NSTEP:
        step += 1
        return template_search(step)
    return

############# ------- HANDLER FNXS
def TEMPLATE_SEARCH(mol=None):
    '''initiate conformational search by ensemble template'''
    if not mol:
        print('No docked MOL specified')
        sys.exit()

    if mol:
        MOL_INIT.MOL = mol.MOL
        MOL_INIT.SAVE_CONFS = None

        PARAMS.BINDING_POS = mol.Vector.BindingPos
        PARAMS.REF_POS = mol.Vector.TipPos
        #print(PARAMS.BINDING_POS, PARAMS.REF_POS)
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
            MOL_OUT = Chem.RWMol(MOL_INIT.SAVE_CONFS)
            MOL_OUT.RemoveAtom(TEMPLATE.MOLS.GetNumAtoms()-1)
            AllChem.SanitizeMol(MOL_OUT)
            MOL_INIT.SAVE_CONFS = MOL_OUT
    else:
        print('SIMULATION FAILED')



    return MOL_INIT.SAVE_CONFS
