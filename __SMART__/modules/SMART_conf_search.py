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
import sys, os, random, math, argparse, time
sys.setrecursionlimit(10000)
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit import ForceField
from rdkit import RDLogger
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')
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
                    PARAMS.SEED = int(PAR[i])
                elif i == 'NSTEP':
                    PARAMS.NSTEP = int(PAR[i])
                elif i == 'MAXROTATION':
                    PARAMS.MAXROTATION = float(PAR[i])
                elif i == 'MINROTATION':
                    PARAMS.MINROTATION = float(PAR[i])
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
PARAMS.NSTEP = 10 #500 for torsional; 20 for template
PARAMS.SEED = None
PARAMS.CLASHTOL = 0.0

class TEMPLATE():
    def __init__(self):
        pass
    def get_template(template):
        TEMPLATE_PATH=os.path.join(os.getcwd(),'templates')
        try:
            # get template ENS
            TEMPLATE.MOLS = None
            for mol in Chem.SDMolSupplier(os.path.join(TEMPLATE_PATH, template), removeHs=False):
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

class MOL_INIT:
    pass
MOL_INIT.FNAME = None
MOL_INIT.TYPE = None
MOL_INIT.MOL = None

class CONFS:
    pass
CONFS.MOL = None
CONFS.PROBES = None
CONFS.CPLX = None
CONFS.CPLX_CONFS = None

############# ------- MAIN FNXS
def save_out(name):
    writer = Chem.SDWriter(name+'_cavity.sdf')
    print('writing conformers to file',name+'_cavity.sdf')
    for cid in [conf.GetId() for conf in CONFS.PROBES.GetConformers()]:#range(1, CONFS.PROBES.GetNumConformers()):
        writer.write(CONFS.PROBES, confId=cid)
    writer = Chem.SDWriter(name+'_mol.sdf')
    writer.write(CONFS.MOL)
    writer = Chem.SDWriter(name+'_cplx.sdf')
    writer.write(CONFS.CPLX)

def ConfToMol(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
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

def clash_check(mol, probe):
    #TOL = 0.01

    for p in range(probe.GetNumAtoms()):
        for s in range(mol.GetNumAtoms()):
                i_ = probe.GetConformer().GetAtomPosition(p)
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

def rotation(probe, cid, tip, rotvec, theta):
    for i in range(probe.GetNumAtoms()):
        a = np.array(probe.GetConformer(cid).GetAtomPosition(i))
        dotproduct = rotvec[0]*(float(a[0]) - float(tip[0])) + rotvec[1]*(float(a[1]) - float(tip[1])) + rotvec[2]*(float(a[2]) - float(tip[2]))
        centre = [float(tip[0]) + dotproduct*rotvec[0], float(tip[1]) + dotproduct*rotvec[1], float(tip[2]) + dotproduct*rotvec[2]]
        v = [float(a[0]) - centre[0], float(a[1]) - centre[1], float(a[2]) - centre[2]]
        d = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        px = v[0]*math.cos(theta) + v[1]*math.sin(theta)*rotvec[2] - v[2]*math.sin(theta)*rotvec[1]
        py = v[1]*math.cos(theta) + v[2]*math.sin(theta)*rotvec[0] - v[0]*math.sin(theta)*rotvec[2]
        pz = v[2]*math.cos(theta) + v[0]*math.sin(theta)*rotvec[1] - v[1]*math.sin(theta)*rotvec[0]
        newv = [px + centre[0], py + centre[1], pz + centre[2]]

        probe.GetConformer(cid).SetAtomPosition(i, newv)

def get_ns(atom, ref, atom_ns=None, atom_save=None, atom_n=None):
    # for a in atom_ns (bound to > 1 atom)
    def it(a):
        ns, save = [], []
        an = a.GetNeighbors()
        if len(an) > 0:
            for j in an:
                if j.GetIdx() == ref:
                    continue
                elif atom_ns and j.GetIdx() in [k.GetIdx() for k in atom_ns]:
                    continue
                elif atom_save and j.GetIdx() in [k.GetIdx() for k in atom_save]:
                    continue

                if j.GetExplicitValence() > 1:
                    ns.append(j)
                    save.append(j)
                elif j.GetExplicitValence() == 1:
                    save.append(j)
        return ns, save

    if atom_ns:
        ns_, save_ = [], []
        for a in atom_ns:
            ns, save = it(a)
            ns_.extend(ns)
            save_.extend(save)
        if len(save_) > 0:
            atom_ns.extend(ns_)
            atom_save.extend(save_)
            return get_ns(atom, ref, atom_ns, atom_save, atom_n)
        else:
            return atom_n, atom_save
    else:
        atom_ns, atom_save = it(atom)
        atom_n = atom_ns[0]
        atom_ns.append(atom)
        return get_ns(atom, ref, atom_ns, atom_save, atom_n)

def torsional_step(step=1):
    frags = Chem.GetMolFrags(CONFS.CPLX, asMols=True)
    probe = frags[0]
    structure = Chem.DeleteSubstructs(CONFS.CPLX, probe)
    if not CONFS.MOL:
        CONFS.MOL = structure

    bonds = []
    for b in probe.GetBonds():
        type = b.GetBondType()
        # only rotate single bonds
        if str(type) != 'SINGLE':
            continue
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()
        # only rotate dihedrals
        if a1.GetExplicitValence() <= 1 or a2.GetExplicitValence() <= 1:
            continue
        a1n = sorted([i.GetExplicitValence() for i in a1.GetNeighbors()])
        a2n = sorted([j.GetExplicitValence() for j in a2.GetNeighbors()])
        a1n.pop()
        a2n.pop()
        # do not rotate CH3/CF3/etc.
        if a1n[-1] <= 1 or a2n[-1] <= 1:
            continue
        bonds.append(b)

    RANDOM_TORSION = random.choice(bonds)
    RANDOM_DISPLACMENT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    RANDOM_PROBABILITY = random.random()

    # find atoms for displacment
    a1 = RANDOM_TORSION.GetBeginAtom()
    a2 = RANDOM_TORSION.GetEndAtom()

    # isolate structure fragments
    # OPTIMIZE THIS STEP
    atom1_n, atom1_save = get_ns(a1, a2.GetIdx())
    atom2_n, atom2_save = get_ns(a2, a1.GetIdx())

    # rotation vector a1->a2
    a1_pos = probe.GetConformer().GetAtomPosition(a1.GetIdx())
    a2_pos = probe.GetConformer().GetAtomPosition(a2.GetIdx())
    d =  a2_pos - a1_pos
    rotvec = d/np.linalg.norm(d)
    if len(atom1_save) < len(atom2_save):
        DISP = atom1_save
        tip = a1_pos
    else:
        DISP = atom2_save
        tip = a2_pos
    vecs = vectorize(probe, tip, DISP, rotvec, RANDOM_DISPLACMENT)

    # update frag
    probe_upd = probe
    for i in probe_upd.GetAtoms():
        if i.GetIdx() in [j.GetIdx() for j in DISP]:
            pos = [j.GetIdx() for j in DISP].index(i.GetIdx())
            pt = Point3D(vecs[pos][0], vecs[pos][1], vecs[pos][2])
            probe_upd.GetConformer().SetAtomPosition(i.GetIdx(), vecs[pos])

    candidate = Chem.CombineMols(probe_upd, structure)
    AllChem.SanitizeMol(candidate)

    # check for clashes
    if clash_check(candidate):
        step += 1
        return torsional_step(step)

    if not CONFS.PROBES:
        CONFS.PROBES = probe_upd
        print('-saving conformer-',1)
    else:
        CONFS.PROBES.AddConformer(probe_upd.GetConformer(), assignId=True)
        print('-saving conformer-',CONFS.PROBES.GetNumConformers())

    if step < PARAMS.NSTEP:
        CONFS.CPLX = candidate
        if not CONFS.CPLX_CONFS:
            CONFS.CPLX_CONFS = candidate
        else:
            CONFS.CPLX_CONFS.AddConformer(candidate.GetConformer(), assignId=True)
        step += 1
        return torsional_step(step)
    return

def custom_template(step=1):
    print('generating template')
    frags = Chem.GetMolFrags(CONFS.CPLX, asMols=True)
    probe = frags[0]

    AllChem.EmbedMultipleConfs(probe, numConfs=200, pruneRmsThresh=0.1, useRandomCoords=True, useMacrocycleTorsions=True, useSmallRingTorsions=True)
    #rdMolAlign.alignMolConformers(probe, atomIds=[0])
    print('searching...')
    TEMPLATE.MOLS = probe
    template_search()

def template_search(step=1):
    RANDOM_ROT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    frags = Chem.GetMolFrags(CONFS.CPLX, asMols=True)
    probe = frags[0]
    structure = Chem.DeleteSubstructs(CONFS.CPLX, probe)
    if not CONFS.MOL:
        CONFS.MOL = structure

    new_mol = Chem.Mol(structure)
    for p in range(TEMPLATE.MOLS.GetNumConformers()):
        probe_upd = ConfToMol(TEMPLATE.MOLS, p)
        translation_v =  probe.GetConformer().GetAtomPosition(0) - probe_upd.GetConformer().GetAtomPosition(0)
        for i in range(probe_upd.GetNumAtoms()):
            pt = np.array(probe_upd.GetConformer().GetAtomPosition(i))
            translate = pt+translation_v
            probe_upd.GetConformer().SetAtomPosition(i, translate)
        # check for probe-structure clashes
        if not clash_check(structure, probe_upd):
            candidate = Chem.CombineMols(probe_upd, structure)
            if not CONFS.PROBES:
                CONFS.PROBES = probe_upd
                print('-saving conformer-',1)
            else:
                CONFS.PROBES.AddConformer(probe_upd.GetConformer(), assignId=True)
                print('-saving conformer-',CONFS.PROBES.GetNumConformers())
            CONFS.CPLX = candidate
            if not CONFS.CPLX_CONFS:
                CONFS.CPLX_CONFS = candidate
            else:
                CONFS.CPLX_CONFS.AddConformer(candidate.GetConformer(), assignId=True)
        # rotate probe conformer
        rot = structure.GetConformer().GetAtomPosition(0) - probe.GetConformer().GetAtomPosition(0)
        rotn = np.linalg.norm(rot)
        rotvec = ( rot / rotn )
        rotation(TEMPLATE.MOLS, p, probe.GetConformer().GetAtomPosition(0), rotvec, RANDOM_ROT)

    # end conditions
    if step <= PARAMS.NSTEP:
        step +=1
        return template_search(step)
    return

############# ------- HANDLER FNXS
def TEMPLATE_SEARCH(mol=None):
    '''initiate conformational search by ensemble template'''
    if not mol:
        print('No docked MOL specified')
        sys.exit()

    CONFS.PROBES = None
    MOL_INIT.MOL = mol
    CONFS.CPLX = MOL_INIT.MOL
    CONFS.CPLX_CONFS = None

    # print search parameters
    if PARAMS.SEED:
        random.seed(PARAMS.SEED)
    print('!! TEMPLATE search algorithm !!')
    print('***************** PARAMETERS *****************')
    print('NSTEP:', PARAMS.NSTEP)
    print('SEED: ', PARAMS.SEED)
    print('MINROTATION', PARAMS.MINROTATION)
    print('MAXROTATION', PARAMS.MAXROTATION)
    print('***************** RUNNING... *****************')
    starttime = time.time()
    template_search()

    print('***************** Finished *****************')
    print("--- Time Elapsed: %s seconds ---" % (time.time() - starttime))
    if CONFS.PROBES != None:
        print('---',CONFS.PROBES.GetNumConformers(), 'conformers saved ---')
    else:
        print('SIMULATION FAILED')

    return CONFS.PROBES, CONFS.MOL, CONFS.CPLX_CONFS

def CUSTOM_TEMPLATE_SEARCH(mol=None):
    '''generate custom probe template and follow TEMPLATE_SEARCH protocol'''
    # initialize variables
    if not mol:
        print('No docked MOL structure specified')
        sys.exit()

    CONFS.PROBES = None
    MOL_INIT.MOL = mol
    CONFS.CPLX = MOL_INIT.MOL
    CONFS.CPLX_CONFS = None

    # print search parameters
    if PARAMS.SEED:
        random.seed(PARAMS.SEED)
    print('!! CUSTOM TEMPLATE search algorithm !!')
    print('***************** PARAMETERS *****************')
    print('NSTEP:', PARAMS.NSTEP)
    print('SEED: ', PARAMS.SEED)
    print('MINROTATION', PARAMS.MINROTATION)
    print('MAXROTATION', PARAMS.MAXROTATION)
    print('***************** RUNNING... *****************')
    # start custom template search
    starttime = time.time()
    custom_template()

    print('***************** Finished *****************')
    print("--- Time Elapsed: %s seconds ---" % (time.time() - starttime))
    if CONFS.PROBES != None:
        print('---',CONFS.PROBES.GetNumConformers(), 'conformers saved ---')
    else:
        print('SIMULATION FAILED')

    return CONFS.PROBES, CONFS.MOL, CONFS.CPLX_CONFS

def TORSIONAL_STEP_SEARCH(mol=None):
    '''initiate conformational search by random torsion selection'''
    # initialize variables
    if not mol:
        print('No docked MOL structure specified')
        sys.exit()

    CONFS.PROBES = None
    MOL_INIT.MOL = mol
    CONFS.CPLX = MOL_INIT.MOL
    CONFS.CPLX_CONFS = None

    # print search parameters
    if PARAMS.SEED:
        random.seed(PARAMS.SEED)
    print('!! TORSIONAL STEP search algorithm !!')
    print('***************** PARAMETERS *****************')
    print('NSTEP:', PARAMS.NSTEP)
    print('SEED: ', PARAMS.SEED)
    print('MINROTATION', PARAMS.MINROTATION)
    print('MAXROTATION', PARAMS.MAXROTATION)
    print('***************** RUNNING... *****************')
    # start stepwise search
    starttime = time.time()
    torsional_step()

    print('***************** Finished *****************')
    print("--- Time Elapsed: %s seconds ---" % (time.time() - starttime))
    if CONFS.PROBES != None:
        print('---',CONFS.PROBES.GetNumConformers(), 'conformers saved ---')
    else:
        print('SIMULATION FAILED')

    return CONFS.PROBES, CONFS.MOL, CONFS.CPLX_CONFS

def main(Fin, Pin, PROTOC, name):
    dock = read_structure(Fin)
    errtags = PARAMS.read_parameter_file(Pin)
    if errtags:
        print('Could not parse parameters in file',errtags)
        sys.exit()

    if PROTOC == 1 or PROTOC.lower() == 'torsional':
        TORSIONAL_STEP_SEARCH(dock)
    elif PROTOC == 2 or PROTOC.lower() == 'template':
        TEMPLATE_SEARCH(dock)
    elif PROTOC == 3 or PROTOC.lower() == 'custom':
        CUSTOM_TEMPLATE_SEARCH(dock)
    else:
        print('Please select search protocol:  [ 1 (torsional), 2 (template), 3 (custom) ]')
        sys.exit()

    save_out(name)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe Conformer-Generation Utility Package.',epilog='Probe conformer search by torsional/template search')
    parser.add_argument('-f', required=True) #structure file
    parser.add_argument('-p', required=True) #parameter file
    parser.add_argument('-method', required=True, choices=[1,2,3,'torsional','template','custom']) #search protocol
    parser.add_argument('-o', required=False, default='SMART_conformers_') #output filename

    args = parser.parse_args()

    starttime = time.time()
    main(args.f, args.p, args.method, args.o)
    print("--- %s seconds ---" % (time.time() - starttime))
