#!/usr/bin/python
# SMART probe conformational search utility
#   version: b0.1 (released 10-10-2023)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#
# Dev Notes:
#  Features of FullMonte (https://github.com/patonlab/FullMonte) inspired development of this code.
#  - FullMonte is an archived public GitHub repository.
#  - Robert Paton (https://github.com/patonlab/)
#
############# ------- DEPENDENCIES
import sys, os, random, math, argparse, time
sys.setrecursionlimit(10000)
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
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
                print(i)
                if i == 'DEFAULT':
                    continue
                    #PARAMS.DEFAULT = PAR[i]
                elif i == 'FIXAT':
                    PARAMS.FIXAT = []
                    for j in PAR[i]:
                        if int(j) < 0:
                            rng = list(range(abs(int(j)), int(PAR[i][PAR[i].index(j)+1])+1))
                            PARAMS.FIXAT.extend(rng)
                        else:
                            if int(j) not in PARAMS.FIXAT:
                                PARAMS.FIXAT.append(int(j))
                elif i == 'SEED':
                    PARAMS.SEED = int(PAR[i])
                elif i == 'FREEAT':
                    continue
                elif i == 'NSTEP':
                    PARAMS.NSTEP = int(PAR[i])
                elif i == 'MAXROTATION':
                    PARAMS.MAXROTATION = float(PAR[i])
                elif i == 'MINROTATION':
                    PARAMS.MINROTATION = float(PAR[i])
                elif i == 'EWIN':
                    PARAMS.EWIN = float(PAR[i])
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

PARAMS.DEFAULT = 0
PARAMS.FREEAT = []
PARAMS.FIXAT = []
PARAMS.MAXROTATION = 330
PARAMS.MINROTATION = 30
PARAMS.NSTEP = 500
#PARAMS.EWIN = 20
PARAMS.SEED = None

class MOL_INIT:
    pass
MOL_INIT.FNAME = None
MOL_INIT.TYPE = None
MOL_INIT.MOL = None
#MOL_INIT.ENERGY = 0
#MOL_INIT.GMIN = 0

class CONFS:
    pass
CONFS.MOL = None
CONFS.CONFS = None
CONFS.ID = None
#CONFS.ENERGY = 0
#CONFS.PROBABILITY = 0

############# ------- MAIN FNXS
def CLASH(conf):
    TOL = 0.01
    #clashes = 0
    for i in PARAMS.FIXAT:
        for j in range(1,conf.GetNumAtoms()+1):
            if i != j and j not in PARAMS.FIXAT:
                i_ = conf.GetConformer().GetAtomPosition(i-1)
                i_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(i-1).GetAtomicNum())
                j_ = conf.GetConformer().GetAtomPosition(j-1)
                j_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(j-1).GetAtomicNum())
                d = (i_-j_).LengthSq()
                r_d = i_r + j_r #+ TOL
                if (d - r_d) <= 0:
                #if d <= 1.25:
                    #print(i,j,d)
                    return True
    return False

#def G_MIN(conf, conf_energy):
#    if conf_energy < MOL_INIT.GMIN:
#        print('-new GMIN found-','NEW ENERGY:',conf_energy)
#        return True
#    else:
#        return False

#def E_DIFF(conf, conf_energy):
#    diff = CONFS.ENERGY - conf_energy
#    if diff >= 0:
#        print('-save by energy decrease-','NEW ENERGY:',conf_energy)
#        return True
#    elif conf_energy - MOL_INIT.GMIN >= PARAMS.EWIN: # check high energy
#        return False
#    return False

#def PROBABILITY(conf_energy, rand):
#    diff = CONFS.ENERGY - conf_energy
    #print(diff)
#    conf_probability = 1 / (1 + math.exp(diff/(298.15*8.314462618*(10**-3))))
    #print('CONF PROB:', conf_probability, '-compare to', rand)
    ##larger value suceeds
#    if conf_probability >= rand:
#        print('-save by probability-', 'NEW ENERGY:', conf_energy)
#        return True
#    else:
#        return False

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

def monte_carlo(step=1):

    frags = Chem.GetMolFrags(CONFS.MOL, asMols=True)
    bonds = []
    for b in frags[0].GetBonds():
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
        # dont rotate CH3/CF3/etc.
        if a1n[-1] <= 1 or a2n[-1] <= 1:
            continue
        bonds.append(b)

    RANDOM_TORSION = random.choice(bonds)
    RANDOM_DISPLACMENT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    RANDOM_PROBABILITY = random.random()

    # find atoms for displacment
    a1 = RANDOM_TORSION.GetBeginAtom()
    a2 = RANDOM_TORSION.GetEndAtom()

    atom1_n, atom1_save = get_ns(a1, a2.GetIdx())
    atom2_n, atom2_save = get_ns(a2, a1.GetIdx())

    # rotation vector a1->a2
    a1_pos = frags[0].GetConformer().GetAtomPosition(a1.GetIdx())
    a2_pos = frags[0].GetConformer().GetAtomPosition(a2.GetIdx())
    d =  a2_pos - a1_pos
    rotvec = d/np.linalg.norm(d)
    if len(atom1_save) < len(atom2_save):
        DISP = atom1_save
        tip = a1_pos
    else:
        DISP = atom2_save
        tip = a2_pos
    vecs = vectorize(frags[0], tip, DISP, rotvec, RANDOM_DISPLACMENT)

    # update frag
    frag_ = frags[0]
    for i in frag_.GetAtoms():
        if i.GetIdx() in [j.GetIdx() for j in DISP]:
            pos = [j.GetIdx() for j in DISP].index(i.GetIdx())
            pt = Point3D(vecs[pos][0], vecs[pos][1], vecs[pos][2])
            frag_.GetConformer().SetAtomPosition(i.GetIdx(), vecs[pos])

    conf = Chem.DeleteSubstructs(CONFS.MOL, frags[0])
    conf = Chem.CombineMols(frag_, conf)
    AllChem.SanitizeMol(conf)
    #conf_energy = energy(conf)
    #print('NEW ENERGY:',conf_energy)

    # PASS/FAIL CONDITIONS
    #check for clashes
    if CLASH(conf):
        CONFS.MOL = conf
        step += 1
        return monte_carlo(step)
    # check for energy decrease
    #elif not E_DIFF(conf, conf_energy):
        #check probability of transition
    #    if not PROBABILITY(conf_energy, RANDOM_PROBABILITY):
    #        CONFS.MOL = conf
    #        step += 1
    #        return monte_carlo(step)

    #check for global minimum energy
    #if G_MIN(conf, conf_energy):
    #    MOL_INIT.GMIN = conf_energy

    if not CONFS.CONFS:
        CONFS.CONFS = conf
        print('-saving conformer-',1)
    else:
        CONFS.CONFS.AddConformer(conf.GetConformer(), assignId=True)
        print('-saving conformer-',CONFS.CONFS.GetNumConformers())

    #CONFS.ENERGY = conf_energy
    if step < PARAMS.NSTEP:
        CONFS.MOL = conf
        step += 1
        return monte_carlo(step)
    return

#def energy(mol):
#    if AllChem.UFFHasAllMoleculeParams(mol):
#        ff = AllChem.UFFGetMoleculeForceField(mol)
        #for i in PARAMS.FIXAT:
        #    ff.UFFAddPositionConstraint(i-1, 0, 1.e4)
        #ff.energy(maxIts=50)
#    else:
#        print("\nFATAL ERROR"%file)
#        sys.exit()

#    return ff.CalcEnergy()

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
    return MOL_INIT.MOL

def save_out(name):

    writer = Chem.SDWriter(name+'.sdf')
    print('writing conformers to file',name+'.sdf')
    for cid in range(1, CONFS.CONFS.GetNumConformers()):

        writer.write(CONFS.CONFS, confId=cid)

def start(mol, name='out'):
    # initialize variables
    CONFS.CONFS = None
    CONFS.ENERGY = 0
    ##########################
    MOL_INIT.MOL = mol
    #MOL_INIT.ENERGY = energy(MOL_INIT.MOL)
    #MOL_INIT.GMIN = MOL_INIT.ENERGY
    #print('INIT NRG: '+str(MOL_INIT.ENERGY))

    CONFS.MOL = MOL_INIT.MOL
    #CONFS.ENERGY = MOL_INIT.ENERGY
    if PARAMS.SEED:
        random.seed(PARAMS.SEED)
    print('***************** PARAMETERS *****************')
    #print('EWIN:', PARAMS.EWIN)
    print('NSTEP:', PARAMS.NSTEP)
    print('SEED: ', PARAMS.SEED)
    print('MINROTATION', PARAMS.MINROTATION)
    print('MAXROTATION', PARAMS.MAXROTATION)
    print('FIXAT: ', PARAMS.FIXAT)
    print('***************** RUNNING... *****************')

    # search
    monte_carlo()

    print('***************** Finished *****************')
    print(CONFS.CONFS.GetNumConformers(), 'conformers saved')
    return CONFS.CONFS

def main(Fin, Pin, name):
    dock = read_structure(Fin)
    errtags = PARAMS.read_parameter_file(Pin)

    if errtags:
        print('Could not parse parameters in file',errtags)
        sys.exit()

    start(dock, name)
    save_out(name)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Pocket generation Utility Package.',epilog='Probe conformer search by modified Monte Carlo algorithm')
    parser.add_argument('-f', required=True) #structure file
    parser.add_argument('-p', required=True) #parameter file
    parser.add_argument('-o', required=False, default='SMART_conformers_') #output filename

    args = parser.parse_args()

    starttime = time.time()
    main(args.f, args.p, args.o)
    print("--- %s seconds ---" % (time.time() - starttime))
