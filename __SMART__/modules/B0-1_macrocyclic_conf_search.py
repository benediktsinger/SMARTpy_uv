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

class TEMPLATE():
    def __init__(self, template=None, TEMPLATE_PATH=os.path.join(os.getcwd(),'templates')):
        try:
            # get template ENS
            self.MOLS = None
            for mol in Chem.SDMolSupplier(os.path.join(TEMPLATE_PATH, template), removeHs=False):
                if not self.MOLS:
                    self.MOLS = mol
                else:
                    self.MOLS.AddConformer(mol.GetConformer(), assignId=True)
        except Exception as e:
            print('Invalid template structure')
            print(e)
            sys.exit()

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

PARAMS.MINROTATION = 30
PARAMS.MAXROTATION = 330
PARAMS.SEED = None
PARAMS.NSTEP = 5

class CONFS:
    pass
CONFS.MOL = None
CONFS.PROBES = None
CONFS.ID = None

def ConfToMol(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
    return new_mol

def load_docked(inp):
# load docked structure
    docked = Chem.MolFromMolFile(inp, removeHs=False)
    return docked

def clashCheck(mol, probe):
    # check for probe-structure clashes
    for p in range(probe.GetNumAtoms()):
        for s in range(mol.GetNumAtoms()):
            i_ = probe.GetConformer().GetAtomPosition(p)
            i_r = Chem.GetPeriodicTable().GetRvdw(probe.GetAtomWithIdx(p).GetAtomicNum())
            j_ = mol.GetConformer().GetAtomPosition(s)
            j_r = Chem.GetPeriodicTable().GetRvdw(mol.GetAtomWithIdx(s).GetAtomicNum())
            d = (i_-j_).LengthSq()
            r_d = i_r + j_r# + TOL
            if d < r_d:

                return True
    return False

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

def save_out(name='out'):

    writer = Chem.SDWriter(name+'_pocket.sdf')
    print('writing conformers to file',name+'.sdf')
    for cid in [conf.GetId() for conf in CONFS.PROBES.GetConformers()]:
        writer.write(CONFS.PROBES, confId=cid)

    writer = Chem.SDWriter(name+'_mol.sdf')
    writer.write(CONFS.MOL)

def search(docked, template, step=1):
    RANDOM_ROT = random.randint(PARAMS.MINROTATION, PARAMS.MAXROTATION)
    frags = Chem.GetMolFrags(docked, asMols=True)
    probe = frags[0]
    #Chem.MolToMolFile(template.MOLS, 'step_probe'+str(step)+'.mol')
    mol = Chem.DeleteSubstructs(docked, frags[0])
    if not CONFS.MOL:
        CONFS.MOL = mol

    new_mol = Chem.Mol(mol)
    for p in range(template.MOLS.GetNumConformers()):
        step_probe = ConfToMol(template.MOLS, p)
        translation_v =  probe.GetConformer().GetAtomPosition(0) - step_probe.GetConformer().GetAtomPosition(0)
        for i in range(step_probe.GetNumAtoms()):
            pt = np.array(step_probe.GetConformer().GetAtomPosition(i))
            translate = pt+translation_v
            step_probe.GetConformer().SetAtomPosition(i, translate)
        if not clashCheck(mol, step_probe):
            new_mol = AllChem.CombineMols(step_probe, new_mol)
            if not CONFS.PROBES:
                CONFS.PROBES = step_probe
                print('-saving conformer-',1)
            else:
                CONFS.PROBES.AddConformer(step_probe.GetConformer(), assignId=True)
                print('-saving conformer-',CONFS.PROBES.GetNumConformers())
            #Chem.MolToMolFile(new_mol, 'step'+str(p)+str(step)+'.mol')
        rot = mol.GetConformer().GetAtomPosition(0) - probe.GetConformer().GetAtomPosition(0)
        rotn = np.linalg.norm(rot)
        rotvec = ( rot / rotn )
        rotation(template.MOLS, p, probe.GetConformer().GetAtomPosition(0), rotvec, RANDOM_ROT)

    if step <= PARAMS.NSTEP:
        step +=1
        search(docked, template, step)
    return

def start(docked, template):
    if PARAMS.SEED:
        random.seed(PARAMS.SEED)
    print('***************** PARAMETERS *****************')
    print('NSTEP:', PARAMS.NSTEP)
    print('SEED: ', PARAMS.SEED)
    print('MINROTATION', PARAMS.MINROTATION)
    print('MAXROTATION', PARAMS.MAXROTATION)
    print('***************** RUNNING... *****************')

    search(docked, template)

    print('***************** Finished *****************')
    print(CONFS.PROBES.GetNumConformers(), 'conformers saved')
    return CONFS.PROBES, CONFS.MOL

def main(Fin, Pin, Probe, name):
    docked = load_docked(Fin)
    template = TEMPLATE(Probe)
    errtags = PARAMS.read_parameter_file(Pin)

    if errtags:
        print('Could not parse parameters in file',errtags)
        sys.exit()

    start(docked, template)
    save_out(name)

if __name__ == '__main__':
    struc = 'test_cyclic_docked.mol'
    probe_name = 'S_SiF2_8_cyclic_MM.sdf'
    docked = load_docked(struc)
    template = TEMPLATE(probe_name)
    main(struc, '', probe_name)
#
