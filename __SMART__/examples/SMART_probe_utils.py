#!/usr/bin/python
# SMART probe addition utility
#   version: b1.1 (released 01-27-2024)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#
############# ------- DEPENDENCIES
import sys, os, argparse, math, mathutils, time
import numpy as np
from scipy.spatial.transform import Rotation as R

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import ForceField
from rdkit import RDLogger
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')
pt = Chem.GetPeriodicTable()

############# ------- UTILITY FUNCTIONS
def octahedral(mol, metalID, bindingAtoms):
    for a in binding_atoms:
        n = 0
        min, max, avg = 360, 0, 0
        for b in binding_atoms:
            if a == b:
                continue
            n += 1
            angle = rdMolTransforms.GetAngleDeg(mol.GetConformer(), a-1, metal_id-1, b-1)
            avg += angle
            if angle > max:
                max = angle
            if angle < min:
                min = angle
        #print('Min:',min,'Max:',max,'Avg:',avg/n)
        if max < (180-TOL):
            return a-1

def tetrahedral(mol, metalID, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[0]-1)) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[1]-1)) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[2]-1)) ) / 3
    return ref_pos

def trigonal(mol, metalID, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[0]-1)) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[1]-1)) ) / 2
    return ref_pos

############# ------- UTILITY CLASSES
class ReadFile:
    # initialize single structure coordinates
    def __init__(self, ext, fname, path=os.getcwd()):
    #file type, name, (path)
        try:
            self.NAME = fname
            self.TYPE = ext

            file = os.path.join(path, fname+'.'+ext)
            if self.TYPE == 'mol2':
                self.MOL = Chem.MolFromMol2File(file, removeHs=False)
            elif self.TYPE == 'mol':
                self.MOL = Chem.MolFromMolFile(file, removeHs=False)
            elif self.TYPE == 'xyz':
                self.MOL = Chem.MolFromXYZFile(file, removeHs=False)
            elif self.TYPE == 'pdb':
                self.MOL = Chem.MolFromPDBFile(file, removeHs=False)
            elif self.TYPE == 'sdf':
                self.MOL = Chem.SDSupplier(file, removeHs=False)[0]
            self.ATOMS = self.MOL.GetAtoms()
            self.NATOMS = self.MOL.GetNumAtoms()
        except Exception as e:
            print('no structure file found -', file)
            print(e)
            sys.exit()

    def Reference_Vector(self, tip, tail=None, plane=None):
        if tail:
            self.TAIL_N = tail
            self.TIP_N = tip
            self.Vector = BindingVector(self.MOL, self.TIP_N, self.TAIL_N)
        else:
            self.TIP_N = tip
            self.PLN1_N = plane[0]
            self.PLN2_N = plane[1]
            self.Vector = BindingVector(self.MOL, self.TIP_N, [self.PLN1_N, self.PLN2_N])

    def Assess_Geometry(self, geom=None, metalID=None, bindingAtoms=None, TOL=20):
        if geom == 'octahedral':
            if not len(bindingAtoms) == 5:
                print('Cannot compute octahedral geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, octahedral(self.MOL, metalID, bindingAtoms))
        elif geom == 'tetrahedral':
            if not len(bindingAtoms) == 3:
                print('Cannot compute tetrahedral geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, tetrahedral(self.MOL, metalID, bindingAtoms))
        elif geom ==  'trigonal':
            if not len(bindingAtoms) == 2:
                print('Cannot compute trigonal geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, trigonal(self.MOL, metalID, bindingAtoms))
        else:
            print('please a specify geometry')
            sys.exit()

class Structure:
    # initialize single structure coordinates
    def __init__(self, name, mol, path=os.getcwd()):
    #file type, name, (path)
        try:
            self.NAME = name
            self.MOL = mol

            self.ATOMS = self.MOL.GetAtoms()
            self.NATOMS = self.MOL.GetNumAtoms()
        except Exception as e:
            print('no structure file found -', file)
            print(e)
            sys.exit()

    def Reference_Vector(self, tip, tail=None, plane=None):
        if tail:
            self.TIP_N = tip
            self.TAIL_N = tail
            self.Vector = BindingVector(self.MOL, self.TIP_N, self.TAIL_N)
        elif plane:
            self.TIP_N = tip
            self.PLN1_N = plane[0]
            self.PLN2_N = plane[1]
            self.Vector = BindingVector(self.MOL, self.TIP_N, [self.PLN1_N, self.PLN2_N])
        else:
            print('specify reference atom(s)')
            sys.exit()

    def Assess_Geometry(self, geom=None, metalID=None, bindingAtoms=None, TOL=20):
        if geom == 'octahedral':
            if not len(bindingAtoms) == 5:
                print('Cannot compute octahedral geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, octahedral(self.MOL, metalID, bindingAtoms))
        elif geom == 'tetrahedral':
            if not len(bindingAtoms) == 3:
                print('Cannot compute tetrahedral geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, tetrahedral(self.MOL, metalID, bindingAtoms))
        elif geom ==  'trigonal':
            if not len(bindingAtoms) == 2:
                print('Cannot compute trigonal geometry from these atoms:')
                print(bindingAtoms)
                sys.exit()
            self.Vector = GeometryVector(self.MOL, self.TIP_N, trigonal(self.MOL, metalID, bindingAtoms))
        else:
            print('please a specify geometry')
            sys.exit()

class Probe:
    # initialize built-in probe
    def __init__(self, probe):
    # probe name
        try:
            PROBE_PATH = os.path.join(os.path.dirname(__file__), 'Probes')
            self.NAME = probe
            file = os.path.join(PROBE_PATH, probe+'.mol2')
            self.MOL = Chem.MolFromMol2File(file, removeHs=False)
            self.ATOMS = self.MOL.GetAtoms()
            self.NATOMS = self.MOL.GetNumAtoms()
            self.TIP_N = 0
            self.TAIL_N = self.NATOMS-1
        except Exception as e:
            print('no probe found -', file)
            print(e)
            sys.exit()
        self.Vector = BindingVector(self.MOL, self.TIP_N, self.TAIL_N)

class GeometryVector:
    def __init__(self, MOL, tip, tail_pos):
        try:
            self.TAIL_POS = tail_pos
            self.TIP_POS = MOL.GetConformer().GetAtomPosition(int(tip))
            v = (self.TIP_POS - self.TAIL_POS)
            n = np.linalg.norm(v)
            self.U = (v / n)
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()

class BindingVector:
    # calculate binding vectors
    def __init__(self, MOL, tip, tail):
    # ReadFiles/Probe.mol, tip atom number, tail atom number
        try:
            if type(tail) == list:
                self.TIP_POS = MOL.GetConformer().GetAtomPosition(int(tip))
                v1 = MOL.GetConformer().GetAtomPosition(int(tail[0])) - MOL.GetConformer().GetAtomPosition(int(tip))
                n1 = np.linalg.norm(v1)
                u1 = (v1 / n1)
                v2 = MOL.GetConformer().GetAtomPosition(int(tail[1])) - MOL.GetConformer().GetAtomPosition(int(tip))
                n2 = np.linalg.norm(v2)
                u2 = (v2 / n2)
                self.U = np.cross(u1,u2)
            else:
                self.TAIL_POS = MOL.GetConformer().GetAtomPosition(int(tail))
                self.TIP_POS = MOL.GetConformer().GetAtomPosition(int(tip))
                v = (self.TIP_POS - self.TAIL_POS)
                n = np.linalg.norm(v)
                self.U = (v / n)
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()

class ExportStructure:
    def __init__(self, mol, outname='out'):
        Chem.MolToMolFile(mol, outname+'.mol')

############# ------- MAIN FNXS
def clashCheck(conf, struc, probe):
    # check for probe-structure clashes
    for p in range(probe.NATOMS):
        if p == probe.TIP_N or p == probe.TAIL_N:
            continue
        for s in range(probe.NATOMS,probe.NATOMS+struc.NATOMS-1):
            if s == struc.TIP_N:
                #print(s,'bind')
                continue
            i_ = conf.GetConformer().GetAtomPosition(p)
            i_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(p).GetAtomicNum())
            j_ = conf.GetConformer().GetAtomPosition(s)
            j_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(s).GetAtomicNum())
            d = (i_-j_).LengthSq()
            r_d = i_r + j_r# + TOL
            if d < r_d:
                return True
    return False

def optimizeFit(struc, probe, incr, dist, theta=5):
    # rotate probe to fit structure
    conf = Chem.CombineMols(probe.MOL, struc.MOL)
    conf_ = Chem.RWMol(conf)
    conf_.RemoveAtom(probe.TAIL_N)
    AllChem.SanitizeMol(conf_)
    #if AllChem.UFFHasAllMoleculeParams(conf_):
    #    ff = AllChem.UFFGetMoleculeForceField(conf_)
    #    ff.UFFAddPositionConstraint(0, 0, 1.e4)
    #    for i in range(probe.NATOMS-1, probe.NATOMS+struc.NATOMS):
    #        ff.UFFAddPositionConstraint(i-1, 0, 1.e4)
    #    ff.Minimize(maxIts=100)
    #else:
    #    print("\nFATAL ERROR"%file)
    #    sys.exit()

    #return conf_
    #Chem.MolToMolFile(conf_, str(theta)+'conf_.mol')

    if not clashCheck(conf_, struc, probe):
        return conf_
    else:
        if theta < 360:
            theta += incr
            #probe_r = rotation(probe, probe.Vector.U, theta)
            probe_r = rotation2(probe, np.array(struc.Vector.TIP_POS), probe.Vector.U, math.radians(theta))
            return optimizeFit(struc, probe_r, incr, dist, theta)
        else:
            print('failed:',struc.NAME)
            return conf_
    return conf_

def rotateAlign(struc, probe):
    r = R.align_vectors([struc.Vector.U], [probe.Vector.U])
    r[0].as_matrix()
    newv = r[0].apply([np.array(probe.MOL.GetConformer().GetAtomPosition(i)) for i in range(probe.NATOMS)])

    for i in range(probe.NATOMS):
        probe.MOL.GetConformer().SetAtomPosition(i, newv[i])

    probe.Vector.TIP_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TIP_N)
    probe.Vector.TAIL_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TAIL_N)
    v = probe.Vector.TIP_POS  - probe.Vector.TAIL_POS
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def rotation2(probe, tip, rotvec, theta):
    for i in range(probe.NATOMS):
        a = np.array(probe.MOL.GetConformer().GetAtomPosition(i))
        dotproduct = rotvec[0]*(float(a[0]) - float(tip[0])) + rotvec[1]*(float(a[1]) - float(tip[1])) + rotvec[2]*(float(a[2]) - float(tip[2]))
        centre = [float(tip[0]) + dotproduct*rotvec[0], float(tip[1]) + dotproduct*rotvec[1], float(tip[2]) + dotproduct*rotvec[2]]
        v = [float(a[0]) - centre[0], float(a[1]) - centre[1], float(a[2]) - centre[2]]
        d = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        px = v[0]*math.cos(theta) + v[1]*math.sin(theta)*rotvec[2] - v[2]*math.sin(theta)*rotvec[1]
        py = v[1]*math.cos(theta) + v[2]*math.sin(theta)*rotvec[0] - v[0]*math.sin(theta)*rotvec[2]
        pz = v[2]*math.cos(theta) + v[0]*math.sin(theta)*rotvec[1] - v[1]*math.sin(theta)*rotvec[0]
        newv = [px + centre[0], py + centre[1], pz + centre[2]]

        probe.MOL.GetConformer().SetAtomPosition(i, newv)

    probe.Vector.TIP_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TIP_N)
    probe.Vector.TAIL_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TAIL_N)
    v = probe.Vector.TIP_POS  - probe.Vector.TAIL_POS
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def rotation(probe, rotvec, theta):
    # rotate 3D about vector
    quat = mathutils.Quaternion(rotvec, math.radians(theta))
    r = R.from_quat(quat)
    r.as_quat()
    newv = r.apply([np.array(probe.MOL.GetConformer().GetAtomPosition(i)) for i in range(probe.NATOMS)])
    for i in range(probe.NATOMS):
        probe.MOL.GetConformer().SetAtomPosition(i, newv[i])

    probe.Vector.TIP_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TIP_N)
    probe.Vector.TAIL_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TAIL_N)
    v = probe.Vector.TIP_POS  - probe.Vector.TAIL_POS
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def translate(struc, probe, dist):
    # translate to binding Point3D
    binding_point = np.array(struc.Vector.TIP_POS) + (float(dist)*np.array(struc.Vector.U))
    translation_v = binding_point - np.array(probe.Vector.TIP_POS)

    for i in range(probe.NATOMS):
        pt = np.array(probe.MOL.GetConformer().GetAtomPosition(i))
        translate = pt+translation_v
        probe.MOL.GetConformer().SetAtomPosition(i, translate)

    probe.Vector.TIP_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TIP_N)
    probe.Vector.TAIL_POS = probe.MOL.GetConformer().GetAtomPosition(probe.TAIL_N)
    v = probe.Vector.TIP_POS  - probe.Vector.TAIL_POS
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def addProbe(struc, probe, dist=2.0, rot=5.0):
    # main utility function
    probe_r = rotateAlign(struc, probe)
    #Chem.MolToMolFile(probe_r.MOL, 'rot.mol')
    probe_t = translate(struc, probe_r, dist)
    #Chem.MolToMolFile(probe_t.MOL, 'transl.mol')
    probe_optimized = optimizeFit(struc, probe_t, rot, dist)
    return probe_optimized

def main(struc_file, probe_name, tip, tail, plane, dist, out):

    name, ext = struc_file.split('.')[0], struc_file.split('.')[-1]
    # define structure
    if tail:
        structure = ReadFile(ext, name, tip, tail=tail)
    elif plane:
        structure = ReadFile(ext, name, tip, plane=plane)
    else:
        print('Please supply tail/plane reference atom(s)')
        sys.exit()
    # define probe
    probe = Probe(probe_name)
    # dock probe
    docked = addProbe(structure, probe, dist=dist)
    # export docked structure as .mol
    ExportStructure(docked, out)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe addition Utility Package.',epilog='Add molecular probes to structures for SMART parameter generation.')
    parser.add_argument('-f', required=True) #structure file
    parser.add_argument('-o', required=False, default='SMART_probe_') #output file
    parser.add_argument('-tip', required=True) # tip atom (0-indexed)
    parser.add_argument('-tail', required=False) # tail atom (0-indexed)
    parser.add_argument('-plane', required=False) # tail plane atoms (2n, 0-indexed)
    parser.add_argument('-p', required=True) #probe name
    parser.add_argument('-dist', required=False, default=2.0) #probe distance

    args = parser.parse_args()

    starttime = time.time()
    main(args.f, args.p, args.tip, args.tail, args.plane, args.dist, args.o)
    print("--- %s seconds ---" % (time.time() - starttime))
