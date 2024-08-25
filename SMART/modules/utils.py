#!/usr/bin/python
# SMART probe addition utility
#   version: b2.1 (released 03-14-2024)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#   GitHub:
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

############# ------- VECTOR UTILITY FUNCTIONS
def _octahedral(mol, tip_pos, bindingAtoms):
    TOL = 20
    for a in bindingAtoms:
        n = 0
        min, max, avg = 360, 0, 0
        for b in bindingAtoms:
            if a == b:
                continue
            n += 1
            angle = rdMolTransforms.GetAngleDeg(mol.GetConformer(), a, id, b)
            avg += angle
            if angle > max:
                max = angle
            if angle < min:
                min = angle
        if max < (180-TOL):
            return np.array(mol.GetConformer().GetAtomPosition(a))
            #return a-1

def _tetrahedral(mol, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[0])) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[1])) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[2])) ) / 3
    return ref_pos

def _trigonal_planar(mol, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(int(bindingAtoms[0]))) +
                np.array(mol.GetConformer().GetAtomPosition(int(bindingAtoms[1]))) ) / 2
    return ref_pos

############# ------- PROTECTED VECTOR CONSTRUCTION CLASSES
class _reference_vector:
    def __init__(self, idPos, refPos, dist):
        try:
            self.TipPos = np.array(idPos)
            self.RefPos = np.array(refPos)
            v = (idPos - refPos)
            n = np.linalg.norm(v)
            self.U = np.array(v / n)
            self.BindingPos = np.array(self.TipPos) + (float(dist)*self.U)
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()
        return self

class _crossproduct_vector:
    # calculate binding vectors
    def __init__(self, idPos, refPos1, refPos2, dist):
        try:
            self.TipPos = idPos
            v1 = refPos1 - self.TipPos
            n1 = np.linalg.norm(v1)
            self.V1 = (v1 / n1)
            v2 = refPos2 - self.TipPos
            n2 = np.linalg.norm(v2)
            self.V2 = (v2 / n2)
            self.U = np.cross(self.V1, self.V2)
            self.BindingPos = np.array(self.TipPos) + (float(dist)*np.array(self.U))
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()
        #return self

############# ------- VECTOR INIT CLASSES
class Reference_Vector(_reference_vector):
    def __init__(self, tip_id=None, ref_id=None, tip_pos=None, ref_pos=None, dist=0.0):
        '''Reference_Vector: probe binding vector based on tail/tip vector definition'''
        if tip_pos:
            self.TipPos = tip_pos
        elif tip_id != None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            print(tip_id)
            sys.exit()

        if ref_pos:
            self.RefPos = tip_pos
        elif ref_id != None:
            self.RefId = ref_id
            self.RefPos = self.MOL.GetConformer().GetAtomPosition(int(self.RefId))
        else:
            print('Please supply ref_id or ref_pos')
            sys.exit()

        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)


class Reference_Angle(_crossproduct_vector):
    def __init__(self, tip_id=None, tip_pos=None, refs_id=None, refs_pos=None, dist=0.0):
        '''Reference_Angle: cross product of input angle as binding reference vector'''
        if tip_pos:
            self.TipPos = tip_pos
        elif tip_id != None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            sys.exit()

        if refs_pos:
            if len(refs_pos) == 2:
                tail1_pos = refs_pos[0]
                tail2_pos = refs_pos[1]
            else:
                print('Please supply 2 ref points')
        elif refs_id!= None:
            if len(refs_id) == 2:
                tail1_pos = self.MOL.GetConformer().GetAtomPosition(int(refs_id[0]))
                tail2_pos = self.MOL.GetConformer().GetAtomPosition(int(refs_id[1]))
            else:
                print('Please supply 2 ref points')
        else:
            print('Please supply refs_id or refs_pos')
            sys.exit()
        self.Vector= _crossproduct_vector.__init__(self, self.TipPos, tail1_pos, tail2_pos)

class Define_Geometry(_reference_vector):
    def __init__(self, tip_id=None, tip_pos=None, bindingAtoms=None, geom=None, dist=0.0):
        '''Define_Geometry: define binding center geometry and reference atoms'''
        self.Geometry = geom
        self.RefAtoms = bindingAtoms
        if tip_pos:
            self.TipPos = tip_pos
        elif tip_id!= None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            sys.exit()

        if not self.RefAtoms:
            print('please specify binding atoms or use Detect_Geometry feature')
            print(bindingAtoms)
            sys.exit()

        if self.Geometry == 'octahedral':
            self.RefPos = self.octahedral(dist)
        elif self.Geometry == 'tetrahedral':
            self.RefPos = self.tetrahedral(dist)
        elif self.Geometry ==  'trigonal':
            self.RefPos = self.trigonal_planar(dist)
        else:
            print('please specify a known geometry: trigonal, tetrahedral, octahedral')
            print(self.Geometry)
            sys.exit()

    def octahedral(self, dist):
        if not len(self.RefAtoms) == 5:
            print('Cannot compute octahedral geometry from these atoms:')
            print(self.RefAtoms)
            sys.exit()
        self.RefPos = _octahedral(self.MOL, self.TipPos, self.RefAtoms)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

    def tetrahedral(self, dist):
        if not len(self.RefAtoms) == 3:
            print('Cannot compute tetrahedral geometry from these atoms:')
            print(self.RefAtoms)
            sys.exit()
        self.RefPos = _tetrahedral(self.MOL, self.RefAtoms)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

    def trigonal_planar(self, dist):
        if not len(self.RefAtoms) == 2:
            print('Cannot compute trigonal planar geometry from these atoms:')
            print(self.RefAtoms)
            sys.exit()
        self.RefPos = _trigonal_planar(self.MOL, self.RefAtoms)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

        #def trigonal_bipyramidal(self):

class Detect_Geometry(Define_Geometry):
    def __init__(self, tip_id=None, covalent=True, searchRadius=None):
        '''Detect_Geometry: automatically detect binding center geometry by binding radii '''
        print('Detect_Geometry: currently not availiable')
        sys.exit()

        if tip_id!= None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id')
            sys.exit()

        atom = self.MOL.GetAtomWithIdx(tip_id)
        if covalent:
            self.RefAtoms = atom.GetNeighbors()
        else:
            self.RefAtoms = []
            atom_ = self.MOL.GetConformer().GetAtomPosition(tip_id)
            atom_r = Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
            for i in self.MOL.GetAtoms():
                i_ = self.MOL.GetConformer().GetAtomPosition(i.GetIdx())
                i_r = Chem.GetPeriodicTable().GetRvdw(i.GetAtomicNum())
                d = (i_-atom_).LengthSq()
                r_d = i_r + atom_r
                if d <= r_d:
                    self.RefAtoms.append(i.GetIdx())

        if len(neighbors) == 2:
            Define_Geometry.trigonal_planar(self.MOL, self.RefAtoms)
        elif len(neighbors) == 3:
            Define_Geometry.tetrahedral(self.MOL, self.RefAtoms)
        #elif len(neighbors) == 4:
        #    Define_Geometry.trigonal_bipyramidal(self)
        elif len(neighbors) == 5:
            Define_Geometry.octahedral(self.MOL, self.TipPos, self.RefAtoms)
        else:
            print('Unknown geometry. Try covalent=False or Define_Geometry')

############# ------- INPUT/OUTPUT UTILITY CLASSES
class PARAMS():
    def __init__(self):
        pass
PARAMS.CLASHTOL = 0.0

class ReadFile(Reference_Vector, Detect_Geometry, Reference_Angle):
    # initialize single structure coordinates
    def __init__(self, fname, path=os.getcwd()):
        '''initialize ReadFile: file.ext as input and base attributes'''
        try:
            os.path.isfile(os.path.join(path, fname))
        except Exception as e:
            print('no structure file found -', path, fname)
            print(e)
            sys.exit()

        self.Name = fname
        self.Type = fname.split('.')[-1]
        self.File = os.path.join(path, fname)

        try:
            if self.Type == 'mol2':
                self.MOL = Chem.MolFromMol2File(self.File, removeHs=False, sanitize=False)
            elif self.Type == 'mol':
                self.MOL = Chem.MolFromMolFile(self.File, removeHs=False, sanitize=False)
            elif self.Type == 'xyz':
                self.MOL = Chem.MolFromXYZFile(self.File, removeHs=False, sanitize=False)
            elif self.Type == 'pdb':
                self.MOL = Chem.MolFromPDBFile(self.File, removeHs=False, sanitize=False)
            elif self.Type == 'sdf':
                self.MOL = Chem.SDMolSupplier(self.File, removeHs=False, sanitize=False)[0]
            self.NumAtoms = self.MOL.GetNumAtoms()
        except Exception as e:
            print('no structure in file -', fname)
            print(e)
            sys.exit()

    def reference_vector(self, tip_id=None, tip_pos=None, ref_id=None, ref_pos=None, dist=0.0):
        Reference_Vector.__init__(self, tip_id, tip_pos, ref_id, ref_pos, dist)

    def reference_angle(self, tip_id=None, tip_pos=None, refs_id=None, refs_pos=None, dist=0.0):
        Reference_Angle.__init__(self, tip_id, tip_pos, refs_id, refs_pos, dist)

    def reference_geometry(self, tip_id=None, tip_pos=None, bindingAtoms=None, geom=None, dist=0.0):
        Define_Geometry.__init__(self, tip_id, tip_pos, bindingAtoms, geom, dist)

    def detect_geometry(self, tip_id=None, tip_pos=None, covalent=True, searchRadius=None):
        Detect_Geometry.__init__(self, tip_id, tip_pos, covalent, searchRadius)

class ReadMol(ReadFile):
    # initialize RDKit MOL as input
    def __init__(self, mol):
        '''ReadMol: RDKit MOL as input and base attributes'''
        self.MOL = mol
        self.NumAtoms = self.MOL.GetNumAtoms()
        self.Name = None

class ReadProbe(Reference_Vector):
    # initialize probe
    def __init__(self, probe, path=None):
        '''ReadProbe: initialize probe from .mol2 file'''
        if path:
            try:
                PROBE_PATH = path
                os.path.isfile(os.path.join(PROBE_PATH, probe))
            except Exception as e:
                print('no probe file found -', path, probe)
                print(e)
                sys.exit()
        else:
            try:
                PROBE_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Probes')
                os.path.isfile(os.path.join(PROBE_PATH, probe))
            except Exception as e:
                print('no probe file found -', PROBE_PATH, probe)

        if not probe.endswith('.mol2'):
            file = os.path.join(PROBE_PATH, probe+'.mol2')
        else:
            file = os.path.join(PROBE_PATH, probe)
        self.Name = probe
        self.MOL = Chem.MolFromMol2File(file, removeHs=False)
        self.NumAtoms = self.MOL.GetNumAtoms()
        self.TipId = 0
        self.RefId = self.NumAtoms-1
        Reference_Vector.__init__(self, tip_id=self.TipId, tip_pos=None, ref_id=self.RefId, ref_pos=None)

class ExportStructure:
    def __init__(self, mol, outname='out'):
        '''ExportStructure: export structure as mol file'''
        Chem.MolToPDBFile(mol, outname+'.pdb')

############# ------- MAIN FNXS
def _clashCheck(conf, struc, probe):
    # check for probe-structure clashes
    for p in range(probe.NumAtoms):
        i_ = np.array(conf.GetConformer().GetAtomPosition(p))
        if (i_ == probe.Vector.TipPos).all():
            continue
        for s in range(probe.NumAtoms,probe.NumAtoms+struc.NumAtoms):
            j_ = np.array(conf.GetConformer().GetAtomPosition(s))
            if (j_ == struc.Vector.TipPos).all() or (j_ == struc.Vector.RefPos).all():
                continue
            i_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(p).GetAtomicNum())
            j_r = Chem.GetPeriodicTable().GetRvdw(conf.GetAtomWithIdx(s).GetAtomicNum())
            d = np.linalg.norm(i_-j_)
            if d <= 4:
                r_d = i_r + j_r + PARAMS.CLASHTOL
                if d <= r_d:
                    return True
    return False

def _optimizeFit(struc, probe, incr, theta=5):
    # rotate probe to fit structure
    conf = Chem.CombineMols(probe.MOL, struc.MOL)
    conf_ = Chem.RWMol(conf)
    conf_.RemoveAtom(probe.RefId)
    AllChem.SanitizeMol(conf_)
    return conf_

    if not _clashCheck(conf_, struc, probe):
        return conf_
    else:
        if theta < 360:
            theta += incr
            probe_r = _rotation(struc, probe, probe.Vector.U, math.radians(theta))
            #probe_r = _rotation2(struc, probe, np.array(struc.Vector.TipPos), probe.Vector.U, math.radians(theta))
            return _optimizeFit(struc, probe_r, incr, theta)
        else:
            print('failed:', struc.Name)
            return conf_
    return conf_

def _rotateAlign(struc, probe):
    r = R.align_vectors([struc.Vector.U], [probe.Vector.U])
    r[0].as_matrix()
    newv = r[0].apply([np.array(probe.MOL.GetConformer().GetAtomPosition(i)) for i in range(probe.NumAtoms)])

    for i in range(probe.NumAtoms):
        probe.MOL.GetConformer().SetAtomPosition(i, newv[i])

    #probe.Vector = _reference_vector(probe.MOL.GetConformer().GetAtomPosition(probe.Id), probe.MOL.GetConformer().GetAtomPosition(probe.Ref))

    probe.Vector.TipPos = probe.Vector.TipPos
    if (probe.Vector.TipPos==struc.Vector.TipPos).all():
        probe.Vector.RefPos = struc.Vector.RefPos
    else:
        probe.Vector.RefPos = struc.Vector.TipPos
    v = probe.Vector.TipPos - probe.Vector.RefPos
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def _rotation2(struc, probe, tip, rotvec, theta):
    for i in range(probe.NumAtoms):
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

    probe.Vector.TipPos = probe.Vector.TipPos
    if (probe.Vector.TipPos==struc.Vector.TipPos).all():
        probe.Vector.RefPos = struc.Vector.RefPos
    else:
        probe.Vector.RefPos = struc.Vector.TipPos
    v = probe.Vector.TipPos - probe.Vector.RefPos
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def _rotation(struc, probe, rotvec, theta):
    # rotate 3D about vector
    quat = mathutils.Quaternion(rotvec, math.radians(theta))
    r = R.from_quat(quat)
    r.as_quat()
    newv = r.apply([np.array(probe.MOL.GetConformer().GetAtomPosition(i)) for i in range(probe.NumAtoms)])
    for i in range(probe.NumAtoms):
        probe.MOL.GetConformer().SetAtomPosition(i, newv[i])

    probe.Vector.TipPos = probe.Vector.TipPos
    if (probe.Vector.TipPos == struc.Vector.TipPos).all():
        probe.Vector.RefPos = struc.Vector.RefPos
    else:
        probe.Vector.RefPos = struc.Vector.TipPos
    v = probe.Vector.TipPos - probe.Vector.RefPos
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

def _translate(struc, probe):
    # translate to binding Point3D
    translation_v = struc.Vector.BindingPos - probe.Vector.TipPos

    for i in range(probe.NumAtoms):
        pt = np.array(probe.MOL.GetConformer().GetAtomPosition(i))
        translate = pt+translation_v
        probe.MOL.GetConformer().SetAtomPosition(i, translate)

    probe.Vector.TipPos = probe.Vector.TipPos
    if (probe.Vector.TipPos == struc.Vector.TipPos).all():
        probe.Vector.RefPos = struc.Vector.RefPos
    else:
        probe.Vector.RefPos = struc.Vector.TipPos
    v = probe.Vector.TipPos - probe.Vector.RefPos
    n = np.linalg.norm(v)
    probe.Vector.U = v / n

    return probe

############# ------- CTRL FNXS
def add_probe(struc, probe, rot=5.0):
    '''addProbe: compute probe translation and rotation before docking to structure'''
    # main utility function
    probe_r = _rotateAlign(struc, probe)
    #Chem.MolToMolFile(probe_r.MOL, 'rot.mol')
    probe_t = _translate(struc, probe_r)
    #Chem.MolToMolFile(probe_t.MOL, 'transl.mol')
    probe_optimized = _optimizeFit(struc, probe_t, rot)
    return probe_optimized

def main(struc_file, probe_name, tip, tail, plane, dist, out='SMART_probe_'):
    '''main: ctrl function for cmd line usage'''
    name, ext = struc_file.split('.')[0], struc_file.split('.')[-1]
    # define structure
    if tail:
        structure = ReadFile(ext, name, tip, tail=tail)
    elif plane:
        structure = ReadFile(ext, name, tip, plane=plane)
    else:
        print('Please supply reference atom(s)')
        sys.exit()
    # define probe
    probe = Probe(probe_name)
    # dock probe
    docked = addProbe(structure, probe, dist=dist)
    # export docked structure as .mol
    ExportStructure(docked, out)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe Addition Utility Package.',epilog='Add molecular probes to structures for SMART parameter generation.')
    parser.add_argument('-f', required=True) #structure file
    parser.add_argument('-p', required=True) #probe name
    parser.add_argument('-id', required=True) # binding (tip) atom (0-indexed)

    parser.add_argument('-ref', required=False) # tail atom (0-indexed)
    parser.add_argument('-angle', required=False) # angle reference atoms (2n, 0-indexed)
    parser.add_argument('-geom', required=False, choices=['tetrahedral','trigonal','octahedral','auto']) # geometry
    parser.add_argument('-dist', required=False, default=2.0) #probe distance
    parser.add_argument('-o', required=False, default='SMART_probe_') #output file

    args = parser.parse_args()

    starttime = time.time()
    main(args.f, args.p, args.tip, args.tail, args.plane, args.dist, args.o)
    print("--- %s seconds ---" % (time.time() - starttime))
