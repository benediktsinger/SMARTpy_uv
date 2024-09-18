#!/usr/bin/python
# SMART __init__()
#   version: b2.1 (released 03-14-2024)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#   GitHub:
############# ------- DEPENDENCIES
import sys, os, math, mathutils
import numpy as np
from scipy.spatial.transform import Rotation as R
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import ForceField
from rdkit import RDLogger
from rdkit.Geometry import Point3D
RDLogger.DisableLog('rdApp.*')
pt = Chem.GetPeriodicTable()

from SMART.utils import octahedral_, tetrahedral_, trigonal_planar_

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
        #print(self.U)
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()
        return self

class _crossproduct_vector:
    # calculate binding vectors
    def __init__(self, idPos, refPos1, refPos2, dist):
        try:
            self.TipPos = np.array(idPos)
            v1 = refPos1 - self.TipPos
            n1 = np.linalg.norm(v1)
            self.V1 = (v1 / n1)
            v2 = refPos2 - self.TipPos
            n2 = np.linalg.norm(v2)
            self.V2 = (v2 / n2)
            self.U = np.cross(self.V1, self.V2)
            self.BindingPos = np.array(self.TipPos) + (float(dist)*self.U)
            print(self.U)
        except np.linalg.LinAlgError:
            print('invalid coordinates')
            sys.exit()
        return self

############# ------- VECTOR INIT CLASSES
class Reference_Vector(_reference_vector):
    def __init__(self, tip_id=None, tip_pos=None, ref_id=None, ref_pos=None, dist=1.0):
        '''Reference_Vector: probe binding vector based on tail/tip vector definition'''
        if not tip_pos is None:
            self.TipPos = tip_pos
        elif tip_id != None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            print(tip_id)
            sys.exit()

        if not ref_pos is None:
            self.RefPos = ref_pos
        elif ref_id != None:
            self.RefId = ref_id
            self.RefPos = self.MOL.GetConformer().GetAtomPosition(int(self.RefId))
        else:
            print('Please supply ref_id or ref_pos')
            sys.exit()

        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

class Reference_Angle(_crossproduct_vector):
    def __init__(self, tip_id=None, tip_pos=None, refs_id=None, refs_pos=None, dist=1.0):
        '''Reference_Angle: cross product of input angle as binding reference vector'''
        if not tip_pos is None:
            self.TipPos = tip_pos
        elif tip_id != None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            sys.exit()

        if not refs_pos is None:
            if len(refs_pos) == 2:
                tail1_pos = refs_pos[0]
                tail2_pos = refs_pos[1]
            else:
                print('Please supply 2 ref points')
        elif refs_id != None:
            if len(refs_id) == 2:
                tail1_pos = self.MOL.GetConformer().GetAtomPosition(int(refs_id[0]))
                tail2_pos = self.MOL.GetConformer().GetAtomPosition(int(refs_id[1]))
            else:
                print('Please supply 2 ref points')
        else:
            print('Please supply refs_id or refs_pos')
            sys.exit()
        self.Vector= _crossproduct_vector.__init__(self, self.TipPos, tail1_pos, tail2_pos, dist)

class Define_Geometry(_reference_vector):
    def __init__(self, tip_id=None, tip_pos=None, refs_ids=None, refs_pos=None, geom=None, dist=1.0):
        '''Define_Geometry: define binding center geometry and reference atoms'''
        self.Geometry = geom
        self.RefsPos = refs_pos
        if not tip_pos is None:
            self.TipPos = tip_pos
        elif tip_id != None:
            self.TipId = tip_id
            self.TipPos = self.MOL.GetConformer().GetAtomPosition(int(self.TipId))
        else:
            print('Please supply tip_id or tip_pos')
            sys.exit()

        if not self.RefsPos:
            print('please specify binding atoms or use Detect_Geometry feature')
            print(refs_pos)
            sys.exit()

        if self.Geometry == 'octahedral':
            self.RefsPos = self.octahedral(dist)
        elif self.Geometry == 'tetrahedral':
            self.RefsPos = self.tetrahedral(dist)
        elif self.Geometry ==  'trigonal':
            self.RefsPos = self.trigonal_planar(dist)
        else:
            print('please specify a known geometry: trigonal, tetrahedral, octahedral')
            print(self.Geometry)
            sys.exit()

    def octahedral(self, dist):
        if not len(self.RefsPos) == 5:
            print('Cannot compute octahedral geometry from these atoms:')
            print(self.RefsPos)
            sys.exit()
        self.RefPos = octahedral_(self.MOL, self.TipPos, self.RefsPos)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

    def tetrahedral(self, dist):
        if not len(self.RefsPos) == 3:
            print('Cannot compute tetrahedral geometry from these atoms:')
            print(sself.RefsPos)
            sys.exit()
        self.RefPos = tetrahedral_(self.MOL, self.RefsPos)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

    def trigonal_planar(self, dist):
        if not len(self.RefsPos) == 2:
            print('Cannot compute trigonal planar geometry from these atoms:')
            print(self.RefsPos)
            sys.exit()
        self.RefPos = trigonal_planar_(self.RefsPos)
        self.Vector = _reference_vector.__init__(self, self.TipPos, self.RefPos, dist)

        #def trigonal_bipyramidal(self):

class Detect_Geometry(Define_Geometry):
    def __init__(self, tip_id=None, covalent=True, searchRadius=None):
        '''Detect_Geometry: automatically detect binding center geometry by binding radii '''
        print('Detect_Geometry: currently not availiable')
        sys.exit()

        if tip_id != None:
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

    def reference_geometry(self, tip_id=None, tip_pos=None, refs_id=None, refs_pos=None, geom=None, dist=0.0):
        Define_Geometry.__init__(self, tip_id, tip_pos, refs_id, refs_pos, geom, dist)

    def detect_geometry(self, tip_id=None, tip_pos=None, covalent=True, searchRadius=None):
        Detect_Geometry.__init__(self, tip_id, tip_pos, covalent, searchRadius)

    def export_alignment(self, out='out_align', dummy='X'):
        Export_Alignment.__init__(self, out, dummy)

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

class Export_Alignment:
    def __init__(self, out, dummy):
        Chem.MolToXYZFile(self.MOL, out+'.xyz')
        with open(out+'.xyz', 'a+') as f:
            f.write('\t'.join([dummy]+[str(i) for i in self.Vector.BindingPos.tolist()]))
            #f.write('\t'.join(['X']+list(self.Vector.BindingPos)))

class ExportStructure:
    def __init__(self, mol, outname='out'):
        '''ExportStructure: export structure as mol file'''
        writer = Chem.SDWriter(outname+'.sdf')
        for cid in [i.GetId() for i in mol.GetConformers()]:
            #print(mol.GetConformer(cid).GetAtomPosition(0))
            writer.write(mol, confId=cid)
        writer.close()
