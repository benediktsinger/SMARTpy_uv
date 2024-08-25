#!/usr/bin/python
# SMART utility functions
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
def octahedral_(mol, tip_pos, bindingAtoms):
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

def tetrahedral_(mol, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[0])) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[1])) +
                np.array(mol.GetConformer().GetAtomPosition(bindingAtoms[2])) ) / 3
    return ref_pos

def trigonal_planar_(mol, bindingAtoms):
    ref_pos = ( np.array(mol.GetConformer().GetAtomPosition(int(bindingAtoms[0]))) +
                np.array(mol.GetConformer().GetAtomPosition(int(bindingAtoms[1]))) ) / 2
    return ref_pos


############# ------- CONF SEARCH UTILITY FNXS
def clash_check_(mol, probe, cid, CLASHTOL=0.0):
    for p in range(1, probe.GetNumAtoms()):
        for s in range(mol.GetNumAtoms()):
                i_ = probe.GetConformer(cid).GetAtomPosition(p)
                i_r = Chem.GetPeriodicTable().GetRvdw(probe.GetAtomWithIdx(p).GetAtomicNum())
                j_ = mol.GetConformer().GetAtomPosition(s)
                j_r = Chem.GetPeriodicTable().GetRvdw(mol.GetAtomWithIdx(s).GetAtomicNum())
                d = (i_-j_).LengthSq()
                r_d = i_r + j_r + CLASHTOL
                if (d - r_d) <= 0:
                    return True
    return False

def vectorize_(frag, tip, atoms, rotvec, theta):
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

def rotation_(mol, theta):
    # rotate 3D about vector
    cid_ = mol.GetConformers()[0].GetId()
    rotvec = mol.GetConformer(cid_).GetAtomPosition(0) - mol.GetConformer(cid_).GetAtomPosition(1) / np.linalg.norm(mol.GetConformer(cid_).GetAtomPosition(0) - mol.GetConformer(cid_).GetAtomPosition(1))
    quat = mathutils.Quaternion(rotvec, theta)
    r = R.from_quat(quat)
    r.as_quat()

    for cid in [conf.GetId() for conf in mol.GetConformers()]:
        newv = r.apply([np.array(mol.GetConformer(cid).GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        for i in range(mol.GetNumAtoms()):
            mol.GetConformer(cid).SetAtomPosition(i, newv[i])
    return mol

def translate_(struc, probe):
    # translate to binding Point3D
    translation_v = struc.Vector.BindingPos - probe.Vector.TipPos

    for i in range(probe.NumAtoms):
        pt = np.array(probe.MOL.GetConformer().GetAtomPosition(i))
        translate = pt+translation_v
        probe.MOL.GetConformer().SetAtomPosition(i, translate)

    return probe

############# ------- GENRAL UTILITY FNXS
def cont_to_mol_(mol, conf):
    id = conf.GetId()
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(mol.GetConformer(id))
    return new_mol
