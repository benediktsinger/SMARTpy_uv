#!/usr/bin/python
# SMART descriptor generators
#   version: b0.1 (released 10-10-2023)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#
############# ------- DEPENDENCIES
import math, argparse, time, sys, os
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFreeSASA, Descriptors3D
#from morfeus import BuriedVolume
import dbstep.Dbstep as db
import pyvista as pv

class read_probe_conformers():
    def __init__(self, file, path=os.getcwd()):
        try:
            self.NAME = file.split('.')[0]
            self.TYPE = file.split('.')[-1]
            file = os.path.join(path, file)
            suppl = Chem.SDMolSupplier(file, removeHs = False)
            self.CONFS = None
            for mol in suppl:
                if not self.CONFS:
                    self.CONFS = mol
                else:
                    self.CONFS.AddConformer(mol.GetConformer(), assignId=True)
        except Exception as e:
            print('no structure file found -', file)
            print(e)
            sys.exit()

############# ------- MAIN FNXS
def SplitFragments(confs):
    # separate probe and struc
    frags = Chem.GetMolFrags(confs, asMols=True)
    probe = frags[0]
    structure = Chem.DeleteSubstructs(confs, frags[0])

    return structure, probe

def CoordsFromMols(confs):
    # get atomic coords
    coords = []
    e = []

    for cid in range(confs.GetNumConformers()):
        for i in range(confs.GetNumAtoms()):
            atom = confs.GetAtomWithIdx(i)
            elem = atom.GetSymbol()
            pos = confs.GetConformer(cid).GetAtomPosition(i)
            e.append(elem)
            coords.append(np.array(pos))

    return coords, e

def DBSTEP_Properties(confs, m_id):
    structure, probe = SplitFragments(confs)

    properties = {}
    coords, elems = CoordsFromMols(probe)
    coords.append(np.array(confs.GetConformer().GetAtomPosition(int(m_id))))
    elems.append(confs.GetAtomWithIdx(int(m_id)).GetSymbol())

    max_radius = 0
    TOL = 1
    for i in range(len(coords)):
        dist = math.dist(coords[i], coords[-1])
        if dist > max_radius:
            max_radius = dist

    sterimol = db.dbstep(confs, atom1=int(m_id)+1+probe.GetNumAtoms(), atom2=1, sterimol=True)
    bv = db.dbstep(probe, atom1=1, volume=True)
    bv_max = db.dbstep(probe, atom1=1, volume=True, r=max_radius+TOL)
    #print(((4/3)*np.pi*((3.50)**3)))
    #print(((4/3)*np.pi*((max_radius+TOL)**3)))
    properties['Sterimol-L'] = sterimol.L #prox probe vol
    properties['Sterimol-B1'] = sterimol.Bmin #dist probe vol
    properties['Sterimol-B5'] = sterimol.Bmax #free sphere vol

    properties['%VBurPROXcavity'] = bv.bur_vol #prox probe vol
    properties['VBurPROXcavity'] = (bv.bur_vol/100)*((4/3)*np.pi*((3.50)**3)) #prox probe vol
    properties['VBurFULLcavity'] = (bv_max.bur_vol/100)*((4/3)*np.pi*((max_radius+TOL)**3)) #full probe vol

    #properties['VBurDISTcavity'] = bv_max.bur_vol-bv.bur_vol  #dist probe vol
    #properties['VBurFULL'] = bv_max.free_volume #free sphere vol

    return properties

def Morfeus_Properties(confs, m_id):
    structure, probe = SplitFragments(confs)

    properties = {}
    coords, elems = CoordsFromMols(confs)
    coords.append(np.array(confs.GetConformer().GetAtomPosition(int(m_id))))
    elems.append(confs.GetAtomWithIdx(int(m_id)).GetSymbol())

    max_radius = 0
    TOL = 2
    for i in range(len(coords)):
        dist = math.dist(coords[i],coords[0])
        if dist > max_radius:
            max_radius = dist


    #sterimol = mf.Sterimol(elems, coords, len(coords)-1, 0)
    bv = BuriedVolume(elems, coords, len(coords)-1).compute_distal_volume(method="buried_volume")
    bv_max = BuriedVolume(elems, coords, len(coords)-1, radius=max_radius+TOL).compute_distal_volume(method="buried_volume")

    properties['VBurPROXcavity'] = bv.buried_volume #prox probe vol
    properties['VBurDISTcavity'] = bv.distal_volume #dist probe vol
    #properties['VBurFULL'] = bv_max.free_volume #free sphere vol
    properties['VBurFULLcavity'] = bv_max.buried_volume #full probe vol

    return properties

def PyVista_Properties(confs, m_id):
    structure, probe = SplitFragments(confs)

    properties = {}
    pcoords, pelems = CoordsFromMols(confs)
    #scoords, selems = CoordsFromMols(structure)

    # gen point cloud from 3D
    p_ptcloud = pv.PolyData(pcoords)
    p_ptcloud.cast_to_pointset()
    #s_ptcloud = pv.PolyData(scoords)
    #s_ptcloud.cast_to_pointset()

    cloud = p_ptcloud.delaunay_3d(alpha=5)
    cloud.extract_geometry()

    properties['pyVista_VOLcavity'] =  cloud.volume
    properties['pyVista_AREAcavity'] = cloud.area

    return properties

def RDKit_Properties(confs, m_id):
    #compute rdkit descriptors
    structure, probe = SplitFragments(confs)

    properties = {}
    tradii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in confs.GetAtoms()]
    pradii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in probe.GetAtoms()]
    sradii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in structure.GetAtoms()]

    #pcoords, pelems = CoordsFromMols(probe)
    #scoords, selems = CoordsFromMols(structure)
    #mcoords, melems = CoordsFromMols(confs)

    properties['SASAcavity'] = rdFreeSASA.CalcSASA(probe, pradii)
    properties['SASAcat'] = rdFreeSASA.CalcSASA(structure, sradii)
    properties['SASAcomplx'] = rdFreeSASA.CalcSASA(confs, tradii)
    properties['CSA'] = ((properties['SASAcat']+properties['SASAcavity'])-properties['SASAcomplx'])/2
    properties['ESA'] = properties['SASAcavity']-properties['CSA']

    properties['Asphericity'] = Descriptors3D.Asphericity(probe)
    properties['Eccentricity'] = Descriptors3D.Eccentricity(probe)
    properties['InertialShapeFactor'] = Descriptors3D.InertialShapeFactor(probe)
    properties['SpherocityIndex'] = Descriptors3D.SpherocityIndex(probe)
    properties['NormalizedInertiaRatio1/3'] = Descriptors3D.NPR1(probe)

    return properties

def main(file, m_id, rdkit, dbstep):

    dfs = {}

    confs = read_probe_conformers(file)

    if rdkit:
        print('\ncomputing RDKit descriptors')
        properties = RDKit_Properties(confs.CONFS, m_id)

        dfs['rdkit_SASAcavity'] = properties['SASAcavity']
        dfs['rdkit_SASAcat']= properties['SASAcat']
        dfs['rdkit_SASAcomplx'] = properties['SASAcomplx']
        dfs['rdkit_CSA'] = properties['CSA']
        dfs['rdkit_ESA'] = properties['ESA']
        dfs['rdkit_Asphericity'] = properties['Asphericity']
        dfs['rdkit_Eccentricity'] = properties['Eccentricity']
        dfs['rdkit_InertialShapeFactor'] = properties['InertialShapeFactor']
        dfs['rdkit_SpherocityIndex'] = properties['SpherocityIndex']
        dfs['rdkit_NormalizedInertiaRatio1/3'] = properties['NormalizedInertiaRatio1/3']
        print(dfs)

    if dbstep:
        print('\ncomputing DBSTEP descriptors')
        properties = DBSTEP_Properties(confs.CONFS, m_id)
        dfs['dbstep_Sterimol-L'] = properties['Sterimol-L']
        dfs['dbstep_Sterimol-B1'] = properties['Sterimol-B1']
        dfs['dbstep_Sterimol-B5'] = properties['Sterimol-B5']
        dfs['dbstep_%VBurPROXcavity'] = properties['%VBurPROXcavity']
        dfs['dbstep_VBurPROXcavity'] = properties['VBurPROXcavity']
        #dfs['dbstep_VBurDISTcavity'] = properties['VBurDISTcavity']
        #dfs['dbstep_VBurFULL'] = properties['VBurFULL']
        dfs['dbstep_VBurFULLcavity'] = properties['VBurFULLcavity']
        print(dfs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Probe addition Utility Package.',epilog='Add molecular probes to structures for SMART parameter generation.')
    parser.add_argument('-f', required=True) #probe conformer file
    parser.add_argument('-m', required=True) #metal ID
    #parser.add_argument('-o', required=False, action='store_true') #save output
    parser.add_argument('--pyvista', required=False, action='store_true', default=False)
    parser.add_argument('--rdkit', required=False, action='store_true', default=False)
    parser.add_argument('--dbstep', required=False, action='store_true', default=False)

    args = parser.parse_args()

    starttime = time.time()
    main(args.f, args.m, args.rdkit, args.dbstep)
    print("--- %s seconds ---" % (time.time() - starttime))
