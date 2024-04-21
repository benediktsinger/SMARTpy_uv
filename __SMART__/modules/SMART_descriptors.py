#!/usr/bin/python
# SMART descriptor generators
#   version: b2.1 (released 03-13-2024)
#   developer: Beck R. Miller (beck.miller@utah.edu)
#
############# ------- DEPENDENCIES
import math, argparse, time, sys, os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFreeSASA, Descriptors3D, AllChem

try:
    from morfeus import BuriedVolume, Sterimol, SASA
except Exception as e:
    print('morfeus is not installed')
try:
    import dbstep.Dbstep as db
except Exception as e:
    print('dbstep is not installed')
try:
    import pyvista as pv
except Exception as e:
    print('pyvista is not installed')

############# ------- UTILITY FNXS
def read_CAVITY_file(cavity_file):
    name = cavity_file.split('.')[0]
    ext = cavity_file.split('.')[1]
    cav = None
    if ext == 'mol2':
        cav = Chem.MolFromMol2File(file, removeHs=False)
    elif ext == 'mol':
        cav = Chem.MolFromMolFile(file, removeHs=False)
    elif ext == 'xyz':
        cav = Chem.MolFromXYZFile(file, removeHs=False)
    elif ext == 'pdb':
        cav = Chem.MolFromPDBFile(file, removeHs=False)
    elif ext == 'sdf':
        suppl = Chem.SDMolSupplier(file, removeHs = False)
        for m in suppl:
            if not mol:
                cav = m
            else:
                cav.AddConformer(m.GetConformer(), assignId=True)
    return cav

def read_STRUCTURE_file(mol_file):
    name = mol_file.split('.')[0]
    ext = mol_file.split('.')[1]
    mol = None
    if ext == 'mol2':
        mol = Chem.MolFromMol2File(file, removeHs=False)
    elif ext == 'mol':
        mol = Chem.MolFromMolFile(file, removeHs=False)
    elif ext == 'xyz':
        mol = Chem.MolFromXYZFile(file, removeHs=False)
    elif ext == 'pdb':
        mol = Chem.MolFromPDBFile(file, removeHs=False)
    elif ext == 'sdf':
        suppl = Chem.SDMolSupplier(file, removeHs = False)
        for m in suppl:
            if not mol:
                mol = m
            else:
                mol.AddConformer(m.GetConformer(), assignId=True)
    return mol

def read_CPLX_file(cplx_file):
    name = cplx_file.split('.')[0]
    ext = cplx_file.split('.')[1]
    cplx = None
    if ext == 'mol2':
        cplx = Chem.MolFromMol2File(file, removeHs=False)
    elif ext == 'mol':
        cplx = Chem.MolFromMolFile(file, removeHs=False)
    elif ext == 'xyz':
        cplx = Chem.MolFromXYZFile(file, removeHs=False)
    elif ext == 'pdb':
        cplx = Chem.MolFromPDBFile(file, removeHs=False)
    elif ext == 'sdf':
        suppl = Chem.SDMolSupplier(file, removeHs = False)
        for m in suppl:
            if not mol:
                cplx = m
            else:
                cplx.AddConformer(m.GetConformer(), assignId=True)
    return cplx

def _coords_from_mols(confs):
    ''' convert CONF object to xyz coordinate, element name lists '''
    #print('gen coords',confs.GetNumConformers())
    coords = []
    elems = []
    for c in confs.GetConformers():
        cid = c.GetId()
        for i, atom in enumerate(confs.GetAtoms()):
            #atom = confs.GetAtomWithIdx(i)
            pos = confs.GetConformer(cid).GetAtomPosition(i)
            coords.append(np.array(pos))
            elem = atom.GetSymbol()
            elems.append(elem)
    return coords, elems

def _coords_from_mol(struc):
    ''' convert MOL object to xyz coordinate, element name lists '''
    #print('gen coords',confs.GetNumConformers())
    coords = []
    elems = []
    for i, atom in enumerate(struc.GetAtoms()):
        #atom = struc.GetAtomWithIdx(i)
        pos = struc.GetConformer().GetAtomPosition(i)
        coords.append(np.array(pos))
        elem = atom.GetSymbol()
        elems.append(elem)
    return coords, elems

def _max_radius(coords, TOL=5):
    ''' compute max radius across list of coords (increased by TOL factor to account for radii) '''
    max_radius = 0
    for i in range(len(coords)):
        dist = math.dist(coords[i], coords[-1])
        if dist > max_radius:
            max_radius = dist
    return max_radius

def export_descriptors(properties, out):
    df = pd.DataFrame(properties)
    df.to_excel(out+'.xlsx')

############# ------- MAIN FNXS
def DBSTEP_Properties(struc, probe, id, prox_radius=3.5):
    ''' Uses dbstep package to compute properties by buried volume and sterimol '''

    print('\nComputing DBSTEP descriptors')
    properties = {}
    coords, elems = _coords_from_mols(probe)

    coords.append(np.array(struc.GetConformer().GetAtomPosition(int(id))))
    elems.append(struc.GetAtomWithIdx(int(id)).GetSymbol())

    max_radius = _max_radius(coords, TOL=5)

    sterimol = db.dbstep(probe, atom1=int(id)+probe.GetNumAtoms(), atom2=0, sterimol=True)
    bv = db.dbstep(probe, atom1=0, volume=True, r=prox_radius)
    bv_max = db.dbstep(probe, atom1=0, volume=True, r=max_radius)

    # VBur cavity properties
    properties['dbstep_VBur_VOL_cavity'] = (bv_max.bur_vol/100)*((4/3)*np.pi*((max_radius)**3))
    properties['dbstep_VBur_%VOL_cavity'] = bv_max.bur_vol
    properties['dbstep_VBur_proxVOL_cavity'] = (bv.bur_vol/100)*((4/3)*np.pi*((prox_radius)**3))
    properties['dbstep_VBur_%proxVOL_cavity'] = bv.bur_vol
    properties['dbstep_VBur_distVOL_cavity'] = properties['dbstep_VBur_VOL_cavity'] - properties['dbstep_VBur_proxVOL_cavity']
    # Sterimol cavity properties
    properties['dbstep_Sterimol_L_cavity'] = sterimol.L
    properties['dbstep_Sterimol_B1_cavity'] = sterimol.Bmin
    properties['dbstep_Sterimol_B5_cavity'] = sterimol.Bmax

    return properties

def Morfeus_Properties(struc, probe, id, prox_radius=3.5):
    ''' Uses morfeus package to compute properties by buried volume, sterimol, and sasa '''

    print('\nComputing Morfeus descriptors\n')
    properties = {}
    pcoords, pelems = _coords_from_mols(probe)
    scoords, selems = _coords_from_mols(struc)
    max_radius = _max_radius(pcoords, TOL=5)

    bv = BuriedVolume(pelems, pcoords, 0, radius=prox_radius).compute_distal_volume(method="buried_volume")
    bv_max = BuriedVolume(pelems, pcoords, 0, radius=max_radius).compute_distal_volume(method="buried_volume")

    probe_sasa = SASA(pelems, pcoords)
    struc_sasa = SASA(selems, scoords)
    cmplx_sasa = SASA(pelems+selems, pcoords+scoords)

    #print('Sterimol')
    #pcoords.append(np.array(struc.GetConformer().GetAtomPosition(int(id))))
    #pelems.append(struc.GetAtomWithIdx(int(id)).GetSymbol())
    #sterimol = Sterimol(pelems, pcoords, len(pcoords)-1, 0)

    # SASA cavity properties
    properties['morfeus_SASA_AREA_cavity'] = probe_sasa.area
    properties['morfeus_SASA_CSA_cavity'] = ( (probe_sasa.area + struc_sasa.area) - cmplx_sasa.area) / 2
    properties['morfeus_SASA_ESA_cavity'] = probe_sasa.area - properties['morfeus_SASA_CSA_cavity']
    properties['morfeus_SASA_%ESA_cavity'] = 100 * ( properties['morfeus_SASA_ESA_cavity'] / probe_sasa.area )
    # VBur cavity properties
    properties['morfeus_VBur_VOL_cavity'] = bv_max.buried_volume
    properties['morfeus_VBur_%VOL_cavity'] = bv_max.percent_buried_volume
    properties['morfeus_VBur_proxVOL_cavity'] = bv.buried_volume
    properties['morfeus_VBur_%proxVOL_cavity'] = bv.percent_buried_volume
    properties['morfeus_VBur_distVOL_cavity'] = bv.distal_volume
    # Sterimol cavity properties
    #properties['morfeus_Sterimol_L_cavity'] = sterimol.L_value
    #properties['morfeus_Sterimol_B1_catvity'] = sterimol.B_1_value
    #properties['morfeus_Sterimol_B5_cavity'] = sterimol.B_5_value

    return properties

def PyVista_Properties(struc, probe, id, prox_radius=3.5, a=0):
    ''' Uses pyvista package to assemble geometric objects and extract properties '''

    print('\nComputing PyVista descriptors\n')
    properties = {}
    pcoords, pelems = _coords_from_mols(probe)
    scoords, selems = _coords_from_mols(struc)
    pt = scoords[id]

    # gen point cloud from 3D coords
    p_ptcloud = pv.PolyData(pcoords)
    p_ptcloud.cast_to_pointset()
    s_ptcloud = pv.PolyData(scoords)
    s_ptcloud.cast_to_pointset()
    cmplx_ptcloud = pv.PolyData(pcoords+scoords)
    cmplx_ptcloud.cast_to_pointset()

    # make shape object from pt cloud (delaunay triangulation)
    cavity = p_ptcloud.delaunay_3d(alpha=a)
    cavity = cavity.extract_geometry()
    cavity = cavity.extract_surface().triangulate()
    struc = s_ptcloud.delaunay_3d(alpha=a)
    struc = struc.extract_geometry()
    struc = struc.extract_surface().triangulate()
    cmplx = cmplx_ptcloud.delaunay_3d(alpha=a)
    cmplx = cmplx.extract_geometry()
    cmplx = cmplx.extract_surface().triangulate()

    # proximal/distal clipping
    sphere = pv.Sphere(radius=prox_radius, center=pt)
    prox = p_ptcloud.clip_surface(sphere, invert=True)
    prox = prox.delaunay_3d(alpha=a)
    prox = prox.extract_geometry()
    prox = prox.extract_surface().triangulate()

    #p = pv.Plotter()
    #p.add_mesh(p_ptcloud, color='gold', label='full', opacity=0.75)
    #p.add_mesh(sphere, color='black', label='full', opacity=0.75)
    #p.add_mesh(prox, color='blue', label='prox', opacity=0.75)
    #p.add_mesh(dist, color='red', label='dist', opacity=0.75)
    #p.add_legend()
    #p.show()

    #print('cavity', cavity.volume, cavity.area)
    #print('struc',  struc.volume, struc.area)
    #print('cmplx',  cmplx.volume, cmplx.area)
    #print('proximal', prox.volume, prox.area)

    # full cavity properties
    properties['pyVista_delaunay_VOL_cavity'] =  float(cavity.volume)
    properties['pyVista_delaunay_AREA_cavity'] = float(cavity.area)
    properties['pyVista_delaunay_CSA_cavity'] =  ( ( float(struc.area) + float(cavity.area) ) - float(cmplx.area) ) / 2
    properties['pyVista_delaunay_ESA_cavity'] = float(cavity.area) - properties['pyVista_delaunay_CSA_cavity']
    properties['pyVista_delaunay_%ESA_cavity'] = 100 * ( properties['pyVista_delaunay_ESA_cavity'] / float(cavity.area) )
    # prox/dist cavity properties
    properties['pyVista_delaunay_proxVOL_cavity'] =  prox.volume
    properties['pyVista_delaunay_proxAREA_cavity'] = prox.area
    properties['pyVista_delaunay_distVOL_cavity'] =  cavity.volume - prox.volume
    #properties['pyVista_delaunay_distAREA_cavity'] = cavity.area - prox.area

    return properties

def RDKit_Properties(struc, probe, cmplx):
    ''' Uses rdkit package to compute properties '''

    print('\nComputing RDKit descriptors\n')
    properties = {}
    probe_radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in probe.GetAtoms()]
    struc_radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in struc.GetAtoms()]
    cmplx_sradii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in cmplx.GetAtoms()]

    area_struc = rdFreeSASA.CalcSASA(struc, struc_radii)
    area_cmplx = rdFreeSASA.CalcSASA(cmplx, cmplx_sradii)

    #properties['RDKit_Mol_VOL_cavity'] = AllChem.ComputeMolVolume(probe)
    # SASA properties
    properties['RDKit_SASA_AREA_cavity'] = rdFreeSASA.CalcSASA(probe, probe_radii)
    #print('computed sasa')
    properties['RDKit_SASA_CSA_cavity'] = ( ( area_struc + properties['RDKit_SASA_AREA_cavity'] ) - area_cmplx ) / 2
    properties['RDKit_SASA_ESA_cavity'] = properties['RDKit_SASA_AREA_cavity'] - properties['RDKit_SASA_CSA_cavity']
    properties['RDKit_SASA_%ESA_cavity'] = 100 * ( properties['RDKit_SASA_ESA_cavity'] / properties['RDKit_SASA_AREA_cavity'] )
    # Miscellaneous shape properties
    properties['RDKit_Asphericity_cavity'] = Descriptors3D.Asphericity(probe)
    properties['RDKit_Eccentricity_cavity'] = Descriptors3D.Eccentricity(probe)
    properties['RDKit_InertialShapeFactor_cavity'] = Descriptors3D.InertialShapeFactor(probe)
    properties['RDKit_SpherocityIndex_cavity'] = Descriptors3D.SpherocityIndex(probe)
    properties['RDKit_NormalizedInertiaRatio1/3_cavity'] = Descriptors3D.NPR1(probe)

    return properties

############# ------- HANDLER FNXS
def get_all_properties(struc, probe, cmplx, id, prox_radius=3.5, alpha=0):
    properties = pd.DataFrame()
    props = {}
    try:
        props = RDKit_Properties(struc, probe, cmplx)
        print(pd.DataFrame(props, index=[0]))
        properties = pd.DataFrame(props, index=[0])
        #properties = pd.concat([properties, pd.DataFrame(props, index=[0])])
    except Exception as e:
        print('Error gathering RDKit properties ...')
        print(e)
    try:
        props = PyVista_Properties(struc, probe, id, prox_radius, alpha)
        print(pd.DataFrame(props, index=[0]))
        if properties.empty:
            properties = pd.DataFrame(props, index=[0])
        else:
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
    except Exception as e:
        print('Error gathering PyVista properties ...')
        print(e)
    try:
        props = Morfeus_Properties(struc, probe, id, prox_radius)
        print(pd.DataFrame(props, index=[0]))
        if properties.empty:
            properties = pd.DataFrame(props, index=[0])
        else:
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
    except Exception as e:
        print('Error gathering Morfeus properties ...')
        print(e)
    try:
        props = DBSTEP_Properties(struc, probe, id, prox_radius)
        print(pd.DataFrame(props, index=[0]))
        if properties.empty:
            properties = pd.DataFrame(props, index=[0])
        else:
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
    except Exception as e:
        print('Error gathering DBSTEP properties ...')
        print(e)

    return properties

def main(cav_file, cplx_file, mol_file, id, rdkit, dbstep, morfeus, pyvista, out='out'):
    properties = pd.DataFrame()
    struc = read_STRUCTURE_file(mol_file)
    probe = read_CAVITY_file(cav_file)
    cplx = read_CPLX_file(cplx_file)

    if not rdkit and not dbstep and not morfeus and not pyvista:
        print('No Descriptor method selected.')
        sys.exit()

    if rdkit:
        try:
            props = RDKit_Properties(struc, probe)
            print(pd.DataFrame(props, index=[0]))
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
        except Exception as e:
            print('Error gathering RDKit properties ...')
            print(e)

    if dbstep:
        try:
            props = DBSTEP_Properties(struc, probe, id, prox_radius)
            print(pd.DataFrame(props, index=[0]))
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
        except Exception as e:
            print('Error gathering DBSTEP properties ...')
            print(e)

    if morfeus:
        try:
            props = Morfeus_Properties(struc, probe, id, prox_radius)
            print(pd.DataFrame(props, index=[0]))
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
        except Exception as e:
            print('Error gathering Morfeus properties ...')
            print(e)

    if pyvista:
        try:
            props = PyVista_Properties(struc, probe, id, prox_radius, alpha)
            print(pd.DataFrame(props, index=[0]))
            properties = properties.merge(pd.DataFrame(props, index=[0]), how='inner', left_index=True, right_index=True)
        except Exception as e:
            print('Error gathering PyVista properties ...')
            print(e)

    export_descriptors(properties, out)

if __name__ == '__main__':
    ##! PARSE CMDLINE INPUT !##
    parser = argparse.ArgumentParser(prog='Spatial Molding for Approchable Rigid Targets (SMART)',description='Molecular Descriptor Calculations.',epilog='Uses multiple open-source python packages to compute SMART molecular descriptors.')
    parser.add_argument('-cav', required=True) # cavity conformer file
    parser.add_argument('-cplx', required=True) # complex conformer file
    parser.add_argument('-struc', required=True) # structure file
    parser.add_argument('-id', required=True) # metal/binding site ID
    parser.add_argument('-r', required=False, default=3.5) # proximal radius (pyvista, morfeus, dbstep)
    parser.add_argument('-a', required=False, default=0) # alpha (pyvista)
    parser.add_argument('-o', required=False, default='out') # output xlsx name
    parser.add_argument('--pyvista', required=False, action='store_true', default=False)
    parser.add_argument('--rdkit', required=False, action='store_true', default=False)
    parser.add_argument('--dbstep', required=False, action='store_true', default=False)
    parser.add_argument('--morfeus', required=False, action='store_true', default=False)

    args = parser.parse_args()

    starttime = time.time()
    main(args.cav, args.cplx, args.struc, args.id, args.rdkit, args.dbstep, args.morfeus, args.pyvista, args.o) #rdkit, dbstep, morfeus, pyvista, out
    print("--- %s seconds ---" % (time.time() - starttime))
