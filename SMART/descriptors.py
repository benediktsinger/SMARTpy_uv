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
    ''' read exported cavity file into CAV object
    '''

    name = cavity_file.split('.')[0]
    ext = cavity_file.split('.')[1]
    cav = None
    if ext == 'mol2':
        cav = Chem.MolFromMol2File(cavity_file, removeHs=False)
    elif ext == 'mol':
        cav = Chem.MolFromMolFile(cavity_file, removeHs=False)
    elif ext == 'xyz':
        cav = Chem.MolFromXYZFile(cavity_file, removeHs=False)
    elif ext == 'pdb':
        cav = Chem.MolFromPDBFile(cavity_file, removeHs=False)
    elif ext == 'sdf':
        suppl = Chem.SDMolSupplier(cavity_file, removeHs = False)
        for m in suppl:
            if not cav:
                cav = m
            else:
                cav.AddConformer(m.GetConformer(), assignId=True)
    return cav

def read_STRUCTURE_file(mol_file):
    ''' read exported structure file into STRUCT object
    '''

    name = mol_file.split('.')[0]
    ext = mol_file.split('.')[1]
    mol = None
    if ext == 'mol2':
        mol = Chem.MolFromMol2File(mol_file, removeHs=False, sanitize=False)
    elif ext == 'mol':
        mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
    elif ext == 'xyz':
        mol = Chem.MolFromXYZFile(mol_file, removeHs=False, sanitize=False)
    elif ext == 'pdb':
        mol = Chem.MolFromPDBFile(mol_file, removeHs=False, sanitize=False)
    elif ext == 'sdf':
        suppl = Chem.SDMolSupplier(mol_file, removeHs = False, sanitize=False)
        for m in suppl:
            if not mol:
                mol = m
            else:
                mol.AddConformer(m.GetConformer(), assignId=True)
    return mol

def ConfToMol(mol, conf):
    ''' convert CONF object to MOL object
    '''

    mol_ = Chem.Mol(mol)
    mol_.RemoveAllConformers()
    mol_.AddConformer(Chem.Conformer(conf))
    return mol_

def _coords_from_mols(confs):
    ''' convert CONF object to xyz coordinate, element name lists
    '''

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
    ''' convert MOL object to xyz coordinate, element name lists
    '''

    coords = []
    elems = []
    for i, atom in enumerate(struc.GetAtoms()):
        #atom = struc.GetAtomWithIdx(i)
        pos = struc.GetConformer().GetAtomPosition(i)
        coords.append(np.array(pos))
        elem = atom.GetSymbol()
        elems.append(elem)
    return coords, elems

def _max_radius(coords, TOL=2):
    ''' compute max radius across list of coords
        (increased by TOL factor, default = 2 Å)
    '''

    max_radius = 0
    for i in range(len(coords)):
        dist = math.dist(coords[i], coords[-1])
        if dist > max_radius:
            max_radius = dist
    return max_radius

def export_descriptors(properties, out, type='csv'):
    ''' Expprt SMART descriptors to .csv or .xlsx
    '''

    if type == 'xlsx':
        df = pd.DataFrame(properties)
        df.to_excel(out+'.xlsx')
    elif type == 'csv':
        df = pd.DataFrame(properties)
        df.to_csv(out+'.csv')

############# ------- MAIN FNXS
def DBSTEP_Properties(struc, probe, id, prox_radius=3.5, vol=True, sterimol=False):
    ''' Uses dbstep package to compute properties by:
            buried volume (dbstep_VBur_),
            sterimol (dbstep_Sterimol_)
    '''

    print('\nComputing DBSTEP descriptors')
    properties = {}
    coords, elems = _coords_from_mols(probe)

    coords.append(np.array(struc.GetConformer().GetAtomPosition(int(id))))
    elems.append(struc.GetAtomWithIdx(int(id)).GetSymbol())

    max_radius = _max_radius(coords, TOL=5)

    sterimol = db.dbstep(probe, atom1=len(coords)-1, atom2=0, sterimol=True)
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

def Morfeus_QO_Properties(struc, probe, id=0, z_axis=[], xz_plane=[], prox_radius=3.5, octant=True, quadrant=False):
    ''' Compute quadrant and octant analysis using the Morfeus package
        Parameters:
        - struct = STRUCT object
        - probe = PROBE object
        - id = atom ID for proximal sphere center (default: probe tether [0])
        - z_axis = list of atom IDs in z-axis (0-indexed)
        - xz_plane = list of atom IDs in xz-plane (0-indexed)
        - prox_radius = proximal radius (default: 3.5 Å)
        Returns:
        - Quadrant volume (O00_, O01_, O10_, O11_)
        - Octant volume (Q(+/-)00_, Q(+/-)01_, Q(+/-)10_, Q(+/-)11_)
    '''

    print('\nComputing Morfeus -OCTANT/QUADRANT- descriptors\n')
    properties = {}
    pcoords, pelems = _coords_from_mols(probe)
    scoords, selems = _coords_from_mols(struc)

    max_radius = _max_radius(pcoords, TOL=2)

    bv = BuriedVolume(pelems+selems, pcoords+scoords, 0, z_axis_atoms=[i+len(pcoords) for i in z_axis], xz_plane_atoms=[i+len(pcoords) for i in xz_plane], radius=prox_radius)#.compute_distal_volume(method="buried_volume")
    bv_max = BuriedVolume(pelems+selems, pcoords+scoords, 0, z_axis_atoms=z_axis, xz_plane_atoms=xz_plane, radius=max_radius)#.compute_distal_volume(method="buried_volume")

    if octant:
        bv.octant_analysis()
        bv_max.octant_analysis()

        octs = bv.octants
        octs_max = bv_max.octants

        properties['morfeus_O(+,+,+)_VOL_cavity'] = octs_max['buried_volume'][0]
        properties['morfeus_O(+,+,+)_%VOL_cavity'] = octs_max['percent_buried_volume'][0]
        properties['morfeus_O(-,+,+)_VOL_cavity'] = octs_max['buried_volume'][1]
        properties['morfeus_O(-,+,+)_%VOL_cavity'] = octs_max['percent_buried_volume'][1]
        properties['morfeus_O(-,-,+)_VOL_cavity'] = octs_max['buried_volume'][2]
        properties['morfeus_O(-,-,+)_%VOL_cavity'] = octs_max['percent_buried_volume'][2]
        properties['morfeus_O(+,-,+)_VOL_cavity'] = octs_max['buried_volume'][3]
        properties['morfeus_O(+,-,+)_%VOL_cavity'] = octs_max['percent_buried_volume'][3]
        properties['morfeus_O(+,+,-)_VOL_cavity'] = octs_max['buried_volume'][4]
        properties['morfeus_O(+,+,-)_%VOL_cavity'] = octs_max['percent_buried_volume'][4]
        properties['morfeus_O(-,+,-)_VOL_cavity'] = octs_max['buried_volume'][5]
        properties['morfeus_O(-,+,-)_%VOL_cavity'] = octs_max['percent_buried_volume'][5]
        properties['morfeus_O(-,-,-)_VOL_cavity'] = octs_max['buried_volume'][6]
        properties['morfeus_O(-,-,-)_%VOL_cavity'] = octs_max['percent_buried_volume'][6]
        properties['morfeus_O(+,-,-)_VOL_cavity'] = octs_max['buried_volume'][7]
        properties['morfeus_O(+,-,+)_%VOL_cavity'] = octs_max['percent_buried_volume'][7]

        properties['morfeus_O(+,+,+)_proxVOL_cavity'] = octs['buried_volume'][0]
        properties['morfeus_O(+,+,+)_%proxVOL_cavity'] = octs['percent_buried_volume'][0]
        properties['morfeus_O(-,+,+)_proxVOL_cavity'] = octs['buried_volume'][1]
        properties['morfeus_O(-,+,+)_%proxVOL_cavity'] = octs['percent_buried_volume'][1]
        properties['morfeus_O(-,-,+)_proxVOL_cavity'] = octs['buried_volume'][2]
        properties['morfeus_O(-,-,+)_%proxVOL_cavity'] = octs['percent_buried_volume'][2]
        properties['morfeus_O(+,-,+)_proxVOL_cavity'] = octs['buried_volume'][3]
        properties['morfeus_O(+,-,+)_%proxVOL_cavity'] = octs['percent_buried_volume'][3]
        properties['morfeus_O(+,+,-)_proxVOL_cavity'] = octs['buried_volume'][4]
        properties['morfeus_O(+,+,-)_%proxVOL_cavity'] = octs['percent_buried_volume'][4]
        properties['morfeus_O(-,+,-)_proxVOL_cavity'] = octs['buried_volume'][5]
        properties['morfeus_O(-,+,-)_%proxVOL_cavity'] = octs['percent_buried_volume'][5]
        properties['morfeus_O(-,-,-)_proxVOL_cavity'] = octs['buried_volume'][6]
        properties['morfeus_O(-,-,-)_%proxVOL_cavity'] = octs['percent_buried_volume'][6]
        properties['morfeus_O(+,-,-)_proxVOL_cavity'] = octs['buried_volume'][7]
        properties['morfeus_O(+,-,-)_%proxVOL_cavity'] = octs['percent_buried_volume'][7]

    if quadrant:
        bv.octant_analysis()
        bv_max.octant_analysis()
        quads = bv.quadrants
        quads_max = bv_max.quadrants

        properties['morfeus_Q(+,+)_VOL_cavity'] = quads_max['buried_volume'][1]
        properties['morfeus_Q(+,+)_VOL%_cavity'] = quads_max['percent_buried_volume'][1]
        properties['morfeus_Q(-,+)_VOL_cavity'] = quads_max['buried_volume'][2]
        properties['morfeus_Q(-,+)_VOL%_cavity'] = quads_max['percent_buried_volume'][2]
        properties['morfeus_Q(-,-)_VOL_cavity'] = quads_max['buried_volume'][3]
        properties['morfeus_Q(-,-)_VOL%_cavity'] = quads_max['percent_buried_volume'][3]
        properties['morfeus_Q(+,-)_VOL_cavity'] = quads_max['buried_volume'][4]
        properties['morfeus_Q(+,-)_VOL%_cavity'] = quads_max['percent_buried_volume'][4]

        properties['morfeus_Q(+,+)_proxVOL_cavity'] = quads['buried_volume'][1]
        properties['morfeus_Q(+,+)_proxVOL%_cavity'] = quads['percent_buried_volume'][1]
        properties['morfeus_Q(-,+)_proxVOL_cavity'] = quads['buried_volume'][2]
        properties['morfeus_Q(-,+)_proxVOL%_cavity'] = quads['percent_buried_volume'][2]
        properties['morfeus_Q(-,-)_proxVOL_cavity'] = quads['buried_volume'][3]
        properties['morfeus_Q(-,-)_proxVOL%_cavity'] = quads['percent_buried_volume'][3]
        properties['morfeus_Q(+,-)_proxVOL_cavity'] = quads['buried_volume'][4]
        properties['morfeus_Q(+,-)_proxVOL%_cavity'] = quads['percent_buried_volume'][4]

    return properties

def Morfeus_Properties(struc, probe, id=0, prox_radius=3.5, vol=True, sasa=False, sterimol=False):
    ''' Uses morfeus package to compute SMART descriptors
        Parameters:
        - struct = STRUCT object
        - probe = PROBE object
        - id = atom ID for proximal sphere center (default: probe tether [0])
        - prox_radius =  ()
        Returns:
        - BuriedVolume (V, prox_V),
        - Sterimol (B1, B5, L),
        - SASA (SASA, C_SASA, E_SASA)
    '''

    print('\nComputing Morfeus descriptors\n')
    properties = {}
    pcoords, pelems = _coords_from_mols(probe)
    scoords, selems = _coords_from_mol(struc)
    ccoords, celems = pcoords+scoords, pelems+selems

    max_radius = _max_radius(pcoords, TOL=5)

    if vol:
        bv = BuriedVolume(pelems, pcoords, 0, radius=prox_radius).compute_distal_volume(method="buried_volume")
        bv_max = BuriedVolume(pelems, pcoords, 0, radius=max_radius).compute_distal_volume(method="buried_volume")

        # VBur cavity properties
        properties['morfeus_VBur_VOL_cavity'] = bv_max.buried_volume
        properties['morfeus_VBur_%VOL_cavity'] = bv_max.percent_buried_volume
        properties['morfeus_VBur_proxVOL_cavity'] = bv.buried_volume
        properties['morfeus_VBur_%proxVOL_cavity'] = bv.percent_buried_volume
        properties['morfeus_VBur_distVOL_cavity'] = bv.distal_volume

    if sterimol:
        # Sterimol cavity properties
        sterimol = Sterimol(celems, ccoords, id, 0)
        properties['morfeus_Sterimol_L_cavity'] = sterimol.L_value
        properties['morfeus_Sterimol_B1_catvity'] = sterimol.B_1_value
        properties['morfeus_Sterimol_B5_cavity'] = sterimol.B_5_value

    if sasa:
        # SASA cavity properties
        probe_sasa = SASA(pelems, pcoords)
        struc_sasa = SASA(selems, scoords)
        cmplx_sasa = SASA(pelems+selems, pcoords+scoords)
        properties['morfeus_SASA_AREA_cavity'] = probe_sasa.area
        properties['morfeus_SASA_CSA_cavity'] = ( (probe_sasa.area + struc_sasa.area) - cmplx_sasa.area) / 2
        properties['morfeus_SASA_ESA_cavity'] = probe_sasa.area - properties['morfeus_SASA_CSA_cavity']
        properties['morfeus_SASA_%ESA_cavity'] = 100 * ( properties['morfeus_SASA_ESA_cavity'] / probe_sasa.area )

    return properties

def PyVista_Properties(struc, probe, id, prox_radius=3.5, a=0):
    ''' Uses pyvista package to assemble geometric objects and extract properties '''

    print('\nComputing PyVista descriptors\n')
    properties = {}
    pcoords, pelems = _coords_from_mols(probe)

    # gen point cloud from 3D coords
    p_ptcloud = pv.PolyData(pcoords)
    p_ptcloud.cast_to_pointset()

    # make shape object from pt cloud (delaunay triangulation)
    cavity = p_ptcloud.delaunay_3d(alpha=a)
    cavity = cavity.extract_geometry()
    cavity = cavity.extract_surface().triangulate()

    print('cavity', cavity.volume, cavity.area)

    # proximal/distal clipping
    sphere = pv.Sphere(radius=prox_radius, center=pcoords[0])
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

    print('proximal', prox.volume, prox.area)

    # full cavity properties
    properties['pyVista_delaunay_VOL_cavity'] =  float(cavity.volume)
    properties['pyVista_delaunay_AREA_cavity'] = float(cavity.area)
    # prox/dist cavity properties
    properties['pyVista_delaunay_proxVOL_cavity'] =  prox.volume
    properties['pyVista_delaunay_proxAREA_cavity'] = prox.area
    properties['pyVista_delaunay_distVOL_cavity'] =  cavity.volume - prox.volume

    return properties

def RDKit_Properties(struc, probe):
    ''' Uses rdkit package to compute properties
        WARNING: Beta version, not guaranteed to work!
    '''

    print('\nComputing RDKit descriptors\n')
    properties = {}
    probe_radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in probe.GetAtoms()]
    struc_radii = [Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in struc.GetAtoms()]
    cmplx_sradii = probe_radii + struc_radii
    cmplx = None
    for conf in probe.GetConformers():
        mol_ = ConfToMol(probe, conf)
        if not cmplx:
            cmplx = Chem.CombineMols(mol_, struc)
        else:
            cmplx = Chem.CombineMols(mol_, cmplx)

    area_struc = rdFreeSASA.CalcSASA(struc, struc_radii)
    area_cmplx = rdFreeSASA.CalcSASA(cmplx, cmplx_sradii)

    # SASA properties
    properties['RDKit_SASA_AREA_cavity'] = rdFreeSASA.CalcSASA(probe, probe_radii)
    #print('computed sasa')
    properties['RDKit_SASA_CSA_cavity'] = ( ( area_struc + properties['RDKit_SASA_AREA_cavity'] ) - area_cmplx ) / 2
    properties['RDKit_SASA_ESA_cavity'] = properties['RDKit_SASA_AREA_cavity'] - properties['RDKit_SASA_CSA_cavity']
    properties['RDKit_SASA_%ESA_cavity'] = 100 * ( properties['RDKit_SASA_ESA_cavity'] / properties['RDKit_SASA_AREA_cavity'] )
    # Miscellaneous shape properties
    #properties['RDKit_Mol_VOL_cavity'] = AllChem.ComputeMolVolume(probe)
    properties['RDKit_Asphericity_cavity'] = Descriptors3D.Asphericity(probe)
    properties['RDKit_Eccentricity_cavity'] = Descriptors3D.Eccentricity(probe)
    properties['RDKit_InertialShapeFactor_cavity'] = Descriptors3D.InertialShapeFactor(probe)
    properties['RDKit_SpherocityIndex_cavity'] = Descriptors3D.SpherocityIndex(probe)
    properties['RDKit_NormalizedInertiaRatio1/3_cavity'] = Descriptors3D.NPR1(probe)

    return properties

############# ------- HANDLER FNXS
def get_all_properties(struc, probe, id, prox_radius=3.5, alpha=0):
    ''' Get all main properties
        Calls functions:
           Morfeus_Properties()
           PyVista_Properties()
           DBSTEP_Properties()
         Returns:
           DataFrame
    '''

    properties = pd.DataFrame()
    props = {}
    #try:
    ##    props = RDKit_Properties(struc, probe)
    #    print(pd.DataFrame(props, index=[0]))
    #    properties = pd.DataFrame(props, index=[0])
        #properties = pd.concat([properties, pd.DataFrame(props, index=[0])])
    #except Exception as e:
    #    print('Error gathering RDKit properties ...')
    #    print(e)
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

def main(cav_file, mol_file, id, rdkit, dbstep, morfeus, pyvista, out='out'):
    properties = pd.DataFrame()
    struc = read_STRUCTURE_file(mol_file)
    probe = read_CAVITY_file(cav_file)
    #cplx = read_CPLX_file(cplx_file)

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
