# Spatial Molding for Approachable Rigid Targets (SMART)
<img width="195" alt="image" src="https://github.com/SigmanGroup/SMART-molecular-descriptors/assets/84196711/6da48469-0a27-4084-bf85-a8d474757127"> <img width="200" alt="image" src="https://github.com/SigmanGroup/SMART-molecular-descriptors/assets/84196711/84905748-9780-4411-a9cb-8d97dcc63e2a">

An open-source Python package for generation of SMART probe ensembles and calculation of SMART molecular descriptors.

## Dependencies
- argparse
- pandas
- scipy
- mathutils
- RDKit
- dbstep (optional)
- pyvista (optional)
  
# Instructions
### Step 1: Select and Add Probe to Structure(s) of Interest
```
import SMART_probe_utils as SMART
file = 'R-TCPTTL_1'
ext = 'mol'
structure = SMART.ReadFile(ext, file, tip=0, tail=1) # initialize structure
probe = SMART.Probe('acyclic_6SiH2') # initialize probe
docked = SMART.addProbe(structure, probe, dist=2.0) # dock probe to structure
SMART.ExportStructure(docked, file+'_probe') # (optional) export docked structure as mol file
```
Where:

Tip = 0-indexed ID of binding atom

Tail = 0-indexed ID of reference atom

Add probe using command line interface:
```
python3 SMART_probe_utils.py -f R-TCPTTL_1.mol -o R-TCPTTL_1_probe -tip 0 -tail 1 -p acyclic_6SiH2 -dist 2.0
```
### Step 2: Generate Probe Conformational Ensemble
```
import SMART_conf_search as SMART_conf
PAR = {'FIXAT':[-22, 100, 1], 'NSTEP':500, 'MAXROTATION':330,'MINROTATION':30} #freeze atoms 22-100, 1
SMART_conf.PARAMS.read_parameters(PAR) # set search parameters
try:
    confs = SMART_conf.start(docked) # run algorithm
except RecursionError:
    print('simulation failed for', file)

SMART_conf.save_out('SMART_confs_'+file) # (optional) export ensemble as mol2 file
```
Run simulation using command line interface:
```
python3 SMART_conf_search.py ... TBD
```
### Step 3: Calculate Molecular Descriptors
```
import SMART_descriptors as SMART_des
binding_ID = 0 # structure tip atom
RDKit_props = SMART_des.RDKit_Properties(confs, binding_ID)
DBSTEP_props = SMART_des.DBSTEP_Properties(confs, binding_ID)
PyVista_props = SMART_des.PyVista_Properties(confs, binding_ID)
```
Calculate descriptors using command line interface:
```
python3 SMART_descriptors.py ... TBD
```
# Citations
This package:
- pending...

Literature Using SMART:
- Cammarota, R. C., Liu, W., Bacsa, J., Davies, H. M. L., & Sigman, M. S. Mechanistically Guided Workflow for Relating Complex Reactive Site Topologies to Catalyst Performance in C-H Functionalization Reactions. Journal of the American Chemical Society, 2022, 144(4), 1881–1898. [https://doi.org/10.1021/jacs.1c12198](https://doi.org/10.1021/jacs.1c12198)
- Lucas W. Souza, Beck R. Miller, Ryan C. Cammarota, Anna Lo, Ixchel Lopez, Yuan-Shin Shiue, Benjamin D. Bergstrom, Sarah N. Dishman, James C. Fettinger, Matthew S. Sigman, and Jared T. Shaw, ACS Catalysis, 2024, 14 (1), 104-115, [4256](https://doi.org/10.1021/acscatal.3c04256)https://doi.org/10.1021/acscatal.3c04256
