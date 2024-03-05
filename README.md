# Spatial Molding for Approachable Rigid Targets (SMART)
<img width="195" alt="image" src="https://github.com/SigmanGroup/SMART-molecular-descriptors/assets/84196711/6da48469-0a27-4084-bf85-a8d474757127"> <img width="200" alt="image" src="https://github.com/SigmanGroup/SMART-molecular-descriptors/assets/84196711/84905748-9780-4411-a9cb-8d97dcc63e2a">

An open-source Python package for generation of SMART probe ensembles and calculation of SMART molecular descriptors.

## Dependencies
- pandas
- scipy
- mathutils
- RDKit
- morfeus (optional)
- dbstep (optional)
- pyvista (optional)
  
# Instructions
### Step 1: Select and Add Probe to Structure(s) of Interest
#### Option 1: Tail-to-Tip Reference Vector
~Insert img here~
```
import SMART_probe_utils as SMART
file = 'R-TCPTTL_1'
structure = SMART.ReadFile(file) # initialize structure
structure.reference_vector(tip=0, tail=1, dist=2.0) # define reference binding vector
probe = SMART.Probe('S_SiF2_12_cyclic.mol2') # initialize selected probe (in Probes/ folder)
docked = SMART.add_probe(structure, probe) # start docking protocol
SMART.ExportStructure(docked, file+'_probe') # (optional) export docked structure as .mol file
```
Where:

tip = 0-indexed ID of binding atom

tail = 0-indexed ID of reference atom

dist = distance of probe binding point from binding atom (Å)

Add probe using command line interface:
```
python3 SMART_probe_utils.py -f R-TCPTTL_1.mol2 -o R-TCPTTL_1_probe -id 0 -ref 1 -p S_SiF2_12_cyclic.mol2 -dist 2.0
```
#### Option 2: Angle Cross Product Reference Vector
Insert img here
```
import SMART_probe_utils as SMART
file = 'R-TCPTTL_1'
structure = SMART.ReadFile(file) # initialize structure
structure.reference_angle(tip=0, tails=[50,21], dist=2.0) # define reference binding vector
probe = SMART.Probe('S_SiF2_12_cyclic.mol2') # initialize selected probe (in Probes/ folder)
docked = SMART.add_probe(structure, probe) # start docking protocol
SMART.ExportStructure(docked, file+'_probe') # (optional) export docked structure as .mol file
```
Add probe using command line interface:
```
python3 SMART_probe_utils.py -f R-TCPTTL_1.mol2 -o R-TCPTTL_1_probe -id 0 -ref 1 -p S_SiF2_12_cyclic.mol2 -dist 2.0
```
#### Option 3: Define Geometry of Binding Center
Insert img here

#### Option 4: Auto-Detect Geometry of Binding Center
Insert img here

### Step 2: Generate Probe Conformational Ensemble
#### Option 1: Template Search
Insert img here

generate conformers by fitting and rotating a probe conformer ensemble template to cavity of interest
```
import SMART_conf_search as SMART_conf
SMART_conf.PARAMS.read_parameters({'NSTEP':20, 'MINROTATION':30, 'MAXROTATION':330}) # initialize custom search parameters
try:
    cavity, structure, cmplx = SMART_conf.CUSTOM_TEMPLATE_SEARCH(docked) # run template algorithm
except Exception as e:
    print('simulation failed for - ', file)
    print(e)
SMART_conf.save_out('SMART_confs_'+file) # (optional) export ensemble as mol2 file
```
Run simulation using command line interface:
```
python3 SMART_conf_search.py -f R-TCPTTL_1_probe.mol -method template -p PARAMS.txt
```
#### Option 2: Torsional Search
Insert img here

*this search method does not support cyclic probes*
```
import SMART_conf_search as SMART_conf
SMART_conf.PARAMS.read_parameters({'NSTEP':500, 'MINROTATION':30, 'MAXROTATION':330}) # initialize custom search parameters
try:
    cavity, structure, cmplx = SMART_conf.TORSIONAL_STEP_SEARCH(docked) # run torsional algorithm
except Exception as e:
    print('simulation failed for - ', file)
    print(e)
SMART_conf.save_out('SMART_confs_'+file) # (optional) export ensemble as mol2 file
```
Run simulation using command line interface:
```
python3 SMART_conf_search.py -f R-TCPTTL_1_probe.mol -method torsional -p PARAMS.txt
```
### Step 3: Calculate Molecular Descriptors
All descriptors can be gathered and compiled using the function ```get_all_properties()```
```
import SMART_descriptors as SMART_des
binding_ID = 0 # structure tip atom
properties = SMART_des.get_all_properties(cavity, structure, binding_ID, prox_radius=3.5, alpha=0)
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
