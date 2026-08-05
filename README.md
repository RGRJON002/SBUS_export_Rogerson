# SBUS_export_Rogerson
Series of  matlab scripts used to sample Eulerian tracer fields from CROCO along Lagrangian particle trajectories from ROFF. Forms part of the work for the paper: 'Carbon transport across the Cape Point Jet in the Southern Benguela Upwelling System'.

The MATLAB scripts call the CROCO processing tools:

(1) CROCO_TOOLS: refer to Penven et al. (2008) or download from https://www.croco-ocean.org/download/

List of Matlab Scripts:

**Eulerian filed extraction**

NPZD_get_flux_slice_V2.m <br>
NPZD_get_tracer_slice.m <br>

**Lagrangian sampling**
TRACK_PARTICLE.m <br>
AGE_TRACER.m <br>
FU3D.m <br>
FUtime.m <br>
Fto3D.m <br>
age_statistics.m <br>
mask_by_age.m <br>
sample_croco.m <br>

**Usage**
These scripts are fairly generic, and should work with other CROCO/ROMS outputs. Top extract a desired tracer field: edit and run NPZD_get_tracer_slice.m or if you want a flux across a depth interface use NPZD_get_flux_slice_V2.m. Both will create a series of .mat files which will be used by the later scripts.

With the desired Eulerian tracer fields extracted, they are horizontally sampled along the backwards Lagrangian trajectories produced by ROFF. The main script is TRACK_PARTICLE.m. It leverages all of the functions to sample the Eulerian fields and produces a data structure which is saved to a .mat file. This .mat contains all the data required for plotting and further analysis.   
