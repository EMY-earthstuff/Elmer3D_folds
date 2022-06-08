Files needed to generate synthetic glaciers, apply surging behaviour, 
and extract fold structure are contained here. Note that the full 
workflow can be combined with paraview macros (can be saved as python
scripts) and bash scripts to automate the paraview particle tracking.

Files available:

-DEMs: folder containing the synthetic DEMs for the ref geometry.

-Bash scripts: buildup_run, junk_fold_run, and triple_fold_run will run
the initial setup + spin-up, discraded surge, and 3x surge simulations
sequentially. Need to update paths.

-Fortran functions: MB.f90, sliding_funk.f90, and sliding_funk_nsurges_q.f90
are used in the .sif files. See bash scripts for compiling them.

-SIF files: synthetic_SS.sif, synthetic_TR_buildup.sif, synthetic_TR_big.sif,
and synthetic_TR_nfoldtest_n5.sif are used to run the different simulations

-Python scripts: Simulate_tributaries_2.0.py creates the initial synthetic 
DEMs, Makegeo.py (obtained on Elmer/Ice courses webpage) creates the
input files for gmsh, time_fix_m.py updates velocity values in .vtu
files for use in paraview particle tracking, particle_advector_ver2_paraview.py 
opens the point cloud files exported from paraview (export with every 
time step = True) and calculates planar geometry, and 
read_structure_outputs_depthdiscriminate.py opens the outputs of the 
previous script and plots them.