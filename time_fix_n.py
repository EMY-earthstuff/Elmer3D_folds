# -*- coding: utf-8 -*-
###script to adjust velocity magnitudes to variable timestep duration for particle tracking
###created by EMY

import vtk
import meshio
import numpy as np

time_range = np.linspace(1,300, 300)

for t in time_range:
    
    if 40 < t < 90 or 130 < t < 180 or 220 < t < 270: 
    
        #some finangling for weird file naming
        if t < 10:
            pad0 = str(00)
        elif 10 <= t < 100:
            pad0 = str(0)
        else:
            pad0 = ""

        # The source file
        file_name = "synthetic_surge__n_fold_n5_t0" + pad0 + str(int(t)) + ".vtu"
        # The output file to be written
        file_out = "synthetic_surge__n_fold_n5_t0" + pad0 + str(int(t)) + ".vtu"


        #read the current timestep .vtu file
        mesh = meshio.read(file_name)
        #extract velocities and create modified array
        new_vels = mesh.point_data['velocity']/12.5
        #update velocities
        mesh.point_data['velocity'] = new_vels
        
        #write the modified vtu
        meshio.Mesh(
            mesh.points,
            mesh.cells,
            mesh.point_data,
            mesh.cell_data
            ).write(file_out)
    
    else:
        pass