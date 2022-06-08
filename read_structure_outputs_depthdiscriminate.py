#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###created by EMY
###open strike and dipe files output by particle_advector and plot


import numpy as np
from scipy.interpolate import griddata as gd
from scipy.interpolate import interp1d as i1d
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import mplstereonet

##Open strike/dip/age/ files and plot specific time/space components###########
    
sim_num = 288

file_prefix = ['dips', 'strikes', 'ages'] 
file_sufix = '.csv'

#horizontal (time of injection) selection bounds
injection_track = [1, 77]
injection_track2 = [78, 167]
injection_track3 = [168, 257]


#vertical (depth) selection bounds
depth_track = [0,21]
#depth intervals for plotting
depth_interv = [0, 5, 10, 15, 20, 23]

t0 = np.linspace(1, 20, 20)
t1 = np.linspace(21, 40, 20)
t2 = np.linspace(40.08, 44, 50)
t3 = np.linspace(45, 64, 20)
t4 = np.linspace(65, 84, 20)
t5 = np.linspace(84.08, 88, 50)
t6 = np.linspace(89, 108, 20)
t7 = np.linspace(109, 128, 20)
t8 = np.linspace(128.08, 132, 50)
t9 = np.linspace(133, 152, 20)
times_interv = np.concatenate((t0, t1, t2, t3, t4, t5, t6, t7, t8, t9))


for a in range(2, sim_num):
    dips = np.loadtxt(file_prefix[0] + str(a) + file_sufix)
    strikes = np.loadtxt(file_prefix[1] + str(a) + file_sufix)
    ages = np.loadtxt(file_prefix[2] + str(a) + file_sufix)
    
    fig = plt.figure(figsize=(8,8))
    
    for d in range(0, len(depth_interv)-1):
    
        dips_of_interest = dips[depth_interv[d]:depth_interv[d+1]][(ages[depth_interv[d]:depth_interv[d+1]] >= injection_track[0]) & (ages[depth_interv[d]:depth_interv[d+1]] <= injection_track[1])]
        strikes_of_interest = strikes[depth_interv[d]:depth_interv[d+1]][(ages[depth_interv[d]:depth_interv[d+1]] >= injection_track[0]) & (ages[depth_interv[d]:depth_interv[d+1]] <= injection_track[1])]
    
        #filter out nans for contour plots
        nanlocs = np.where(np.isnan(dips_of_interest) == 1)
        dips_unnaned = np.delete(dips_of_interest, nanlocs)
        strikes_unnaned = np.delete(strikes_of_interest, nanlocs)
        error_locs = np.where(dips_unnaned == 0.)
        dips_corrected = np.delete(dips_unnaned, error_locs)
        strikes_corrected = np.delete(strikes_unnaned, error_locs)
    

        
        #Plot structure for timestep
        ax = fig.add_subplot(111, projection='stereonet')
        labelstring = str(depth_interv[d]) + "-" + str(depth_interv[d+1])
        ax.pole(strikes_corrected, dips_corrected, 'o', markersize=5, label=labelstring)

        

    #ax.density_contourf(strikes_corrected, dips_corrected, measurement='poles', cmap='BuGn')
    ax.grid()
    ax.legend(title="particle distance from bed")
    timestep = "Timestep: " + str(times_interv[a])
    ax.text(150, -300, timestep)
    



    
    if a < 10:
        pre_a = "00"
    else:
        pre_a = "0" 
    
      
    plt.savefig('fold1_stereo' + pre_a + str(a) +'.png')
    plt.close()
    


































