#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###Created by EMY
###script to extract plane attitudes from point cloud and output for every timestep

import numpy as np
from scipy.interpolate import griddata as gd
from scipy.interpolate import interp1d as i1d
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import mplstereonet

#find nearest function: https://stackoverflow.com/questions/2566412/find-
#nearest-value-in-numpy-array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


#function for exporting .csv file with (z,l,strike,dip) data as x(t)
def measure_structured_planes(X_s, Y_s, Z_s, ts):
    
    #define normals for vertical y plane(0-360 great circle) and equator
    N_y = np.asarray((1, 0, 0))
    N_hor = np.asarray((0, 0, 1))
    
    #create out array of shape z*l(t) and len(id, z, l(t), strike, dip)
    out_array = np.zeros(((len(Z_s)-1)*(ts-1),4))
    kk = 0
    #create function return array of shape Z_s
    strikes = np.zeros(Z_s.shape)
    dips = np.zeros(Z_s.shape)
    dip_dirs = np.zeros(Z_s.shape)
    #idents = np.zeros(Z_s.shape).tolist()
    #idents_truths = np.zeros(Z_s.shape)
    
    #number of planes = n-1 z * m-1 l(t)
    #iterate through z-1
    for ii in range(0, len(Z_s)-1):
        #iterate through l(t)-1
        for jj in range(0, ts-1):
            #print(X_s[ii,jj])
          
            #get p1...p4
            P1 = (X_s[ii,jj], Y_s[ii,jj], Z_s[ii,jj])
            P2 = (X_s[ii,jj+1], Y_s[ii,jj+1], Z_s[ii,jj+1])
            P3 = (X_s[ii+1,jj], Y_s[ii+1,jj], Z_s[ii+1,jj])
            P4 = (X_s[ii+1,jj+1], Y_s[ii+1,jj+1], Z_s[ii+1,jj+1])
            
            #get vectors defining plane
            V12 = np.asarray(P4) - np.asarray(P2)
            V13 = np.asarray(P3) - np.asarray(P2)
            #Xproduct to obtain normal
            N_surf = np.cross(V12, V13)
            #N_surf = np.array([0.5,-1,1])
            #surf = N_surf * P1
            #Get normal to vertical strike plane
            N_surf_strike = N_surf*(1.,1.,0.)
            N_surf_dip = np.cross(V12, V13)
            
            
            #Get angle between reference planes and plane of iterest
            alpha_y = np.arccos((np.dot(N_surf_strike, N_y)/(np.linalg.norm(N_surf_strike)*np.linalg.norm(N_y))))
            alpha_z = np.arccos((np.dot(N_surf_dip, N_hor)/(np.linalg.norm(N_surf_dip)*np.linalg.norm(N_hor))))
            
            #convert to degrees
            alpha_y_deg = (alpha_y*360.)/(2*np.pi)
            alpha_z_deg = (alpha_z*360.)/(2*np.pi)
            
            #use plane sign to homogenize conjugate normals to same dip direction surfaces
            #dd NE quadrant
            if np.sum(np.sign(N_surf) == [-1,-1,-1]) == 3 or np.sum(np.sign(N_surf) == [1,1,1]) == 3:
                if alpha_y_deg > 90:
                    alpha_y_deg = alpha_y_deg - 90.   
                else:
                    alpha_y_deg = 90. - alpha_y_deg  
                #break
            #dd NW quadrant
            elif np.sum(np.sign(N_surf) == [1,-1,-1]) == 3 or np.sum(np.sign(N_surf) == [-1,1,1]) == 3:
                if alpha_y_deg > 90:
                    alpha_y_deg = (270. - alpha_y_deg) +180.       
                else:
                    alpha_y_deg = 270. + alpha_y_deg   
                #break
            #dd SW quadrant
            elif np.sum(np.sign(N_surf) == [1,1,-1]) == 3 or np.sum(np.sign(N_surf) == [-1,-1,1]) == 3:
                if alpha_y_deg > 90:
                    alpha_y_deg = 90. - alpha_y_deg       
                else:
                    alpha_y_deg = (90. - alpha_y_deg) + 180.   
                #break
            #dd SE quadrant
            elif np.sum(np.sign(N_surf) == [1,-1,1]) == 3 or np.sum(np.sign(N_surf) == [-1,1,-1]) == 3:
                if alpha_y_deg > 90:
                    alpha_y_deg = 270. - alpha_y_deg      
                else:
                    alpha_y_deg = alpha_y_deg  + 90.
                #break
        
            ##add cardinal points rather than trust sketchy indexing based on patterns
            #due North
            elif  np.sum(np.sign(N_surf) == [0,-1,-1]) == 3 or np.sum(np.sign(N_surf) == [0,1,1]) == 3:
                alpha_y_deg = 0.
                #break
            #due South
            elif  np.sum(np.sign(N_surf) == [0,1,-1]) == 3 or np.sum(np.sign(N_surf) == [0,-1,1]) == 3:
                alpha_y_deg = 180.
                #break
            #due East
            elif  np.sum(np.sign(N_surf) == [1,0,1]) == 3 or np.sum(np.sign(N_surf) == [-1,0,-1]) == 3:
                alpha_y_deg = 90.
                #break
            #due East
            elif  np.sum(np.sign(N_surf) == [-1,0,1]) == 3 or np.sum(np.sign(N_surf) == [1,0,-1]) == 3:
                alpha_y_deg = 270.
                #break
            else:
                alpha_y_deg = np.nan
                     

            #append values to output
            out_array[kk, 0] = ii; out_array[kk, 1] = jj; out_array[kk, 2] = alpha_y_deg; out_array[kk, 3] = alpha_z_deg;  
            
            #convert normal azimuth to strike
            if alpha_y_deg < 90:
                alpha_y_deg = 360. - (90. - alpha_y_deg)
            else:
                alpha_y_deg = alpha_y_deg - 90.
            
            strikes[ii,jj] = alpha_y_deg
            dips[ii,jj] = alpha_z_deg
            dip_dirs[ii,jj] = alpha_y_deg + 90. 
            #idents[ii][jj] = np.sign(N_surf)
            
            kk += 1
            
            #print(alpha_y_deg, ii, jj)
    
    

    
    
    
        
    #correct for RHR using plane identity signs as index for quadrant
    dip_correct_locs = np.where(dips > 90)
    dips[dip_correct_locs] = 180. - dips[dip_correct_locs]   


    
    #save outputs                    
    #np.savetxt('planes' + str(ts) + '.csv', out_array)
    #create ages array 
    ages_array = np.linspace(0, ts, ts+1)
    birthday = np.zeros(dips.shape)
    for aa in range(0, len(birthday)):
        birthday[aa] = ages_array
    
    np.savetxt('dips' + str(ts) + '.csv', dips)
    np.savetxt('strikes' + str(ts) + '.csv', strikes)
    np.savetxt('ages' + str(ts) + '.csv', birthday)
    
    #fig = plt.figure()
    #plt.plot(strikes.flatten(), 'x')
    #plt.savefig('striketest' + pre_a + str(a) +'.png')
    #plt.close()


       
    return strikes, dips, dip_dirs


###load particle data##########################################################

sim_num = 287

file_prefix = 'particles.'
file_sufix = '.csv'

template = file_prefix + str(sim_num-1) + file_sufix
vertical_levels = np.unique(np.loadtxt(template, delimiter =',', skiprows = 1)[:,38])

#empty list to save attitudes for plotting from memory
Strike_list = []
Dip_list = []
dip_dir_list = []

for a in range(0, sim_num):
    file_in = file_prefix + str(a) + file_sufix
    timestep = np.loadtxt(file_in, delimiter =',', skiprows = 1)
    
    pages = np.absolute(timestep[:,41] - a)
    pids = timestep[:,38]
    pxs = timestep[:,45]
    pys = timestep[:,46]
    pzs = timestep[:,47]
    
    vertical_addresses = np.unique(pids)
    X_s = np.zeros((len(vertical_levels), a+1))*np.nan
    Y_s = np.zeros((len(vertical_levels), a+1))*np.nan
    Z_s = np.zeros((len(vertical_levels), a+1))*np.nan
    
    for idx in vertical_levels:
        data_locs = np.where(pids == idx)
        address = np.where(vertical_addresses == idx)
        iid = int(idx)
        
        if len(address[0]) == 0:
                pass
        else:
            X_s[iid][pages[data_locs].astype(int)] = pxs[data_locs]
            Y_s[iid][pages[data_locs].astype(int)] = pys[data_locs]
            Z_s[iid][pages[data_locs].astype(int)] = pzs[data_locs]



    #naming adjustments for generating .gif
    if a < 10:
        pre_a = "00"
    else:
        pre_a = "0"    


    if a < 2:
        pass
    else:
        strikes, dips, dip_dirs = measure_structured_planes(X_s, Y_s, Z_s, a)
        Strike_list.append(strikes)
        Dip_list.append(dips)
        dip_dir_list.append(dip_dirs)
        
        #Plot stereoplots
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='stereonet')
        #ax.line(90. - dips.flatten(), dip_dirs.flatten(), 'ko', markersize=5)
        
        #filter out nans for contour plots
        nanlocs_dip = np.where(np.isnan(dips.flatten()) == 1)
        nanlocs_dir = np.where(np.isnan(dip_dirs.flatten()) == 1)
        nanlocs_strike = np.where(np.isnan(strikes.flatten()) == 1)
        dips_unnaned = np.delete(dips.flatten(), nanlocs_dip)
        dirs_unnaned = np.delete(dip_dirs.flatten(), nanlocs_dir)
        strikes_unnaned = np.delete(strikes.flatten(), nanlocs_dip)
        error_locs = np.where(dips_unnaned == 0.)
        dips_corrected = np.delete(dips_unnaned, error_locs)
        dirs_corrected = np.delete(dirs_unnaned, error_locs)
        strikes_corrected = np.delete(strikes_unnaned, error_locs)
        
        #actual plotting
        #ax.line(90. - dips_corrected, dirs_corrected, 'ko', markersize=5)
        ax.pole(strikes_corrected, dips_corrected, 'ko', markersize=5)

        #ax.density_contourf(strikes_corrected, dips_corrected, measurement='poles', cmap='Reds')
        ax.grid()
        
        plt.savefig('stereo' + pre_a + str(a) +'.png')
        plt.close()
    
    #collapse into vectors for plotting
    Xpts = X_s[:,1:-1].flatten()
    Ypts = Y_s[:,1:-1].flatten()
    Zpts = Z_s[:,1:-1].flatten()
    
    # Plot X,Y,Z
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Xpts, Ypts, Zpts, c='red')

    plt.savefig('step' + pre_a + str(a) +'.png')
    plt.close()
    
    

###synthetic planes to test quadrant attribution and angle results        
# # ###TEST arrays for structure function
#  synx = np.zeros((2,2))
#  syny = np.zeros((2,2))
#  synz = np.zeros((2,2))

# # #270/45 N
#  synx = np.array(([0, 100], [0, 100]))
#  syny = np.array(([0, 0], [100, 100]))
#  synz = np.array(([50, 50], [0, 0]))
# # #090/45 S
#  synx = np.array(([0, 100], [0, 100]))
#  syny = np.array(([0, 0], [-100, -100]))
#  synz = np.array(([100, 100], [0, 0]))
# # #180/45
#  synx = np.array(([0, 0], [-100, -100]))
#  syny = np.array(([0, 100], [0, 100]))
#  synz = np.array(([100, 100], [0, 0]))
#  #000/45
#  synx = np.array(([0, 0], [100, 100]))
#  syny = np.array(([0, 100], [0, 100]))
#  synz = np.array(([50, 50], [0, 0]))
 
 
#   #225/45 dd NW
#  synx = np.array(([0, 10], [-10, 0]))
#  syny = np.array(([0, 10], [0, 10]))
#  synz = np.array(([10, 10], [0, 0]))
#    #315/45 dd NE
#  synx = np.array(([0, -10], [10, 0]))
#  syny = np.array(([0, 10], [10, 20]))
#  synz = np.array(([10, 10], [0, 0]))
#     #135/45 dd SW
#  synx = np.array(([0, -10], [-10, -20]))
#  syny = np.array(([0, 10], [-10, 0]))
#  synz = np.array(([10, 10], [0, 0]))
#      #045/45 dd SE
#  synx = np.array(([0, -10], [0, 10]))
#  syny = np.array(([0, -10], [-20, -10]))
#  synz = np.array(([10, 10], [0, 0]))

#  X_s = synx
#  Z_s = synz
#  Y_s = syny

        

            
            
            
            
            