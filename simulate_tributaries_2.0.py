# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 16:39:23 2020
@author: EMY
"""
#create synthetic bed and sfc DEM. Output the DEM info needed to use with Elmer
#/Ice griddatainterpolator.


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import griddata
#from shapely.geometry import Polygon

###geometry parameters#########################################################
#slope
alphar = math.radians(1.75)
alphal = math.radians(1.75)
alpha = math.radians(1.75)
#surface angle
sfc_angle = 2.
#tributary angle
beta = math.radians(90.)
#trib_right height
hr = 1.
#trib_left height
hl = 1.
#trib_right width
#dr = 5
#trib_left width
#dl = 5.



#option for depth dependent ds and scaled widths for tributaries
dr_fac = 5.
dl_fac = 5.
dr = hr*dr_fac
dl = hl*dl_fac

#main trunk width, tied to width of smaller tributaries before scaling
dr_prime = dr*math.sin(beta)
dl_prime = dr*math.sin(beta)
D = 2*dl*math.sin(beta/2)

#main trunk height, tied to area of smaller tributaries to avoid mismatch
#get 'real' areas for scaled dr/dl
area_r = ((dr)*(hr)*np.pi)/2
area_l = ((dl)*(hl)*np.pi)/2
area_M = area_r + area_l
H = 2*area_M/(np.pi*D)




#smalltrib length
Tln = 20.
#Main trunk length
Mln = 40.
#Z scaling factor
depth_factor = 600.
#horizontal scaling factor*smooth
H_factor = 800.
#interpolation smoothing value
smooth = 10.
#factor used to separate main trunk from tributaries, increase for high angle 
#combined with high dr/dl widths to avoid interpolation artefacts
X_buffer = 0.
#interp method used
interp_method = 'linear'
#mask factor for proportion of channel lengths used in contour
Contour_trim = 3./4.
#mask factor for transect end locations
Transect_trim = 1.05
#cutoff for setting ice thickness to surface elevation
cut = 50.
#factor to exaggerate valley walls
valley_fac = 10.
#factor to exaggerate sfc slope: standard = 35, steep = 50
sfc_slope_fac = 35.


###Create Trib_r###############################################################
#subdivisions of pi, change to 2*np.pi for full cylinder
pis = np.linspace(0, -np.pi, 64)
#length coordinates
les = np.linspace(0, Tln, 20)
#form surface
Pis, Les = np.meshgrid(pis, les)
#convert to 'real' coordinates
xs_tr = dr * np.cos(Pis)
zs_tr = hr * np.sin(Pis)
ys_tr = Les
#get location-wise dh based on alpha to 'rotate' surfaces
alpha_surface = np.zeros(Les.shape)
for i in range(0, len(alpha_surface[:,])):
    dh = Les[i,0]*math.tan(alphar)
    alpha_surface[i] = alpha_surface[i] + dh
#add 'rotation to surface
zs_rot_tr = alpha_surface + zs_tr

#plot for debug
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xs_tr, ys_tr, zs_rot_tr, color='b')
#plt.title('Trib_r')
#plt.show()

###Create Trib_l###############################################################
#subdivisions of pi, change to 2*np.pi for full cylinder
pis = np.linspace(0, 2 * -np.pi/2, 64)
#length coordinates
les = np.linspace(0, Tln, 20)
#form surface
Pis, Les = np.meshgrid(pis, les)
#convert to 'real' coordinates
xs_tl = dl * np.cos(Pis)
zs_tl = hl * np.sin(Pis)
ys_tl = Les
#get location-wise dh based on alpha to 'rotate' surfaces
alpha_surface = np.zeros(Les.shape)
for i in range(0, len(alpha_surface[:,])):
    dh = Les[i,0]*math.tan(alphal)
    alpha_surface[i] = alpha_surface[i] + dh
#add 'rotation to surface
zs_rot_tl = alpha_surface + zs_tl

#plot for debug
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xs_tl, ys_tl, zs_rot_tl, color='b')
#plt.title('Trib_l')
#plt.show()

###Create Main T###############################################################
#subdivisions of pi, change to 2*np.pi for full cylinder
pis = np.linspace(0, 2 * -np.pi/2, 64)
#length coordinates
les = np.linspace(0, Mln, 20)
#form surface
Pis, Les = np.meshgrid(pis, les)
#convert to 'real' coordinates
xs_mt = D * np.cos(Pis)
zs_mt = H * np.sin(Pis)
ys_mt = Les
#get location-wise dh based on alpha to 'rotate' surfaces
alpha_surface = np.zeros(Les.shape)
for i in range(0, len(alpha_surface[:,])):
    dh = Les[i,0]*math.tan(alphal)
    alpha_surface[i] = alpha_surface[i] + dh
#get max Z and lower to 0 to coincide with tributaries
max_z = np.max(alpha_surface)
#add 'rotation to surface
zs_rot_mt = alpha_surface + zs_mt - max_z

#plot for debug
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xs_mt, ys_mt, zs_rot_mt, color='b')
#plt.title('Main')
#plt.show()

####make a box for them########################################################
#box dimensions XYZ: Y is lateral, X is downglacier, z is up/down
#get trib y span, and Y is 2*trib span + D, use complementary angles and such
#get trib x span, and X is trib span + Mln, use complementary angles and such
#***ONLY HOLDS TRUE FOR SYMMETRIC TRIBUTARIES
phi = (np.pi-beta)/2
trib_xs_up = Tln*math.cos(phi)
Phi = np.pi/2 - phi
trib_xs_dn = dr*math.cos(Phi)
trib_xs = trib_xs_up + trib_xs_dn
#get trib x span, and Y is 2*trib span + D, use complementary angles and such
#phi_p = np.arctan(dl/Tln)
#trib_xs = dl/math.sin(Phi)
phi_p = np.arcsin(trib_xs_up/Tln)
trib_ys = Tln*math.cos(phi_p)
Y = 2*trib_ys + D
X = trib_xs + Mln + X_buffer
Z = np.max(zs_rot_tl)-np.min(zs_rot_mt)
Xs = np.linspace(0, X, int(X*smooth))
Ys = np.linspace(0, Y, int(Y*smooth))
#Get X and Y surfaces
XX, YY = np.meshgrid(Xs, Ys)

####Locate tribs in box########################################################

#Get middle coordinates to start (x7, y7)
x0r = trib_xs_up
y0r = Y/2.
#get angle and interval of Tln change
dTln = Tln/len(xs_tr[:,0])
Tln_angle = np.pi/2 - beta/2
#calc magnitude of change componentwise based on grid spacing
dy_Tlnr = dTln*math.sin(Tln_angle)
dx_Tlnr = dTln*math.cos(Tln_angle)
#rinse wash repat for dd
dd = dr/len(xs_tr[0,:])
dd_angle = beta/2.
#again
dy_dr = dd*math.sin(dd_angle)
dx_dr = dd*math.cos(dd_angle)
#create empty arrays for new x and y coordinates
xs_r_new = np.zeros(xs_tr.shape)
ys_r_new = np.zeros(ys_tr.shape)
#always + for y t_right based on chosen starting coord (x7,y7), -/+ for x dtln/dd
for i in range(0, len(xs_r_new[:,0])):
    xi = x0r - (i*dx_Tlnr)
    yi = y0r + (i*dy_Tlnr)
    for j in range(0, len(xs_r_new[0,:])):
        xs_r_new[i,j] = xi + (j*dx_dr)
        ys_r_new[i,j] = yi + (j*dy_dr)
        
###Trib_l
#Get middle coordinates to start (x7, y7): shared with other tribuatry
x0l = trib_xs_up
y0l = Y/2.
#get angle and interval of Tln change
dTln = Tln/len(xs_tl[:,0])
Tln_angle = np.pi/2 - beta/2
#calc magnitude of change componentwise based on grid spacing
dy_Tlnl = dTln*math.sin(Tln_angle)
dx_Tlnl = dTln*math.cos(Tln_angle)
#see trib r
dd = dl/len(xs_tl[0,:])
dd_angle = beta/2.
#see trib l
dy_dl = dd*math.sin(dd_angle)
dx_dl = dd*math.cos(dd_angle)
#create empty arrays for new x and y coordinates
xs_l_new = np.zeros(xs_tl.shape)
ys_l_new = np.zeros(ys_tl.shape)
#always - for y t_left based on chosen starting coord (x7,y7), -/+ for x dtln/dd
#dd direction are reversed!!!! because of array starting configuration.
for i in range(0, len(xs_l_new[:,0])):
    xi = x0l - (i*dx_Tlnl)
    yi = y0l - (i*dy_Tlnl)
    for j in range(0, len(xs_l_new[0,:])):
        xs_l_new[i,j] = xi + (j*dx_dl)
        ys_l_new[i,j] = yi - (j*dy_dl)
        
###Main T
y0m = Y - trib_ys
x0m = X
dMln = Mln/len(xs_mt[:,0])
dD = D/len(xs_mt[0,:])
xs_mt_new = np.zeros(xs_mt.shape)
ys_mt_new = np.zeros(ys_mt.shape)
for i in range(0, len(xs_mt_new[:,0])):
    xi = x0m - (i*dMln)
    yi = y0m
    for j in range(0, len(xs_mt_new[0,:])):
        xs_mt_new[i,j] = xi
        ys_mt_new[i,j] = yi - (j*dD)

        

###merge glaciers into a single mesh###########################################
X_msh = np.concatenate((xs_r_new, xs_l_new, xs_mt_new))
Y_msh = np.concatenate((ys_r_new, ys_l_new, ys_mt_new))
Z_msh = np.concatenate((zs_rot_tr, zs_rot_tl, zs_rot_mt))

interp_input = np.zeros((len(X_msh.flatten()),2))
interp_input[:,0] = X_msh.flatten()
interp_input[:,1] = Y_msh.flatten()
surf_interp_funk_griddata = griddata(interp_input, Z_msh.flatten(), (XX, YY), method=interp_method)

#make bed sfc positive
Z_min = np.min(Z_msh)
Zs = (surf_interp_funk_griddata + np.abs(Z_min)) * depth_factor
Z_max = np.nanmax(Zs)

###get the corners for contour file and mask###################################
#pretty crude boundaries implementation, could be made dynamic

#start point is confluence middle
x7 = (trib_xs_up - 1.) 
y7 = (int(Y)/2.)

x1 = (x7 + (dx_dl*len(xs_l_new[0,:])) * Transect_trim) 
y1 = (y7 -(dy_dl*len(xs_l_new[0,:])) * Transect_trim)

x4 = (x7 + (dx_dr*len(xs_r_new[0,:])) * Transect_trim)
y4 = (y7 + (dy_dr*len(xs_r_new[0,:])) * Transect_trim)

x8 = (x7 - (dx_Tlnl*Tln) * Contour_trim)
y8 = (y7 - (dy_Tlnl*Tln) * Contour_trim)

x6 = (x7 - (dx_Tlnr*Tln) * Contour_trim)
y6 = (y7 + (dy_Tlnr*Tln) * Contour_trim)

x0 = (x7 + (dx_dl*len(xs_l_new[0,:]) * Transect_trim) - (dx_Tlnl*Tln) * Contour_trim)
y0 = (y7 - (dy_dl*len(ys_l_new[0,:]) * Transect_trim) - (dy_Tlnl*Tln) * Contour_trim)

x5 = (x7 - (dx_Tlnr*Tln * Contour_trim) + (dx_dr*len(xs_r_new[0,:])) * Transect_trim)
y5 = (y7 + (dy_Tlnr*Tln * Contour_trim) + (dy_dr*len(xs_r_new[0,:])) * Transect_trim)

x2 = X * Contour_trim
y2 = trib_ys * 1./Transect_trim
x3 = X * Contour_trim
y3 = (trib_ys + D) * Transect_trim

x9 = x0
y9 = y0

#makelist
Xs = [x0, x1, x2, x3, x4, x5, x6, x7, x8, x9]
Ys = [y0, y1, y2, y3, y4, y5, y6, y7, y8, y9]

###Create surface array and make mask for bed##################################
surface = np.copy(XX)
for i in range(0, len(surface[0,:])):
    if XX[0,i] > X - Mln:
        angle = alpha
    else:
        angle = alphar
    surface[:,i] = (Z_max + valley_fac) - (surface[0,i]*sfc_slope_fac)/math.cos(angle)
    
    #surface[:,i] = 

    
#mask bed DEM with surface
H_ice = surface - Zs 
sfc_id = np.where(H_ice < cut)
Zs[sfc_id] = surface[sfc_id]
    



###export contour and surf/bed synthetic DEM files#############################
#save the contour
contour = np.asarray((Xs,Ys))*H_factor
np.savetxt('contour.dat', contour.transpose())

#save the DEMs by flattening the arrays
DEM_bed_out = np.zeros((len(XX.flatten()),3))
DEM_bed_out[:, 0] = XX.flatten()*H_factor
DEM_bed_out[:, 1] = YY.flatten()*H_factor
DEM_bed_out[:, 2] = Zs.flatten()
#replace nans with surface elevation 
DEM_bed_out[np.where(np.isnan(DEM_bed_out) == 1)] = -9999.
        
np.savetxt('bed_DEM.dat', DEM_bed_out)

DEM_out = np.zeros((len(XX.flatten()),3))
DEM_out[:, 0] = XX.flatten()*H_factor
DEM_out[:, 1] = YY.flatten()*H_factor
DEM_out[:, 2] = surface.flatten()
#replace nans with surface elevation 
DEM_out[np.where(np.isnan(DEM_out) == 1)] = -9999.
        
np.savetxt('sfc_DEM.dat', DEM_out)

DEM_empty_out = DEM_bed_out.copy()
DEM_empty_out[:,2][np.where((DEM_bed_out[:,2] ==-9999.) == 0.)] + 50.
np.savetxt('sfc_DEM_empty.dat', DEM_empty_out)


#various plots for monitoring/debug

fig = plt.figure()
plt.imshow((surface - Zs).transpose())
plt.colorbar()
plt.scatter(np.asarray(Ys)*smooth, np.asarray(Xs)*smooth)
plt.title('H')
plt.xlabel("x" + str(H_factor/smooth) + 'm')
plt.ylabel("x" + str(H_factor/smooth) + 'm')

fig = plt.figure()
plt.imshow(Zs.transpose())
plt.colorbar()
plt.scatter(np.asarray(Ys)*smooth, np.asarray(Xs)*smooth)
plt.title('bed')
plt.xlabel("x" + str(H_factor/smooth) + 'm')
plt.ylabel("x" + str(H_factor/smooth) + 'm')


fig = plt.figure()
plt.imshow(surface.transpose())
plt.colorbar()
plt.scatter(np.asarray(Ys)*smooth, np.asarray(Xs)*smooth)
plt.title('sfc')
plt.xlabel("x" + str(H_factor/smooth) + 'm')
plt.ylabel("x" + str(H_factor/smooth) + 'm')


fig = plt.figure()
plt.imshow(DEM_empty_out.transpose())
plt.colorbar()
plt.scatter(np.asarray(Ys)*smooth, np.asarray(Xs)*smooth)
plt.title('sfc = 0')
plt.xlabel("x" + str(H_factor/smooth) + 'm')
plt.ylabel("x" + str(H_factor/smooth) + 'm')

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xs_r_new, ys_r_new, zs_rot_tr, color='b')
#ax.plot_surface(xs_l_new, ys_l_new, zs_rot_tl, color='b')
#ax.plot_surface(xs_mt_new, ys_mt_new, zs_rot_mt, color='b')
#ax.plot_surface(X_msh, Y_msh, Z_msh, color='b')
#plt.xlim([0, 60])
#plt.ylim([0, 40])
#plt.title('all')
#plt.show()

#fig = plt.figure()
#ax1 = fig.add_subplot(111, projection='3d')
#ax1.plot_surface(XX, YY, Zs, rstride=3, cstride=3, linewidth=1,
#antialiased=True,cmap='inferno', alpha = 0.5)
#ax1.plot_surface(XX, YY, surface, rstride=3, cstride=3, linewidth=1,
#antialiased=True,cmap='inferno', alpha = 1)


####print DEM info#############################################################
#these are the DEM infos needed to use the grid data interpolator script in .sif file
print('DEM size is:')
print(XX.shape)
print('X length and Y length are:')
print(str(np.max(XX.flatten()*H_factor)) + '   ' + str(np.max(YY.flatten()*H_factor)))