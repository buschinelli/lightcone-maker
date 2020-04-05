import astropy.units as u
from astropy.cosmology import Planck15 as p15
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
import matplotlib.pyplot as plt


######## Stiching boxes function #########

##### Pseudo code 

def make_lightcone(nu_size, N, spatial_spacing):
    ## nu_size: integer - size of the frequency axis of lightcone
    ## N: integer - size of the boxes used to produce the lightcone
    ## spatial_spacing: 1D 3 entries array - spatial spacing of the boxes
    #### Initial and final frequency - choose them based on the box you have
    nu_i = 200 # with units
    nu_f = 100 # with units
    nu_spacing = -(nu_f-nu_i)/nu_size
    lightcone = np.zeros((N,N,nu_size))
    nu_emitted = 1420
    for layer in range(nu_size+1):
        nu = nu_spacing*(layer) + nu_f
        z_given_nu = nu_emitted/nu - 1
        zi = nu_emitted/(nu_i) - 1
        com_dist_z = cosmo.comoving_distance(z_given_nu).value
        com_dist_zi = cosmo.comoving_distance(zi).value
        intZ = int(z_given_nu - 0.5)
        #print(intZ)
        data, redshift = data_library(np.array([intZ,13]), N)
        red1 = redshift[0]
        red2 = redshift[1]
        loc_layer = int((com_dist_z - com_dist_zi)//spatial_spacing)
        if loc_layer >= N:
            while loc_layer >= N:
                loc_layer -= N
        #print(red1,red2)
        lightcone[:,:,layer-1] = (data[1,loc_layer] - data[0,loc_layer])*((z_given_nu - red1)/(red2 - red1)) + data[1,loc_layer]
        
    return lightcone

def data_library(redshift_range, N):
    ## redshift_range: 1D 2 entries array - range of redshift you want this code to find the boxes you want to use
    ## N: integer - size of the boxes
    name = 'delta_T_v3_z' #### Standard name to look for the boxes: change the boxes names in your library to this + the redshift they are.
                          #### Example: if the box is at redshift 5.5, its name must be delta_T_v3_z5.5 in your library
    z_run = redshift_range[0]
    zero_str = '0'
    z0 = redshift_range[0]
    z1 = redshift_range[1]
    Range = z1 - z0
    stop = False
    count = 0
    times = 0
    data1 = np.zeros((N,N,N))
    deta2 = np.zeros((N,N,N))
    name_z = str(z0)
    redshift = np.zeros((2))
    while stop == False:
        times += 1
        z_run = float(z_run)
        if count == 2 or times >= 100:
            stop = True
        try:
            #print(name+name_z)
            if count == 0:
                data1 = np.fromfile(name+name_z,dtype=np.float32)
                red1 = z_run
            elif count == 1:
                data2 = np.fromfile(name+name_z,dtype=np.float32)
                red2 = z_run
            z_run += 0.5
            z_run = round(z_run,4)
            #print('data imported for redshift ', z_run)
            count += 1
        except:
            z_run += 0.5
            z_run = round(z_run,4)
        z_run = str(z_run)
        name_z = str(z_run+zero_str)
    redshift[0] = red1
    redshift[1] = red2
    data1 = data1.reshape((N,N,N))
    data2 = data2.reshape((N,N,N))
    data = np.stack((data1,data2),axis=0)
    return data, redshift