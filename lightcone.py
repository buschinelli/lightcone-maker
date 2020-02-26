import astropy.units as u
from astropy.cosmology import Planck15 as p15
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo
import numpy as np


def make_lightcone(nu_i, nu_f, nu_spacing, N, spatial_spacing):
    nu_size = nu_f - nu_i
    lightcone = np.zeros((nu_size,N,N))
    nu_emitted = 420
    for layer in range(nu_size):
        nu = nu_spacing*(layer + nu_i)
        z_given_nu = nu_emitted/nu - 1
        com_dist = cosmo.comoving_distance(z_given_nu).value
        intZ = int(z_given_nu - 1)
        print(intZ)
        data, redshift = data_library(np.array([intZ,13]), N)
        red1 = redshift[0]
        red2 = redshift[1]
        com_dist_red1 = cosmo.comoving_distance(red1).value
        com_dist_red2 = cosmo.comoving_distance(red2).value
        loc_layer_data1 = int((com_dist_red1/com_dist)//spatial_spacing)
        loc_layer_data2 = int((com_dist_red2/com_dist)//spatial_spacing)
        lightcone[layer] = (data[1,loc_layer_data2] - data[0,loc_layer_data1])*((z_given_nu - red1)/(red2 - red1)) + data[1,loc_layer_data1]

    return lightcone

def data_library(redshift_range, N):
    name = 'delta_T_v3_z'
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
            print(name+name_z)
            if count == 0:
                data1 = np.fromfile(name+name_z,dtype=np.float32)
                red1 = round(float(z_run),4)
            elif count == 1:
                data2 = np.fromfile(name+name_z,dtype=np.float32)
                red2 = round(float(z_run),4)
            z_run += 0.5
            z_run = round(z_run,4)
            print('data imported for redshift ', z_run)
            count += 1
        except:
            z_run += 0.5
            z_run = round(z_run,4)
        z_run = str(z_run)
        name_z = str(z_run+zero_str)
        
    data1 = data1.reshape((N,N,N))
    data2 = data2.reshape((N,N,N))
    data = np.stack((data1,data2),axis=0)
    redshift[0] = red1
    redshift[1] = red2
    return data, redshift