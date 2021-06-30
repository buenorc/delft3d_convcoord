# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 15:07:54 2021

@author: rafae
"""


import numpy as np



def cm(lon):
    
    lon = np.rad2deg(lon)   # coordenadas geodesicas (em radianos)
    fuso = round(((180-lon)/6)+1)


    MC   = -1*(183-6*fuso)
    
    
    return np.deg2rad(MC)


lat   = 41.97247        # coordenadas geodesicas (em decimal degrees)
lon   = 2.388670        # coordenadas geodesicas (em decimal degrees)


def geotoutm(lon,lat):

    lat = np.deg2rad(lat)   # coordenadas geodesicas (em radianos)
    lon = np.deg2rad(lon)   # coordenadas geodesicas (em radianos)

    lon_mc = cm(lon)             # meridiano central (em radianos)
    achat = 3.35281066e-03       # achatamento do elipsoide (WGS84)
    semi_eixo = 6.3781370e+06    # semi_eixo_maior do elipsoide (WGS84)
 
    if (lat > 0):
        offy = 0.

    else:
        offy = 10000000 

    k0 = 1. - (1./2500.)
    equad = 2.*achat - (achat**2)
    elinquad = equad/(1. - equad)

    aux1 = equad**2
    aux2 = aux1*equad
    aux3 = np.sin(2*lat)
    aux4 = np.sin(4*lat)
    aux5 = np.sin(6*lat)
    aux6 = (1. - equad/4. - 3.*aux1/64. - 5.*aux2/256.)*lat
    aux7 = (3.*equad/8. + 3.*aux1/32. + 45.*aux2/1024.)*aux3
    aux8 = (15.*aux1/256. + 45.*aux2/1024.)*aux4
    aux9 = (35.*aux2/3072.)*aux5

    n = semi_eixo/np.sqrt(1-equad*np.sin(lat)**2)
    t = np.tan(lat)**2
    c = elinquad*np.cos(lat)**2
    ag = (lon-lon_mc)*np.cos(lat)
    m = semi_eixo*(aux6 - aux7 + aux8 - aux9)

    aux10 = (1.-t+c)*ag*ag*ag/6.
    aux11 = (5.-18.*t+t*t+72.*c-58.*elinquad)*(ag**5)/120.
    aux12 = (5.-t+9.*c+4.*c*c)*ag*ag*ag*ag/24.
    aux13 = (61.-58.*t+t*t+600.*c-330.*elinquad)*(ag**6)/720.  

    x = 500000. + k0*n*(ag + aux10 + aux11)                       # UTM (m)
    y = offy + k0*(m + n*np.tan(lat)*(ag*ag/2. + aux12 + aux13))  # UTM (m)

    return x,y


input_path ='input_degree/'

bath = input_path + 'sau_reservoir_samples_mod'
land = input_path + 'sau-reservoir'
grid = input_path + 'grid80'


samples = np.loadtxt(bath+'.xyz')

x_utm, y_utm = [], []

for i in range(len(samples[:,0])):
    x, y = geotoutm(samples[i,0],samples[i,1])
    
    x_utm.append(x)
    y_utm.append(y)
    
np.savetxt(bath+'_utm.xyz',np.column_stack((x_utm,y_utm,samples[:,2])),fmt='%.10f')


with open(land+'.ldb', 'r') as f:
    
    i = 0
    
    line_x = f.readline()
    line_y = f.readline()
    line_z = f.readline()
    code   = f.readline()
    size   = f.readline()
    
    lon = []
    lat = []
    
    for line in f:
        
        line = line.strip()
        camp = line.split(' ') 
        
        x,y = geotoutm(float(camp[0]), float(camp[1]))
        
        lon.append(x)
        lat.append(y)

lon = np.array(lon)
lat = np.array(lat)

        
with open(land+'_utm.ldb', 'w') as f:
    
    f.write(line_x)
    f.write(line_y)
    f.write(line_z)
    f.write(code)
    f.write(size)
    
    for i in range(len(lon)):
        f.write(' '+str(lon[i])+' '+str(lat[i])+' 0\n')

        
        
        
with open(grid+'.grd', 'r') as f:
    
    i = 0
    
    system    = f.readline()
    grid_size = f.readline()
    dimension = f.readline()
    
    line = grid_size.strip()
    num  = line.split('\t')  
       
    values_idx, values_idy = [], []
    
    for line in f:
        
        line = line.strip()
        camp = line.split(' ')
        
        val = np.array(camp[2:],float)
        
        if(i < int(num[1])):
            values_idx.append(val)
        else:
            values_idy.append(val)
        
        i = i + 1

values_idx = np.array(values_idx)        
values_idy = np.array(values_idy)        
        
for i in range(len(values_idx[:,0])):
    for j in range(len(values_idx[0,:])):
        
        if values_idx[i,j]!=0 and values_idy[i,j]!=0:
            values_idx[i,j], values_idy[i,j] = geotoutm(values_idx[i,j], values_idy[i,j])
            
            
with open(grid+'_utm.grd', 'w') as f:
    
    f.write(system)
    f.write(grid_size)
    f.write(dimension)


    for ni in range(int(num[1])):
        f.write('Eta= {}'.format(ni+1))
        for mi in range(int(num[0])):
            f.write(' {}'.format(values_idx[ni,mi]))
        f.write('\n')
    
    for ni in range(int(num[1])):
        f.write('Eta= {}'.format(ni+1))
        for mi in range(int(num[0])):
            f.write(' {}'.format(values_idy[ni,mi]))
        f.write('\n')    

        
            
        
        
    
    
    


