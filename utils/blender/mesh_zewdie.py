#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from stl import mesh


'''This is a helper script to generate a zewdie isosurface given a set of parameters.
	the isosurface is written to zewdie.stl and is read in by import_dump.py when 
	vizualizing 1CPN in blender'''

import pdb

# define zewdie parameters
# these are parameters from first fit to 3spn data
ps0   = 45.3796835435872         
ps000 = 0.905835010046335        
pscc2 = -0.458169670802355       
ps220 = 1.74891489040291         
ps222 = -1.33065657833228        
ps224 = 1.03998788841889         

pe0   = 43.2151374522855         
pe000 = 0.128604190793286        
pecc2 = -0.431844707607592       
pe220 = 0.915662953592864        
pe222 = 1.72160108874421         
pe224 = 4.74995077023779         

# these are parameters from Kepper2008, Stehr2008
#ps0 = 1.0; ps000 = 1.6957; pscc2 = -0.7641; ps220 = -0.1480; ps222 = -0.2582; ps224 = 0.5112;
#pe0 = 1.0; pe000 = 1.9152; pecc2 =  2.7322; pe220 =  1.2633; pe222 =  2.3440; pe224 = 1.0101;

#Zewdie
#ps0 = 1.0; ps000 = 2.34  ; pscc2 = -1.52  ; ps220 = -0.64  ; ps222 = -0.69  ; ps224 = 1.97  ;
#pe0 = 1.0; pe000 = 2.24  ; pecc2 =  3.58  ; pe220 =  3.16  ; pe222 =  4.29  ; pe224 = 1.30  ;


def zewdie(rn,f0,f1,f2):
  # Business time
  S000 = 1.;
  S202 = (3*f1*f1 - 1)/(2.*np.sqrt(5));
  S022 = (3*f2*f2 - 1)/(2.*np.sqrt(5));
  S220 = (3*f0*f0 - 1)/(2.*np.sqrt(5));
  S222 = (2 - 3*f1*f1 - 3*f2*f2 - 3*f0*f0 + 9*f0*f1*f2)/np.sqrt(70);
  S224 = (1 + 2*f0*f0 - 5*f1*f1 - 5*f2*f2 - 20*f0*f1*f2 + 35*f1*f1*f2*f2)/(4.*np.sqrt(70));
  sig = ps0*(S000*ps000 + S220*ps220 + S222*ps222 + S224*ps224 + (S022 + S202)*pscc2);
  eps = pe0*(S000*pe000 + S220*pe220 + S222*pe222 + S224*pe224 + (S022 + S202)*pecc2);
  U = 4*eps*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));
  return U

def find_zero(phi1,phi2,theta):
  '''Return r that zeros the function at given phi'''
  # parameters
  r=ps0*5.0;
  dr = 1;
  factor = 0.1
  TOLERENCE = 1e-6

  #find zero
  f1 = np.cos(phi1)
  f2 = np.cos(phi2)
  f0 = np.cos(theta)
  eprev = -1;
  while True:
    e = zewdie(r,1.0,f1,f2)
    if (e > 0):
        if (e < TOLERENCE):
          break 
        else:
          r = r+dr       # go back a step
          dr = dr*factor # shrink step size
    r = r-dr
  return r

def main():
   nphi = 40
   ntheta = 20
   surf_xy  = np.zeros((nphi,2))
   n= (nphi-2)*(ntheta) + 2
   #n= ntheta*nphi
   surf_xyz = np.zeros((n,3))
   phis = np.linspace(0,np.pi,nphi)
   #thetas   = np.linspace(0,2*np.pi,ntheta,endpoint=False)
   thetas   = np.linspace(0,2*np.pi,ntheta) # this removed striping
   count = 0;
   for i in range(nphi):
      phi = phis[i]
      r = find_zero(phi,phi,0)
      r *= 0.5
      surf_xy[i,0] = r*np.cos(phi)
      surf_xy[i,1] = r*np.sin(phi)
      for j in range(ntheta):
        theta = thetas[j]
        surf_xyz[count,0] = r*np.cos(phi)
        surf_xyz[count,1] = r*np.sin(phi)*np.cos(theta)
        surf_xyz[count,2] = r*np.sin(phi)*np.sin(theta)
        count += 1
        if i == 0 or i == nphi-1:
            break

   ## now get facets
   facets = [] 
   for i in range(nphi-1):
      for j in range(ntheta):
        if i == 0: #triangles at pole
            if j != ntheta -1: 
              a = 0
              b = a+j+1
              c = b+1
            #else:
            #  a = 0
            #  b = a+j+1
            #  c = a+0+1
            facets.append([a,b,c])
        elif i==nphi-2: # triangles at pole
            if j != ntheta -1: 
              a = (i-1)*ntheta + j + 1 + 1
              b = a - 1
              c = n-1
            #else:
            #  a = (i-1)*ntheta + j + 1
            #  b = (i-1)*ntheta + 0 + 1
            #  c = n-1
            facets.append([a,b,c])
        elif j == ntheta - 1: # wrap in shape
            b = (i-1)*ntheta + j + 1
            a = (i-1)*ntheta + 0 + 1
            c = i*ntheta + j + 1
            d = i*ntheta + 0 + 1
            #facets.append([a,b,d]) #draw 2 triangles
            #facets.append([b,c,d])
        else: 
            a = (i-1)*ntheta + j + 1
            b = a + 1
            c = i*ntheta + j + 1 + 1
            d = c - 1
            facets.append([a,b,c]) #draw 2 triangles
            facets.append([a,c,d])

   #np.savetxt('surf_xy.dat',surf_xy)
   #np.savetxt('surf_xyz.dat',surf_xyz)

   # write stl mesh
   data = np.zeros(len(facets),dtype=mesh.Mesh.dtype)
   for i in range(len(facets)):
     a = surf_xyz[facets[i][0]]
     b = surf_xyz[facets[i][1]]
     c = surf_xyz[facets[i][2]]
     tmp = np.append(np.append(a,b),c)
     tmp = np.reshape(tmp,(3,3))
     data['vectors'][i] = tmp
   mymesh = mesh.Mesh(data)
   mymesh.save('zewdie.stl')

   #plot 
   #x,y,z = surf_xyz[:,0], surf_xyz[:,1], surf_xyz[:,2]
   #fig = plt.figure()
   #ax = fig.add_subplot(111,projection='3d')
   #ax.scatter(x,y,z)
   #ax.set_xlabel('x')
   #ax.set_ylabel('y')
   #plt.show()

if __name__ == "__main__":
    main()
  
