#!/usr/bin/env python
import numpy as np
import math as m
import sys
from molecule import *
from vect_quat_util import *

# current lh parameters come from Luque et al (2014)

class LinkerHistone(object):
    __slots__ = ('ctd_mass', 'gh_mass', 'num_in_gh', 'init_factor',
                'lequil', 'beta', 'salt_scale', 'lnucllh','salt_prefactor',
                'ctd_shape', 'ctd_charges', 'ctd_beads','ldyadlh',
                'ctd_bond_length','bonds','angles','dihedrals',
                'kgh','kctd','kbendgh','ktorsgh','ghupsi',
                'gh_data','ctd_charges','lnuclctd','ldyadctd')

    def __init__(self,d):
      self.bonds = []
      self.angles = []
      self.dihedrals = []
      #Data taken from Luque et al 2014 and translated to fit our reference system
      self.gh_data={'gh1': {'pos': np.array([52.228, 38.609,10.7562]), 'charge': {5:-4.19, 15:-4.12, 80:3.64, 150:-3.29} },
          'gh2': {'pos': np.array([48.939, 30.429, -6.935]), 'charge':{5:2.08, 15:2.13, 80:3.11, 150:4.22} },
          'gh3': {'pos': np.array([57.043, 32.130, 5.458]),  'charge':{5:6.90, 15:6.83, 80:7.47, 150:8.48} },
          'gh4': {'pos': np.array([45.038, 37.282, 3.817]),  'charge':{5:-0.28, 15:-0.18, 80:-0.09, 150:0.28} },
          'gh5': {'pos': np.array([41.897, 39.769, -8.150]), 'charge':{5:1.10, 15:1.38, 80:1.89,  150:2.08} },
          'gh6': {'pos': np.array([57.580, 41.316, -4.953]), 'charge':{5:1.21, 15:1.53, 80:2.53,  150:3.27} } }

      # the charges are taken from the sequence of the H1.4 tails
      self.ctd_charges = [0,2,2,3,0,4,0,2,2,4,0,2,3,2,2,2,2,2,2,2,2,3]; # array of H1.4 ctd CG charges
      # set the rest of the parameters
      self.kgh = 50.0     # gh-gh bond strength
      self.kctd = 0.1     # ctd-ctd bond strength
      self.kbendgh = 10.0 # angle strength
      self.ktorsgh = 10.0 # dihedral strength
      self.ghupsi = 180.0 # dihedral between gh and ctd to prevent rotation
      self.salt_scale = { 5:1.1063,
                          15:1.1644,
                          80:1.4815,
                          150:1.68000}; # salt scaling of ctd charges at 150 mM
      # From Luque2014 SI "the 1.5 is a prefactor that was included in the original model to reproduce mesoscopic chromatin structure"
      self.salt_prefactor = 1.0 # 1.5
      self.num_in_gh = 6
      self.ctd_beads = 22
      #self.init_factor = 1.69 # initialization factor for lhist
      self.init_factor = 1.0 # initialization factor for lhist
      self.lequil = 15.0 # bond equil length
      self.ldyadlh = 33.0 # dyad-lh length
      self.ldyadctd = 45.0 # ctd-dyad length
      self.lnucllh = self.ldyadlh + d # nucl-lh length
      self.lnuclctd = self.ldyadctd + d # nucl-ctd length
      self.beta = 110.0; # beta for the ctd
      self.ctd_bond_length = 1.0 # gh ctd equil length

    # function to calculate the bonded topology for the linker histone
    def calculate_topology(self):
      # bonds first
      self.bonds.append(calc_bonds(self.gh_data['gh1']['pos'],self.gh_data['gh3']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh1']['pos'],self.gh_data['gh4']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh2']['pos'],self.gh_data['gh3']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh2']['pos'],self.gh_data['gh6']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh4']['pos'],self.gh_data['gh5']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh5']['pos'],self.gh_data['gh6']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh1']['pos'],self.gh_data['gh6']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh2']['pos'],self.gh_data['gh5']['pos']))
      self.bonds.append(calc_bonds(self.gh_data['gh3']['pos'],self.gh_data['gh4']['pos']))

      # angles seconds
      self.angles.append(calc_angles(self.gh_data['gh3']['pos'],self.gh_data['gh1']['pos'],self.gh_data['gh4']['pos']))
      self.angles.append(calc_angles(self.gh_data['gh1']['pos'],self.gh_data['gh4']['pos'],self.gh_data['gh5']['pos']))
      self.angles.append(calc_angles(self.gh_data['gh3']['pos'],self.gh_data['gh2']['pos'],self.gh_data['gh6']['pos']))
      self.angles.append(calc_angles(self.gh_data['gh2']['pos'],self.gh_data['gh5']['pos'],self.gh_data['gh6']['pos']))
      self.angles.append(calc_angles(self.gh_data['gh1']['pos'],self.gh_data['gh6']['pos'],self.gh_data['gh5']['pos']))
      self.angles.append(calc_angles(self.gh_data['gh2']['pos'],self.gh_data['gh3']['pos'],self.gh_data['gh4']['pos']))

      # finally dihedrals
      self.dihedrals.append(calc_dihedrals(self.gh_data['gh3']['pos'],self.gh_data['gh1']['pos'],self.gh_data['gh4']['pos'],self.gh_data['gh5']['pos']))
      self.dihedrals.append(calc_dihedrals(self.gh_data['gh3']['pos'],self.gh_data['gh2']['pos'],self.gh_data['gh6']['pos'],self.gh_data['gh5']['pos']))
      return

def write_lhist_variables(fnme,lhist,salt):
  fnew = open (fnme,"w")
  fnew.write("#Error message to make sure lammps specified salt is valid\n")
  fnew.write("if \"${salt} != %f\" then \"print 'Error! Linker Histone charges are only valid at %.2fmM. Generate a new in.lammps if you want to simulate at a different salt!'\" quit\n\n" % (salt,salt))
  fnew.write( "# Defining linker histone variables for LAMMPS\n")
  fnew.write( "variable gha equal %f\n" % lhist.bonds[0])
  fnew.write( "variable ghb equal %f\n" % lhist.bonds[1])
  fnew.write( "variable ghc equal %f\n" % lhist.bonds[2])
  fnew.write( "variable ghd equal %f\n" % lhist.bonds[3])
  fnew.write( "variable ghe equal %f\n" % lhist.bonds[4])
  fnew.write( "variable ghf equal %f\n" % lhist.bonds[5])
  fnew.write( "variable ghg equal %f\n" % lhist.bonds[6])
  fnew.write( "variable ghh equal %f\n" % lhist.bonds[7])
  fnew.write( "variable ghi equal %f\n" % lhist.bonds[8])
  fnew.write( "variable kbondgh equal %f\n" % lhist.kgh)
  fnew.write( "variable kbondctd equal %f\n" % lhist.kctd)
  fnew.write( "variable kbendgh equal %f\n" % lhist.kbendgh)
  fnew.write( "variable ktorsgh equal %f\n" % lhist.ktorsgh)
  fnew.write( "variable ghalpha equal %f\n" % lhist.angles[0])
  fnew.write( "variable ghbeta equal %f\n" % lhist.angles[1])
  fnew.write( "variable ghgamma equal %f\n" % lhist.angles[2])
  fnew.write( "variable ghdelta equal %f\n" % lhist.angles[3])
  fnew.write( "variable ghzeta equal %f\n" % lhist.angles[4])
  fnew.write( "variable ghomega equal %f\n" % lhist.angles[5])
  fnew.write( "variable ghphi equal %f\n" % int(lhist.dihedrals[0]))
  fnew.write( "variable ghpsi equal %f\n" % int(lhist.dihedrals[1]))
  fnew.write( "variable ghupsi equal %f\n" % lhist.ghupsi)
  fnew.write( "variable beta equal %f\n" % lhist.beta)
  fnew.write( "variable ctda equal %f\n" % lhist.lequil)
  fnew.write( "variable ldgh equal %f\n" % lhist.ldyadlh)
  fnew.write( "variable lngh equal %f\n" % lhist.lnucllh)
  fnew.write( "variable ldctd equal %f\n" % lhist.ldyadctd)
  fnew.write( "variable lnctd equal %f\n" % lhist.lnuclctd)


# returns a bond from the positions
def calc_bonds(pos1,pos2):
    bond = np.linalg.norm(pos2-pos1)
    return bond

# returns an equilibrium angle from the positions
def calc_angles(pos1,pos2,pos3):
    v1 = pos1 - pos2
    v2 = pos3 - pos2
    cos_angle = np.dot(v1,v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.rad2deg(m.acos(cos_angle))
    return angle

# returns an equilibrium dihedral from the positions
def calc_dihedrals(pos1,pos2,pos3,pos4):
    v1 = -1.0*(pos2 - pos1)
    v2 = pos3 - pos2
    v2 /= np.linalg.norm(v2)
    v3 = pos4 - pos3
    v = v1 - np.dot(v1,v2)*v2
    w = v3 - np.dot(v3,v2)*v2

    x = np.dot(v,w)
    y = np.dot(np.cross(v2,v),w)
    psi = np.degrees(np.arctan2(y,x))
    # this keeps the dihedrals with consistent sign conventions as lammps
    # not my favorite implementation
    return -psi

#Modular function that adds in the linker histones if selected
def add_linker_histones(molecule,lhist,param):
  # mape site names to numbers
  typemap = {'nucl': 1,
             'bead': 2,
             'ghost': 3,
             'gh' : 4,
             'ctd' : 5}
  if param.salt not in [5,15,80,150]:
    print "Error! Salt must be 5,15,80,150 for linker histone. Salt=%0.2f is invalid!" % param.salt
    exit()

  # calculate the bonds
  lhist.calculate_topology()

  iellipsoid = len(molecule.ellipsoids)
  molcount = 2
  for ellipsoid in molecule.ellipsoids:
    # after the nucleosome is set we need to add in the linker histone
    if ellipsoid.mytype == 1:
        quat = ellipsoid.quat
        pos  = ellipsoid.pos
        fvu  = quat_fvu_rot(np.eye(3),quat)
        fvu0 = np.eye(3)
        molid = ellipsoid.molid
        molcount += 1

        # the globular head gets added first
        for igh in range(lhist.num_in_gh):
          # print gh bead to string
          typestr = 'gh%i' % (igh+1)
          lvec = lhist.gh_data[typestr]['pos']
          # positions need to be rotated about the axis to match 1CPN notation
          quat_v_rot = tu2rotquat(m.pi/2,fvu0[1])
          # the 54.0 degrees is to transfer from Schlick group nucleosome reference system to 1CPN fvu
          quat_f_rot = tu2rotquat(m.pi*(0.5+54.0/180),fvu0[0])
          quat_fv_rot = quat_multiply(quat_f_rot, quat_v_rot)
          lh_quat = quat_multiply(quat, quat_fv_rot)
          #lh_quat = quat_normalize(lh_quat)
          lh_pos = quat_vec_rot(lvec,lh_quat) # rotate linker vector by nucl quat
          lh_pos = np.add(lh_pos,pos) # translate linker position by nucl pos

          mytype = typemap['gh']
          molid = molcount
          molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,lh_pos/param.lengthscale, quat,lhist.gh_data[typestr]["charge"][param.salt],lhist.ctd_shape,molid))
          iellipsoid += 1

        # for the other AAs
        for ictd in range(lhist.ctd_beads):
          beta_init = lhist.beta*m.pi/180.0
          if ictd == 0:
            lh_pos = np.add(lh_pos,np.multiply(lhist.ctd_bond_length,fvu[2])) #eq dist is 10 angstroms
            q = tu2rotquat(beta_init,fvu[0]) #rotate by 110/1.69 degrees (for the initial structure) about f
            lh_quat = quat_multiply(q,quat)
            lh_fvu = quat_fvu_rot(fvu0,lh_quat)  # update linker fvu
          else:
            # update the rotation of the next lh bead so that the zigzag pattern appears
            if ictd % 2 == 0:
                q = tu2rotquat(beta_init, fvu[0])
            else:
                q = tu2rotquat(-beta_init, fvu[0])
            lh_quat = quat_multiply(q,quat)   # update quat
            lh_fvu = quat_fvu_rot(fvu0,lh_quat)  # update fvu
            lh_mag = lhist.lequil/lhist.init_factor

            # set position using f,v,u and pre-calculated lengths
            lh_pos = np.add(lh_pos,np.multiply(lh_mag, -lh_fvu[2]))

          mytype = typemap["ctd"]
          # the salt_scale is the DiSCO calculated scaling for 150mM
          charge = lhist.ctd_charges[ictd]*lhist.salt_scale[param.salt]*lhist.salt_prefactor
          molid = molcount
          molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,lh_pos/param.lengthscale,lh_quat,charge,lhist.ctd_shape,molid))
          iellipsoid += 1
