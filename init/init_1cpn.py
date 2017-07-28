#!/usr/bin/env python

# Subversion Info
#$Revision$
#$LastChangedDate$

import sys
import math as m
import pdb
import numpy as np
import argparse
import sys
import os.path

from molecule import *
from vect_quat_util import *

# mape site names to numbers
typemap = {'nucl': 1,
           'bead': 2,
           'ghost': 3,
           'gh' : 4,
           'ctd' : 5}

def write_args(args,fnme):
   f=open(fnme,"w")
   command = ' '.join(sys.argv)
   f.write('# Inputs issued to command line\n')
   f.write("%s\n" %command)
   f.write('# All Parameters used\n')
   f.write("%s\n" %args)
   f.write("SVN $Revision$\n")
   f.write("SVN $LastChangedDate$\n")
   f.close()

def write_lammps_variables(fnme,param,geom):
  '''Using model lengths, determine the corect spring constants and write to lammps variables'''
  l0 = param.rise_per_bead
  fnew = open (fnme,"w")
  fnew.write( "# Defining model variables for LAMMPS\n")
  fnew.write( "variable la equal %f\n" % geom.a)
  fnew.write( "variable lb equal %f\n" % l0 )
  fnew.write( "variable lc equal %f\n" % geom.c)
  fnew.write( "variable ld equal %f\n" % (geom.d))
  fnew.write( "variable le equal %f\n" % (geom.e))
  #fnew.write( "variable li equal %f\n" % geom.i)
  #fnew.write( "variable lf equal %f\n" % geom.f)
  #fnew.write( "variable alpha equal %f\n" % (geom.alpha0*180.0/np.pi))
  fnew.write( "variable omega0 equal %f\n" % (param.twist_per_bead*180.0/np.pi))
  fnew.write( "variable mu equal %f\n" % (geom.mu*180.0/np.pi))
  fnew.write( "variable zeta equal %f\n" % (geom.zeta*180.0/np.pi))
  fnew.write( "variable eta equal %f\n" % (geom.eta*180.0/np.pi))
  fnew.write( "variable theta equal %f\n" % (geom.theta*180.0/np.pi))
  #fnew.write( "variable li equal %f\n" % (i))

  #      stehr      kJ->J   kJ->kcal  ->mol
  kbond= 1.10e-18 / 1e3     /4.18     * 6.022e23 / (l0/10.0)**3
  fnew.write( "variable kbond equal %f\n" % kbond)

  #       stehr      kJ->J   kJ->kcal  ->mol
  kangle= 2.06e-19 / 1e3     / 4.18    * 6.022e23 / (l0 / 10.0)
  fnew.write( "variable kangle equal %f\n" % kangle)

  #       stehr      kJ->J   kJ->kcal  ->mol
  ktwist= 2.67e-19 / 1e3     / 4.18    * 6.022e23 / (l0 / 10.0)
  fnew.write( "variable ktwist equal %f\n" % ktwist)
  fnew.write( "variable kalign equal %f\n" % (ktwist*10))

  #              stehr      kJ->J   kJ->kcal  ->mol
  #ktwistnucldna= 2.67e-19 / 1e3     / 4.18    * 6.022e23 / ( param.lnucldna / 10.0)
  #fnew.write( "variable ktwistnucldna equal %f\n" % ktwistnucldna)

  #           stehr      kJ->J   kJ->kcal  ->mol
  kbond_nucl= 1.10e-18 / 1e3     / 4.18    * 6.022e23 / (geom.c / 10.0)
  fnew.write( "variable kbondnucl equal %f\n" % (10*kbond_nucl))

  fnew.close()

def gen_nucl_nucl_bonds(molecule,offset):
   # nucl and nucl + offset are bonded together

   #add aditional bond type
   molecule.nbond_type = 4

   # find nucl indicies
   natom = len(molecule.ellipsoids)
   nucleosomes = []
   for i in range(natom):
       if (molecule.ellipsoids[i].mytype == 1):
          nucleosomes.append(i+1)

   nnucl = len(nucleosomes)

   #bond them
   for i in range(nnucl-offset):
        bondtype = 4
        atom1 = nucleosomes[i]
        atom2 = nucleosomes[i+offset]
        molecule.bonds.append(Bond(bondtype,atom1,atom2))

def set_bonded_interactions(molecule,lhbool):
  # =================
  # key of types
  # =================
  # atom_type 1  nucleosome
  # atom_type 2  bead
  # atom_type 3  stem (ghost site)
  # atom_type 4  globular head of linker histone
  # atom_type 5  c-terminal domain of linker histone
  if lhbool:
      molecule.natom_type = 5
  else:
      molecule.natom_type = 3

  # bond_type 1  bead-nucleosome, length c
  # bond_type 2  bead-bead
  # bond_type 3  bead-bead (enter and exit nucl), length e
  # bond_type 4  bead-centerghost, length j
  # bond_type 5  nucl-centerghost, length i
  # bond_type 6-14 globular head interactions
  # bond_type 15  ctd-ctd
  # bond_type 16  gh-ctd
  # bond_type 17  gh-nucl
  if lhbool:
    molecule.nbond_type = 17
  else:
    molecule.nbond_type = 5

  # angle_type 1 bead-bead (lp)
  # angle_type 2  bead-bead twist, align
  # angle_type 3 wlctwistend for last DNA bead before nucl (or end of molecule)
  # angle_type 4 bead(entry/exit)-nucl orient f relative to r, uses zeta
  # angle_type 5 align nucl u with centerghost
  # angle_type 12 align nucl f with centerghost
  # angle_type 6 angle bewtween entering/exiting DNA and nucl, angle alpha
  # angle_type 7 angle bewtween u of entering DNA and nucl, angle eta
  # angle_type 8 angle bewtween u of exiting DNA and nucl, angle mu
  # angle_type 9 bead entry exit orient f, related to how nucl dna can be twisted
  # angle_type 10 bead(entry/exit)-nucl orient f_dna relative to f_nucl, uses kappa
  # angle_type 11 bead(entry/exit)-nucl orient f relative to r, uses zeta (similar to 4, but angle types are different)
  # angle_type 19 ctd-ctd angle
  # angle_type 13-18 globular head interactions
  if lhbool:
    molecule.nangle_type = 19
  else:
    molecule.nangle_type = 12

  # no dihedrals yet (might use if H1 is present)
  # we are now using dihedrals to preserve the structure of the globular head of the linker histone
  # dihedral_type 1 beads 2-0-3-4 on the globular head
  # dihderal_type 2 beads 2-1-5-4 on the globular head
  # dihderal_type 2 beads 3-6-CTD-CTD to prevent rotation of the GH
  if lhbool:
    molecule.ndihedral_type = 3;
  else:
    molecule.ndihedral_type = 0;

  n = len(molecule.ellipsoids)
  # this allows us to bond the nucleosome to the LH
  nucl_array = []
  inucl = 0
  for a1 in range(n-1):
    ta2 = ta3 = ta4 = ta5 = ta6 = ta7 = -1
    a2 = a3 = a4 = a5 = a6 = a7 = -1

    ta1 = molecule.ellipsoids[a1].mytype
    a2 = a1+1
    ta2 = molecule.ellipsoids[a2].mytype
    if (a2 < n-1):
      a3 = a2+1
      ta3 = molecule.ellipsoids[a3].mytype
    if (a3 < n-1):
      a4 = a3+1
      ta4 = molecule.ellipsoids[a4].mytype
    if (a4 < n-1):
      a5 = a4+1
      ta5 = molecule.ellipsoids[a5].mytype
    if (a5 < n-1):
      a6 = a5+1
      ta6 = molecule.ellipsoids[a6].mytype
    if (a6 < n-1):
      a7 = a6+1
      ta7 = molecule.ellipsoids[a7].mytype

    # bonds
    # note: entire stem must exist. i.e. there must be an entering and exiting DNA, and two ghost sites
    if (ta1 == typemap['ghost']):
        nucl_array.append(a1)
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']):
      molecule.bonds.append(Bond(2,a1,a2)) #dna-dna bond
    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']) and (ta3 == typemap['ghost']) and (ta4 == typemap['bead']):
      molecule.bonds.append(Bond(1,a2,a1)) #dna-nucl bond
      molecule.bonds.append(Bond(1,a2,a4)) #dna-nucl bond
      molecule.bonds.append(Bond(3,a1,a4)) #enterdna-exitdna, length e
      molecule.bonds.append(Bond(4,a3,a1)) #bead-centerghost, length 0.5e
      molecule.bonds.append(Bond(4,a3,a4)) #bead-centerghost, length 0.5e
      molecule.bonds.append(Bond(5,a2,a3)) #bead-centerghost, length 0.5e
    if ((ta1 == typemap['ctd']) or (ta1 == typemap['ghost']) or (ta1 == typemap['bead'])) and (ta2 == typemap['gh']): # all of the below interactions are gh-gh
      molecule.bonds.append(Bond(6,a2,a4))
      molecule.bonds.append(Bond(7,a2,a5))
      molecule.bonds.append(Bond(8,a3,a4))
      molecule.bonds.append(Bond(9,a3,a7))
      molecule.bonds.append(Bond(10,a5,a6))
      molecule.bonds.append(Bond(11,a6,a7))
      molecule.bonds.append(Bond(12,a2,a7))
      molecule.bonds.append(Bond(13,a3,a6))
      molecule.bonds.append(Bond(14,a4,a5))
      molecule.bonds.append(Bond(17,nucl_array[inucl],a2))
      molecule.bonds.append(Bond(17,nucl_array[inucl],a4))
      molecule.bonds.append(Bond(17,nucl_array[inucl],a5))
      inucl = inucl + 1
    if (ta1 == typemap['gh']) and (ta2 == typemap['ctd']): #gh - ctd
      molecule.bonds.append(Bond(16,a1,a2))
    if (ta1 == typemap['ctd']) and (ta2 == typemap['ctd']): #ctd - ctd
      molecule.bonds.append(Bond(15,a1,a2))
    #if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == typemap['bead']) and (ta4 == typemap['nucl']):
    #  molecule.bonds.append(Bond(5,a1,a4)) # morse bead-nucl
    #  molecule.bonds.append(Bond(6,a1-3,a4)) # morse bead-nucl
    #  #molecule.bonds.append(Bond(5,a2,a4)) # morse bead-nucl
    #if (ta1 == typemap['nucl']) and (ta2 == typemap['ghost']) and (ta3 == typemap['bead']) and (ta4 == typemap['bead']) and (ta5 == typemap['bead']):
    #  molecule.bonds.append(Bond(5,a1,a5)) # morse bead-nucl
    #  molecule.bonds.append(Bond(6,a1,a5+3)) # morse bead-nucl
    #  #molecule.bonds.append(Bond(5,a1,a5)) # morse bead-nucl

    # angles
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']):
      molecule.angles.append(Angle(2,a1,a2,a2)) #dna twist
      #molecule.angles.append(Angle(13,a1,a2,a2)) #dna twist

    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']):
      molecule.angles.append(Angle(4,a1,a2,a2)) #orient bead and nucl
      molecule.angles.append(Angle(10,a1,a2,a2)) #orient bead and nucl
      molecule.angles.append(Angle(7,a1,a2,a2)) #orient bead and nucl

    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == typemap['bead']):
      molecule.angles.append(Angle(1,a1,a2,a3)) #dna bend

    if (ta1 == typemap['nucl']) and (ta2 == typemap['ghost']) and (ta3 == typemap['bead']):
      molecule.angles.append(Angle(11,a1,a3,a3)) #orient bead and nucl
      molecule.angles.append(Angle(10,a1,a3,a3)) #orient bead and nucl
      molecule.angles.append(Angle(8,a1,a3,a3)) #orient bead and nucl

    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == typemap['nucl']):
      molecule.angles.append(Angle(3,a1,a2,a2)) #wlctwistend
      molecule.angles.append(Angle(6,a1,a2,a3)) #dna-dna-nucl stem angle
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == -1):
      molecule.angles.append(Angle(3,a1,a2,a2)) #wlctwistend for end of molecule

    if (ta1 == typemap['nucl']) and (ta2 == typemap['ghost']) and (ta3 == typemap['bead']) and (ta4 == typemap['bead']):
      molecule.angles.append(Angle(6,a1,a3,a4)) #nucl-dna-dna stem angle

    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']) and (ta3 == typemap['ghost']) and (ta4 == typemap['bead']):
      molecule.angles.append(Angle(5,a2,a3,a3)) #orient nucl u with center ghost
      molecule.angles.append(Angle(12,a2,a3,a3)) #orient nucl u with center ghost
      molecule.angles.append(Angle(9,a1,a4,a4)) #orient enter exit dna with orient f

    #These are the angles for the globular head of the linker histone
    if ((ta1 == typemap['ctd']) or (ta1 == typemap['ghost']) or (ta1 == typemap['bead'])) and (ta2 == typemap['gh']): # all of the below interactions are gh-gh
      molecule.angles.append(Angle(13,a4,a2,a5)) #gh angle 1
      molecule.angles.append(Angle(14,a2,a5,a6)) #gh angle 2
      molecule.angles.append(Angle(15,a4,a3,a7)) #gh angle 3
      molecule.angles.append(Angle(16,a3,a6,a7)) #gh angle 4
      molecule.angles.append(Angle(17,a2,a7,a6)) #gh angle 5
      molecule.angles.append(Angle(18,a3,a4,a5)) #gh angle 6

    if (ta1 == typemap['ctd']) and (ta2 == typemap['ctd']) and (ta3 == typemap['ctd']):
      molecule.angles.append(Angle(19,a1,a2,a3))

    #dihedrals for the linker histone
    if ((ta1 == typemap['ctd']) or (ta1 == typemap['ghost']) or (ta1 == typemap['bead'])) and (ta2 == typemap['gh']): # all of the below interactions are gh-gh
      molecule.dihedrals.append(Dihedral(1,a4,a2,a5,a6)) #gh dihedral 1
      molecule.dihedrals.append(Dihedral(2,a4,a3,a7,a6)) #gh dihedral 2
    if (ta4 == typemap['gh']) and (ta5 == typemap['ctd']):
      molecule.dihedrals.append(Dihedral(3,a1,a4,a5,a6)) #gh-ctd dihedral

def align_with_1kx5(molecule,param):

    natoms = len(molecule.ellipsoids)

    #find first nucleosome, store its index in inuc
    for i in range(natoms):
        if molecule.ellipsoids[i].mytype == 1:
            inuc = i;
            break

	positions_1kx5={
        9: {'in': (21.001449, 139.245216, 18.249657), 'out': (78.792145, 131.412435, -14.919424), 'nucl': (47.248498, 91.081386, 6.483035)} ,
        10: {'in': (6.672159, 113.122551, 22.479820), 'out': (88.215994, 102.847036, -14.452594), 'nucl': (47.248498, 91.081386, 6.483035)} ,
        11: {'in': (5.231514, 110.051192, 22.080299), 'out': (89.404634, 99.799229, -13.277508), 'nucl': (47.248498, 91.081386, 6.483035)} }

	pos1kx5 = np.zeros((3,3))
	pos1kx5[0] = positions_1kx5[param.nucl_bp_unwrap]['in']
	pos1kx5[1] = positions_1kx5[param.nucl_bp_unwrap]['out']
	pos1kx5[2] = positions_1kx5[param.nucl_bp_unwrap]['nucl']

    # pos, quat, fvu of chosen nucl
    pos0 = molecule.ellipsoids[inuc].pos
    quat0 = molecule.ellipsoids[inuc].quat
    fvu0 = quat_fvu_rot(np.eye(3),quat0)  # update fvu

    #target orientation of first nucl
    # f - aligned with z, u - aligned with y, v - aligend with x
    qgoal = tu2rotquat(-m.pi*2./3.,[1,1,1])

    #rotation quat to apply to all sites
    qrot = quat_multiply(qgoal,quat_conj(quat0))

    # zero position
    for i in range(natoms):
        molecule.ellipsoids[i].pos = molecule.ellipsoids[i].pos - pos0

    # rotate everyone
    for i in range(natoms):
        posnew = quat_vec_rot(molecule.ellipsoids[i].pos,qrot)
        molecule.ellipsoids[i].pos = posnew

        qtmp = molecule.ellipsoids[i].quat
        qnew = quat_multiply(qrot,qtmp)
        molecule.ellipsoids[i].quat = qnew

    #pos0_1kx5 = np.array([46.0960386 ,  90.45177426,   6.50217265]) #1kx5 center of mass
    for i in range(natoms):
        molecule.ellipsoids[i].pos = molecule.ellipsoids[i].pos + pos1kx5[2]

def calculate_nrl_dna_unwrap(param):
  '''Compute NRL and nucleosomal DNA that can unwrap'''

  param.dna_linker_length = param.nrl - param.dna_in_nucl
  param.dna_linker_length_ends = param.nrlends - param.dna_in_nucl

  # how much nucl dna that can 'unbind' depends on the linker length
  #  this is so any arbitrary linker length can be generated keeping a fixed discritization
  if param.nnucleosomes == 1:
    get_nucl_bp_offset = {0:9,1:11,2:10} # this only accounts for one bp_unwrap
  else:
    get_nucl_bp_offset = {0:9,1:10,2:11} # this accounts for two bp_unwrap

  modulo = param.dna_linker_length % param.basepair_per_bead
  param.nucl_bp_unwrap = get_nucl_bp_offset[modulo]

  # correct nrlends if necessary
  modulo = (param.dna_linker_length_ends + param.nucl_bp_unwrap) % param.basepair_per_bead
  l = round((param.dna_linker_length_ends + param.nucl_bp_unwrap) / param.basepair_per_bead) * param.basepair_per_bead
  nrltmp = l - param.nucl_bp_unwrap + param.dna_in_nucl
  ltmp = l - param.nucl_bp_unwrap
  if modulo != 0:
    print "Warning! NRL of ends spefied (%d) was rounded to %d in order to be consistent with param.nucl_bp_unwrap (%d)  obtained from NRL (%d)" % ( param.nrlends, nrltmp, param.nucl_bp_unwrap, param.nrl)
  param.nrlends = nrltmp
  param.dna_linker_length_ends = ltmp


def calculate_geom(geom,param):
  '''Calculate Lengths and Angles corresponding to 1CPN geometry. The image file diagram.svg explains what each of these values are '''
  # first dict value is offset, others are lengths

  lengths_1kx5={
        0: {'c': 55.796923, 'i': 33.169081, 'a': 67.091765} ,
        9: {'c': 48.178119, 'i': 36.932414, 'a': 90.105419} ,
        10: {'c': 48.017333, 'i': 35.357808, 'a': 91.871604} ,
        11: {'c': 47.175048, 'i': 32.370058, 'a': 91.850785} }

  geom.c = lengths_1kx5[param.nucl_bp_unwrap]['c']  # distance of stem site from nucl center
  geom.i = lengths_1kx5[param.nucl_bp_unwrap]['i']  # vertical distance between entering and exiting DNA sites
  geom.a = lengths_1kx5[param.nucl_bp_unwrap]['a']  #real space distance between entering and exiting DNA sites

  geom.f = np.sqrt(geom.a*geom.a - geom.i*geom.i); #horizontal distance between entering and exiting DNA
  geom.g = np.sqrt(geom.c*geom.c - geom.i*geom.i/2.0/2.0)
  geom.h = np.sqrt(geom.g*geom.g - geom.f*geom.f/2.0/2.0)
  geom.epsilon = np.arcsin(0.5*geom.i/geom.c)
  geom.gamma = 2*np.arcsin(0.5*geom.f/geom.g) # then rotate by asin around f to point u at dyad
  print geom.gamma
  # alpha_prime is not alpha, but they're VERY close.
  # For now I'll leave it, but eventually I should do this exactly correct
  #geom.alpha_prime = p_alpha*np.cos(epsilon)
  geom.alpha_prime = geom.alpha
  geom.kappa = (2*np.pi - geom.alpha_prime - geom.gamma) + (0.5*geom.gamma)
  geom.mu = np.pi - geom.alpha_prime - 0.5*geom.gamma
  geom.zeta = np.pi/2.0 + np.arcsin(0.5*geom.i / geom.c) # should be obtuse angle
  geom.eta = 2*np.pi - geom.kappa

  geom.theta = (param.twist_per_bp * param.nucl_bp_unwrap) % (2*np.pi) # angle between f vector of enter/exit DNA and nucl

  #calc j
  if (geom.d < geom.h):
    print "Warning! specified i (%f) is less than calculated h (%f). Setting i=h" % (geom.d,geom.h)
    geom.d = geom.h
  geom.e = np.sqrt((geom.d - geom.h)**2 + (0.5*geom.a)**2)

# Plain Old Data Classes to store data
class Geometry(object):
    __slots__ = ('alpha','beta','gamma','epsilon',
                 'alpha_prime' ,'theta','mu',
                 'zeta','eta','kappa','alpha0',
                 'c','d','e','f','g','h','i','a')
    def __init__(self):
        self.alpha0 = 90.0 * m.pi / 180. # alpha to use in force field

class Parameters(object):
    __slots__ = ('nucl_bp_unwrap',
                'charge_per_nucleosome',
                'bp_mass','charge_per_bp','twist_per_bp','rise_per_bp',
                'bead_mass','charge_per_bead','twist_per_bead','rise_per_bead',
                'nucl_mass',
                'ghost_mass',
                'basepair_per_bead',
                'nrl','nrlends','lh',
                'lnucldna',
                'lengthscale',
                'dna_linker_length','dna_linker_length_ends',
                'dna_in_nucl',
                'bead_shape','nucl_shape','ghost_shape',
                'nucl_nucl_bond_offset',
                'nnucleosomes',
                'directory',)
    def __init__(self):
        self.dna_in_nucl = 147;
        self.directory = "."
        self.basepair_per_bead = 3 # 3 basepair per bead
        self.charge_per_nucleosome = 0;
        self.charge_per_bp = -1;
        self.twist_per_bp = 2.0*m.pi / 10.0 # radian
        self.rise_per_bp = 3.3 # Angstroms

# current lh parameters come from Luque et al (2014)
class LinkerHistone(object):
    __slots__ = ('ctd_mass', 'gh_mass', 'num_in_gh', 'linit',
                'lequil', 'beta', 'salt_scale', 'lnucllh',
                'ctd_shape', 'ctd_charges', 'ctd_beads',
                'ctd_bond_length',)
    def __init_(self):
        self.lequil = 15.0; # bond equil length
        self.ctd_beads = 22; # number of beads in ctd

#Modular function that adds in the linker histones if selected
def add_linker_histones(molecule,lhist,param):

  #Data taken from Lugue et al 2014 and translated to fit our reference system
  gh_data={'gh1': {'pos': np.array([52.228, 38.609,10.7562]), 'charge': -3.29 },
      'gh2': {'pos': np.array([48.939, 30.429, -6.935]), 'charge': 4.22 },
      'gh3': {'pos': np.array([57.043, 32.130, 5.458]), 'charge': 8.48 },
      'gh4': {'pos': np.array([45.038, 37.282, 3.817]), 'charge': 0.28 },
      'gh5': {'pos': np.array([41.897, 39.769, -8.150]), 'charge': 2.08 },
      'gh6': {'pos': np.array([57.580, 41.316, -4.953]), 'charge': 3.27 }}

  #Change the linker histone to be centered at the dyad in our reference system
  #gh_change=np.array([2.546,36.589,0.0069])
  gh_change=np.array([0,0,0])
  ctd_charges = [0,2,2,3,0,4,0,2,2,4,0,2,3,2,2,2,2,2,2,2,2,3]; # array of H1.4 ctd CG charges
  iellipsoid = len(molecule.ellipsoids)

  for ellipsoid in molecule.ellipsoids:
    # after the nucleosome is set we need to add in the linker histone
    if ellipsoid.mytype == 1:
        quat = ellipsoid.quat
        pos  = ellipsoid.pos
        fvu  = quat_fvu_rot(np.eye(3),quat)
        fvu0 = np.eye(3)
        molid = ellipsoid.molid

        # the globular head gets added first
        for igh in range(lhist.num_in_gh):
          # print gh bead to string
          typestr = 'gh%i' % (igh+1)
          lvec = gh_data[typestr]['pos']
          # positions need to be rotated about the axis to match 1CPN notation
          quat_v_rot = tu2rotquat(m.pi/2,fvu0[1])
          # quat_f_rot = tu2rotquat(m.pi*(1.0-54.0/180),fvu0[0])
          # the 54.0 degrees is to transfer from Schlick group nucleosome reference system to 1CPN fvu
          quat_f_rot = tu2rotquat(m.pi*(0.5+54.0/180),fvu0[0])
          quat_fv_rot = quat_multiply(quat_f_rot, quat_v_rot)
          lh_quat = quat_multiply(quat, quat_fv_rot)
          #lh_quat = quat_normalize(lh_quat)
          lh_pos = quat_vec_rot(lvec,lh_quat) # rotate linker vector by nucl quat
          lh_pos = np.add(lh_pos,pos) # translate linker position by nucl pos

          mytype = typemap['gh']
          molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,lh_pos/param.lengthscale, quat,gh_data[typestr]["charge"],lhist.ctd_shape,molid))
          iellipsoid += 1

        # for the other AAs
        for ictd in range(lhist.ctd_beads):
          if ictd == 0:
            lh_pos = np.add(lh_pos,np.multiply(lhist.ctd_bond_length,fvu[2])) #eq dist is 22.64 angstroms
            q = tu2rotquat(lhist.beta,fvu[0]) #rotate by 110 degrees about f
            lh_quat = quat_multiply(q,quat)
            lh_fvu = quat_fvu_rot(fvu0,lh_quat)  # update linker fvu
          else:
            # update the rotation of the next lh bead so that the zigzag pattern appears
            if ictd % 2 == 0:
                q = tu2rotquat(lhist.beta, fvu[0])
            else:
                q = tu2rotquat(-lhist.beta, fvu[0])
            lh_quat = quat_multiply(q,quat)   # update quat
            lh_fvu = quat_fvu_rot(fvu0,lh_quat)  # update fvu

            # set position using f,v,u and pre-calculated lengths
            lh_pos = np.add(lh_pos,np.multiply(lhist.linit, -lh_fvu[2]))

          mytype = typemap["ctd"]
          # the 1.5 is a prefactor that was included in the original model to reproduce mesoscopic chromatin structure
          # the salt_scale is the DiSCO calculated scaling for 150mM
          charge = ctd_charges[ictd]*lhist.salt_scale*1.5
          molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,lh_pos/param.lengthscale,lh_quat,charge,lhist.ctd_shape,molid))
          #print(quat)
          iellipsoid += 1
def main():

  parser = argparse.ArgumentParser()
  #parser.add_argument('outputdir',type=str, help="Directory to write output files")
  parser.add_argument('-a','--stemangle',default=80.0,type=float, help='Angle of stem (in degrees)')
  parser.add_argument('-nrl','--nrl',default=187, type=int, help='Nucleosome repeat length')
  parser.add_argument('-nrlends','--nrlends',type=int, help='Nucleosome repeat length of ends')
  parser.add_argument('-n','--nnucl',default=1,type=int, help='Number of nucleosomes')
  parser.add_argument('-len','--lengthscale',default=1,type=float,help='Normalization length scale. 10 - nanometers, 1 - angstroms, 55 - norm by nucl height')
  parser.add_argument('--nuclbonds',default=0,type=int, help='Generate bonds between neighboring nucleosomes (for joshs cv)')
  parser.add_argument('-d','--d',default=35,type=float,help='length from ghost site to nucl center')
  parser.add_argument('-lh','--linkerhistone',default=False,type=bool,help='turn on linker histones')
  #parser.add_argument('-bp','--bpbead',default=3, type=int, help='Number of basepairs per bead')
  #parser.add_argument('-c','--stemlength',default=48.18,type=float)
  #parser.add_argument('-d','--stemheight',default=36.93,type=float)
  #parser.add_argument('-e','--stemdiag',default=90.11,type=float)
  args = parser.parse_args()

  # structures to store data
  geom = Geometry()
  param = Parameters()

  #------------------------------------------------
  # Input Values
  #------------------------------------------------
  geom.alpha = args.stemangle * m.pi / 180. # alpha that config is initialized to
  param.nnucleosomes = args.nnucl
  p_nucl_nucl_bond_offset = args.nuclbonds; # which nucl do I artificially bond together
  param.lengthscale = args.lengthscale # the normalization factor of output distances
  param.nrl = args.nrl
  param.lh = args.linkerhistone # linker histone boolean
  geom.d = args.d

  if args.nrlends == None:
    param.nrlends = args.nrl
  else:
    param.nrlends = args.nrlends

  if param.lh:
    lhist = LinkerHistone()

  # calculate some parameters
  calculate_nrl_dna_unwrap(param)

  # this is a bit of a hack, if there are zero nucleosomes, just draw DNA length of nrl
  if param.nnucleosomes == 0:
    param.dna_linker_length_ends = param.nrl
    param.nucl_bp_unwrap = 10

  # I dont use this currently, but I might again...
  param.lnucldna = (param.dna_in_nucl - 2*param.nucl_bp_unwrap) * param.rise_per_bp

  # calculate geometrical parameters
  calculate_geom(geom,param)

  print param.nucl_bp_unwrap
  print geom.theta

  #------------------------------------------------
  # Calculate per bead parameters
  #------------------------------------------------
  param.charge_per_bead = param.basepair_per_bead * param.charge_per_bp
  param.twist_per_bead = param.basepair_per_bead * param.twist_per_bp
  param.rise_per_bead = param.basepair_per_bead * param.rise_per_bp

  #------------------------------------------------
  # reference position and orientation
  #------------------------------------------------
  pos0 = np.array([0,0,0])
  quat0 = tu2rotquat(1e-5,[1,0,0])
  q = [0] * 4
  # fvu0 is the reference coordinate system for nucl and beads
  # to get current reference system, just rotate fvu0 by a quat
  fvu0 = np.eye(3)

  #------------------------------------------------
  # Setup molecule
  #------------------------------------------------
  molecule = Molecule()

  # manually specify some molecule parameters
  param.nucl_shape = [55,110,110]
  param.bead_shape = [20,20,20]
  param.ghost_shape = [10,10,10]
  param.bp_mass = 650 #g/mol
  param.nucl_mass = 107616 + param.bp_mass*(147-param.nucl_bp_unwrap)   # Histone = 107616 g/mol, basepairs = 147
  param.bead_mass = param.bp_mass* param.basepair_per_bead               # Mass = MW of bp (650g/mol)*num bp
  param.ghost_mass = 10*param.bead_mass

  # for now the parameters are hard coded
  if param.lh:
      lhist.ctd_shape = [18,18,18]
      lhist.gh_mass = 110.0*80.0/6.0
      lhist.ctd_mass = 110.0*5.0
      # I should move these to just be calculated from the xyz values
      lhist.num_in_gh = 6
      lhist.salt_scale = 1.68000; # salt scaling of ctd charges at 150 mM
      lhist.ctd_beads = 22
      lhist.linit = 7.0; # bond init length
      lhist.lequil = 15.0 # bond equil length
      lhist.lnucllh = 33.0 # ghost-lh length
      lhist.beta = 110.0 * m.pi / 180.; # beta for the ctd
      lhist.ctd_bond_length = 10.0

  # set masses
  molecule.atom_types.append(AtomType(1,param.nucl_mass))
  molecule.atom_types.append(AtomType(2,param.bead_mass))
  molecule.atom_types.append(AtomType(3,param.ghost_mass))
  # check for the linker histone before adding the atom types
  if param.lh:
      molecule.atom_types.append(AtomType(4,lhist.gh_mass))
      molecule.atom_types.append(AtomType(5,lhist.ctd_mass))
  boxl = 5000/param.lengthscale;
  molecule.set_box(-boxl, boxl, -boxl, boxl, -boxl,boxl)



  # ====================================
  # Main generation loop
  # ====================================

  # counts
  iellipsoid = 1

  # INFO
  # for each new site, assume pos, fvu and quat from the previous site are set
  for inuc in range(param.nnucleosomes+1):
    #------------------------------------
    # generate DNA first
    #------------------------------------
    if (inuc == 0):
      pos = pos0
      quat = quat0
      bead_pos0 = pos
      bead_quat0 = quat
    else:
      # set position of bead following nucl
      pos = np.add(pos,-0.5*geom.i*fvu[0])
      pos = np.add(pos,0.5*geom.f*fvu[1])
      pos = np.add(pos,geom.h*fvu[2])

      q = tu2rotquat(geom.mu, fvu[0]) # update orientation
      quat = quat_multiply(q,quat)   # update quat
      fvu = quat_fvu_rot(fvu0,quat)  # update fvu

      #now need to rotate via kappa
      q = tu2rotquat(-geom.theta, fvu[2]) # update orientation
      quat = quat_multiply(q,quat)   # update quat
      fvu = quat_fvu_rot(fvu0,quat)  # update fvu

      bead_pos0 = pos
      bead_quat0 = quat

    # vary number of beads if 1st or last nucl
    if (inuc==0) or (inuc == param.nnucleosomes):
      nbp = param.dna_linker_length_ends + param.nucl_bp_unwrap
    else:
      nbp = param.dna_linker_length + 2* param.nucl_bp_unwrap

    n_bead_in_linker = int(m.ceil(nbp / param.basepair_per_bead))

    #pdb.set_trace()

    for ibead in range(n_bead_in_linker):
      if ibead == 0:
        pos = bead_pos0
        quat = bead_quat0

        #make stem bead (ghost)
        #if (inuc != 0):
        #  fvu = quat_fvu_rot(fvu0,quat)
        #  posstem = np.add(pos,p_d*fvu[0]);
        #  mytype = typemap['ghost']
        #  molid = inuc #first bead is part of prev nucl molecule
        #  molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,posstem/p_lengthscale,quat,0,ghost_shape,molid))
        #  iellipsoid += 1

        mytype = typemap['bead']
        if (inuc != 0):
          molid = inuc #first bead is part of prev nucl molecule
        else:
          molid = 0      #all non-nucl dna gets molid 0
        molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,pos/param.lengthscale,quat,param.charge_per_bead,param.bead_shape,molid))
        iellipsoid += 1

      else:
        fvu = quat_fvu_rot(fvu0,quat)
        pos = np.add(pos,np.multiply(param.rise_per_bead,fvu[2]))    # increment position along u
        q = tu2rotquat(param.twist_per_bead, fvu[2]) # rotate around u
        quat = quat_multiply(q,quat) # update quat

        mytype = typemap['bead']
        if (inuc != param.nnucleosomes) and (ibead == (n_bead_in_linker - 1)):
          molid = inuc+1 #last bead is part of nucl molecule
        else:
          molid = 0 #all non-nucl dna gets molid 0
        molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,pos/param.lengthscale,quat,param.charge_per_bead,param.bead_shape,molid))
        iellipsoid += 1



    #------------------------------------
    # now generate nucleosome
    #------------------------------------
    if (inuc < param.nnucleosomes):

      #make stem bead (ghost)
      #fvu = quat_fvu_rot(fvu0,quat)
      #posstem = np.add(pos,-p_d*fvu[0]);
      #mytype = typemap['ghost']
      #molid = inuc+1 #first bead is part of nucl molecule
      #molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,posstem/p_lengthscale,quat,0,ghost_shape,molid))
      #iellipsoid += 1

      # set nucleosome position
      fvu = quat_fvu_rot(fvu0,quat) # update fvu

      # rotate kappa around u
      q = tu2rotquat(-geom.theta, fvu[2])
      quat = quat_multiply(q,quat)   # update quat
      fvu = quat_fvu_rot(fvu0,quat)  # update fvu

      # rotate first, this makes position update easier
      q = tu2rotquat(geom.kappa, fvu[0])
      quat = quat_multiply(q,quat)   # update quat
      fvu = quat_fvu_rot(fvu0,quat)  # update fvu

      # set position using f,v,u and pre-calculated lengths
      pos = np.add(pos,-0.5*geom.i*fvu[0])
      pos = np.add(pos,0.5*geom.f*fvu[1])
      pos = np.add(pos,-geom.h*fvu[2])


      mytype = typemap['nucl']
      molid = inuc + 1
      molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,pos/param.lengthscale,quat,param.charge_per_nucleosome,param.nucl_shape,molid))
      iellipsoid += 1

      # make the central stem bead
      fvu = quat_fvu_rot(fvu0,quat)
      posstem = np.add(pos,geom.d*fvu[2]);
      mytype = typemap['ghost']
      molid = inuc+1 #first bead is part of nucl molecule
      molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,posstem/param.lengthscale,quat,0,param.ghost_shape,molid))
      iellipsoid += 1

  # add in the linker histones
  if param.lh:
    add_linker_histones(molecule,lhist,param)

  # call funciton to set all bonded interactions
  set_bonded_interactions(molecule,param.lh)

  # function to position first nucleosome at origin
  if (param.nnucleosomes != 0):
    align_with_1kx5(molecule,param)

  # josh for bonding all nucl together
  if p_nucl_nucl_bond_offset:
    gen_nucl_nucl_bonds(molecule,param.nucl_nucl_bond_offset)

  if (not os.path.exists(param.directory)):
    os.mkdir(param.directory)

  write_lammps_variables('in.variables',param,geom)
  molecule.write_dump("in.dump")
  #molecule.write_xyz("in.xyz")
  molecule.write_lammps_input("in.lammps")
  write_args(args,"%s/in.args" %   param.directory)
  molecule.write_psf("%s/in.psf" % param.directory)

if __name__ == "__main__":
    main()
