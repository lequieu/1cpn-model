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

from lhistone import *
from molecule import *
from vect_quat_util import *

# mape site names to numbers
typemap = {'nucl': 1,
           'bead': 2,
           'dyad': 3,
           'gh' : 4,
           'ctd' : 5}

def write_args(args,param,geom,fnme):
   f=open(fnme,"w")
   command = ' '.join(sys.argv)
   f.write('# Inputs issued to command line\n')
   f.write("%s\n" %command)

   f.write('\n# All Arguements\n')
   f.write("%s\n" %args)

   f.write("\n# All paramaters and geometry values (just in case)\n")
   for name in param.__slots__:
      value = eval("param.%s" % name)
      if type(value) == int:
        f.write("param.%s = %d\n" %(name,value))
      elif type(value) == float or type(value) == np.float64:
        f.write("param.%s = %f\n" %(name,value))
      elif type(value) == str:
        f.write("param.%s = %s\n" %(name,value))
      elif type(value) == bool:
        f.write("param.%s = %d\n" %(name,value))
      elif (type(value) == list) and (len(value) == 3):
        f.write("param.%s = [%f,%f,%f]\n" %(name,value[0],value[1],value[2]))
      else:
        print "Variable %s has invalid output type %s" % (name, type(value))

   for name in geom.__slots__:
      value = eval("geom.%s" % name)
      if type(value) == int:
        f.write("geom.%s = %d\n" %(name,value))
      elif type(value) == float or type(value) == np.float64:
        f.write("geom.%s = %f\n" %(name,value))
      elif type(value) == str:
        f.write("geom.%s = %s\n" %(name,value))
      elif type(value) == bool:
        f.write("geom.%s = %d\n" %(name,value))
      elif (type(value) == list) and (len(value) == 3):
        f.write("geom.%s = [%f,%f,%f]\n" %(name,value[0],value[1],value[2]))
      else:
        print "Variable %s has invalid output type %s" % (name, type(value))


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
  fnew.write( "variable omega0 equal %f\n" % (param.twist_per_bead*180.0/np.pi))
  fnew.write( "variable mu equal %f\n" % (geom.mu*180.0/np.pi))
  fnew.write( "variable zeta equal %f\n" % (geom.zeta*180.0/np.pi))
  fnew.write( "variable eta equal %f\n" % (geom.eta*180.0/np.pi))
  fnew.write( "variable theta equal %f\n" % (geom.theta*180.0/np.pi))

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

def set_bonded_interactions(molecule,lhbool):
  # =================
  # key of types
  # =================
  # atom_type 1  nucleosome
  # atom_type 2  bead
  # atom_type 3  stem (dyad site)
  # atom_type 4  globular head of linker histone
  # atom_type 5  c-terminal domain of linker histone
  if lhbool:
      molecule.natom_type = 5
  else:
      molecule.natom_type = 3

  # bond_type 1  bead-nucleosome, length c
  # bond_type 2  bead-bead
  # bond_type 3  bead-bead (enter and exit nucl), length e
  # bond_type 4  bead-centerdyad, length j
  # bond_type 5  nucl-centerdyad, length i
  # bond_type 6-14 globular head interactions
  # bond_type 15  ctd-ctd
  # bond_type 16  gh-ctd
  # bond_type 17  gh-dyad
  # bond_type 18  gh-nucl
  if lhbool:
    molecule.nbond_type = 17
  else:
    molecule.nbond_type = 5

  # angle_type 1 bead-bead (lp)
  # angle_type 2  bead-bead twist, align
  # angle_type 3 wlctwistend for last DNA bead before nucl (or end of molecule)
  # angle_type 4 bead(entry/exit)-nucl orient f relative to r, uses zeta
  # angle_type 5 align nucl u with centerdyad
  # angle_type 12 align nucl f with centerdyad
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
  # dihderal_type 3 beads 3-6-CTD-CTD to prevent rotation of the GH - REMOVED
  if lhbool:
    molecule.ndihedral_type = 2;
  else:
    molecule.ndihedral_type = 0;

  n = len(molecule.ellipsoids)
  # this allows us to bond the nucleosome to the LH
  nucl_array = []
  dyad_array = []
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
    #if (a6 < n-1):
    #  a7 = a6+1
    #  ta7 = molecule.ellipsoids[a7].mytype

    # bonds
    # note: entire stem must exist. i.e. there must be an entering and exiting DNA, and two dyad sites
    if (ta1 == typemap['nucl']):
        nucl_array.append(a1)
        dyad_array.append(a2)
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']):
      molecule.bonds.append(Bond(2,a1,a2)) #dna-dna bond
    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']) and (ta3 == typemap['dyad']) and (ta4 == typemap['bead']):
      molecule.bonds.append(Bond(1,a2,a1)) #dna-nucl bond
      molecule.bonds.append(Bond(1,a2,a4)) #dna-nucl bond
      molecule.bonds.append(Bond(3,a1,a4)) #enterdna-exitdna, length e
      molecule.bonds.append(Bond(4,a3,a1)) #bead-centerdyad, length 0.5e
      molecule.bonds.append(Bond(4,a3,a4)) #bead-centerdyad, length 0.5e
      molecule.bonds.append(Bond(5,a2,a3)) #bead-centerdyad, length 0.5e
    if ((ta1 == typemap['gh']) and (ta6 == typemap['gh'])): # all of the below interactions are gh-gh
      molecule.bonds.append(Bond(6,a1,a3))
      molecule.bonds.append(Bond(7,a1,a4))
      molecule.bonds.append(Bond(8,a2,a3))
      molecule.bonds.append(Bond(9,a2,a6))
      molecule.bonds.append(Bond(10,a4,a5))
      molecule.bonds.append(Bond(11,a5,a6))
      molecule.bonds.append(Bond(12,a1,a6))
      molecule.bonds.append(Bond(13,a2,a5))
      molecule.bonds.append(Bond(14,a3,a4))
      # three bonds for both the dyad and nucl keep the LH bound to the dyad position
      # any removed and you get off-dyad binding
      molecule.bonds.append(Bond(16,dyad_array[inucl],a1)) # dyad-GH1
      molecule.bonds.append(Bond(16,dyad_array[inucl],a3)) # dyad-GH3
      molecule.bonds.append(Bond(16,dyad_array[inucl],a4)) # dyad-GH4
      molecule.bonds.append(Bond(17,nucl_array[inucl],a1)) # nucl-GH1
      molecule.bonds.append(Bond(17,nucl_array[inucl],a3)) # nucl-GH3
      molecule.bonds.append(Bond(17,nucl_array[inucl],a4)) # nucl-GH4
    if (ta1 == typemap['gh']) and (ta2 == typemap['ctd']): #gh - ctd
      molecule.bonds.append(Bond(15,a1,a2))
    if (ta1 == typemap['ctd']) and (ta2 == typemap['ctd']): #ctd - ctd
      molecule.bonds.append(Bond(15,a1,a2))

    # angles
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']):
      molecule.angles.append(Angle(2,a1,a2,a2)) #dna twist

    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']):
      molecule.angles.append(Angle(4,a1,a2,a2)) #orient bead and nucl
      molecule.angles.append(Angle(10,a1,a2,a2)) #orient bead and nucl
      molecule.angles.append(Angle(7,a1,a2,a2)) #orient bead and nucl

    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == typemap['bead']):
      molecule.angles.append(Angle(1,a1,a2,a3)) #dna bend

    if (ta1 == typemap['nucl']) and (ta2 == typemap['dyad']) and (ta3 == typemap['bead']):
      molecule.angles.append(Angle(11,a1,a3,a3)) #orient bead and nucl
      molecule.angles.append(Angle(10,a1,a3,a3)) #orient bead and nucl
      molecule.angles.append(Angle(8,a1,a3,a3)) #orient bead and nucl

    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == typemap['nucl']):
      molecule.angles.append(Angle(3,a1,a2,a2)) #wlctwistend
      molecule.angles.append(Angle(6,a1,a2,a3)) #dna-dna-nucl stem angle
    if (ta1 == typemap['bead']) and (ta2 == typemap['bead']) and (ta3 == -1):
      molecule.angles.append(Angle(3,a1,a2,a2)) #wlctwistend for end of molecule

    if (ta1 == typemap['nucl']) and (ta2 == typemap['dyad']) and (ta3 == typemap['bead']) and (ta4 == typemap['bead']):
      molecule.angles.append(Angle(6,a1,a3,a4)) #nucl-dna-dna stem angle

    if (ta1 == typemap['bead']) and (ta2 == typemap['nucl']) and (ta3 == typemap['dyad']) and (ta4 == typemap['bead']):
      molecule.angles.append(Angle(5,a2,a3,a3)) #orient nucl u with center dyad
      molecule.angles.append(Angle(12,a2,a3,a3)) #orient nucl u with center dyad
      molecule.angles.append(Angle(9,a1,a4,a4)) #orient enter exit dna with orient f

    #These are the angles for the globular head of the linker histone
    if ((ta1 == typemap['gh']) and (ta6 == typemap['gh'])):
      molecule.angles.append(Angle(13,a3,a1,a4)) #gh angle 1
      molecule.angles.append(Angle(14,a1,a4,a5)) #gh angle 2
      molecule.angles.append(Angle(15,a3,a2,a6)) #gh angle 3
      molecule.angles.append(Angle(16,a2,a5,a6)) #gh angle 4
      molecule.angles.append(Angle(17,a1,a6,a5)) #gh angle 5
      molecule.angles.append(Angle(18,a2,a3,a4)) #gh angle 6

    if (ta1 == typemap['ctd']) and (ta2 == typemap['ctd']) and (ta3 == typemap['ctd']):
      molecule.angles.append(Angle(19,a1,a2,a3))

    #dihedrals for the linker histone
    if ((ta1 == typemap['gh']) and (ta6 == typemap['gh'])):
      molecule.dihedrals.append(Dihedral(1,a3,a1,a4,a5)) #gh dihedral 1
      molecule.dihedrals.append(Dihedral(2,a3,a2,a6,a5)) #gh dihedral 2
    # taking away the third dihedral
    if (ta1 == typemap['gh']) and (ta2 == typemap['ctd']):
      #molecule.dihedrals.append(Dihedral(3,nucl_array[inucl],dyad_array[inucl],a1,a2)) #gh-ctd dihedral
      # the nucl updater is moved to here
      inucl = inucl + 1

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
    print "Note: For nrl=%d, n_bp_unwrap=%d, nrlends= %d (changed from %d!)" % (param.nrl, param.nucl_bp_unwrap, nrltmp,param.nrlends)
  else:
    print "Note: For nrl=%d, n_bp_unwrap=%d, nrlends= %d (as specified)" % (param.nrl, param.nucl_bp_unwrap, nrltmp)

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
  #print geom.gamma
  # alpha_prime is not alpha, but they're VERY close. Since this value is only the initial config (and not used in potential), this approximation is okay
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
    __slots__ = ('alpha','gamma','epsilon',
                 'alpha_prime' ,'theta','mu',
                 'zeta','eta','kappa','alpha0',
                 'c','d','e','f','g','h','i','a')
    def __init__(self):
        self.alpha0 = 90.0 * m.pi / 180. # alpha to use in force field
        self.d = 35.0  #length from dyad site to nucl center

class Parameters(object):
    __slots__ = ('nucl_bp_unwrap',
                'charge_per_nucleosome',
                'bp_mass','charge_per_bp','twist_per_bp','rise_per_bp',
                'bead_mass','charge_per_bead','twist_per_bead','rise_per_bead',
                'nucl_mass',
                'dyad_mass',
                'basepair_per_bead',
                'nrl','nrlends','lh','salt',
                'lnucldna',
                'lengthscale',
                'dna_linker_length','dna_linker_length_ends',
                'dna_in_nucl',
                'bead_shape','nucl_shape','dyad_shape',
                'nnucleosomes','boxsize_factor',
                'directory',)
    def __init__(self):
        # some hardcoded parameters
        self.dna_in_nucl = 147;
        self.directory = "."
        self.basepair_per_bead = 3 # 3 basepair per bead
        self.charge_per_nucleosome = 0;
        self.charge_per_bp = -1;
        self.twist_per_bp = 2.0*m.pi / 10.0 # radian
        self.rise_per_bp = 3.3 # Angstroms
        self.lengthscale = 1.0 #args.lengthscale # the normalization factor of output distances, 10 - nanometers, 1 - angstroms, 55 - norm by nucl height

def main():

  parser = argparse.ArgumentParser()
  #parser.add_argument('outputdir',type=str, help="Directory to write output files")
  parser.add_argument('-a','--stemangle',default=90.0,type=float, help='Angle of stem (in degrees)')
  parser.add_argument('-nrl','--nrl',default=187, type=int, help='Nucleosome repeat length')
  parser.add_argument('-nrlends','--nrlends',type=int, help='Nucleosome repeat length of ends')
  parser.add_argument('-n','--nnucl',default=1,type=int, help='Number of nucleosomes')
  parser.add_argument('-lh','--linkerhistone',dest='lh',action='store_true',help='turn on linker histones')
  parser.set_defaults(lhbool=False)
  parser.add_argument('-salt','--salt',default=150,type=float,help='salt to use for linker histone')
  parser.add_argument('-boxfactor','--boxfactor',default=1.5,type=float,help='Factor to scale the simulation box relative to the initial molecule size')
  args = parser.parse_args()

  # structures to store data
  geom = Geometry()
  param = Parameters()

  #------------------------------------------------
  # Input Values
  #------------------------------------------------
  geom.alpha = args.stemangle * m.pi / 180. # alpha that config is initialized to
  param.nnucleosomes = args.nnucl
  param.nrl = args.nrl
  param.lh = args.lh # linker histone boolean
  param.salt = args.salt;
  param.boxsize_factor = args.boxfactor #scaling factor of boxsize relative to molecule initialization size

  if args.nrlends == None:
    param.nrlends = args.nrl
  else:
    param.nrlends = args.nrlends

  if param.lh:
    lhist = LinkerHistone(geom.d)

  # calculate some parameters
  calculate_nrl_dna_unwrap(param)

  # If there are zero nucleosomes, just draw DNA length of nrl
  if param.nnucleosomes == 0:
    param.dna_linker_length_ends = param.nrl
    param.nucl_bp_unwrap = 10

  # I dont use this currently, but I might again...
  param.lnucldna = (param.dna_in_nucl - 2*param.nucl_bp_unwrap) * param.rise_per_bp

  # calculate geometrical parameters
  calculate_geom(geom,param)

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
  param.dyad_shape = [10,10,10]
  param.bp_mass = 650 #g/mol
  param.nucl_mass = 107616 + param.bp_mass*(147-param.nucl_bp_unwrap)   # Histone = 107616 g/mol, basepairs = 147
  param.bead_mass = param.bp_mass* param.basepair_per_bead               # Mass = MW of bp (650g/mol)*num bp
  param.dyad_mass = 10*param.bead_mass

  # for now the parameters are hard coded
  if param.lh:
      lhist.ctd_shape = [18,18,18]
      lhist.gh_mass = 110.0*80.0/6.0
      lhist.ctd_mass = 110.0*5.0

  # set masses
  molecule.atom_types.append(AtomType(1,param.nucl_mass))
  molecule.atom_types.append(AtomType(2,param.bead_mass))
  molecule.atom_types.append(AtomType(3,param.dyad_mass))
  # check for the linker histone before adding the atom types
  if param.lh:
      molecule.atom_types.append(AtomType(4,lhist.gh_mass))
      molecule.atom_types.append(AtomType(5,lhist.ctd_mass))


  # ====================================
  # Main generation loop
  # ====================================

  # counts
  iellipsoid = 1

  # Note: Each new site assumes that pos, fvu and quat from the previous site are set
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

    for ibead in range(n_bead_in_linker):
      if ibead == 0:
        pos = bead_pos0
        quat = bead_quat0

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


      # make nucleosome site
      mytype = typemap['nucl']
      molid = inuc + 1
      molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,pos/param.lengthscale,quat,param.charge_per_nucleosome,param.nucl_shape,molid))
      iellipsoid += 1

      # make the dyad site
      fvu = quat_fvu_rot(fvu0,quat)
      posstem = np.add(pos,geom.d*fvu[2]);
      mytype = typemap['dyad']
      molid = inuc+1 #first bead is part of nucl molecule
      molecule.ellipsoids.append(Ellipsoid(iellipsoid,mytype,posstem/param.lengthscale,quat,0,param.dyad_shape,molid))
      iellipsoid += 1

  # add in the linker histones
  if param.lh:
    add_linker_histones(molecule,lhist,param)

  # call funciton to set all bonded interactions
  set_bonded_interactions(molecule,param.lh)

  # function to position first nucleosome at origin
  if (param.nnucleosomes == 1):
    align_with_1kx5(molecule,param)

  # zero center of mass to (0,0,0)
  molecule.set_com((0,0,0))

  #auto set box size
  molecule.auto_set_box(param.boxsize_factor);

  if (not os.path.exists(param.directory)):
    os.mkdir(param.directory)

  if param.lh:
      write_lhist_variables('in.var-lh',lhist,param.salt)
  write_lammps_variables('in.variables',param,geom)
  molecule.write_dump("in.dump")
  #molecule.write_xyz("in.xyz")
  molecule.write_lammps_input("in.lammps")
  write_args(args, param, geom, "%s/in.args" %   param.directory)
  molecule.write_psf("%s/in.psf" % param.directory)

if __name__ == "__main__":
    main()
