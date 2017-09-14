
#!/usr/bin/python
import numpy as np
import math as m
import sys
from vect_quat_util import *

class Molecule:
    def __init__(self):
        self.ellipsoids = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.atom_types = []
        self.bond_types = []
        self.angle_types = []
        self.dihedral_types = []
        self.natom_type = 0;
        self.nbond_type = 0;
        self.nangle_type = 0;
        self.ndihedral_type = 0;

    def set_box(self,xlo,xhi,ylo,yhi,zlo,zhi):
        self.boxl = [xlo,xhi,ylo,yhi,zlo,zhi]
    def auto_set_box(self,factor):
        ''' Automatically Set Box Size to square that wraps using atom positions '''
        mymin = 0 #np.zeros((3,))
        mymax = 0 #np.zeros((3,))
        for i in range(len(self.ellipsoids)):
            pos = self.ellipsoids[i].pos
            for j in range(3):
                if pos[j] > mymax: mymax = pos[j]
                if pos[j] < mymin: mymin = pos[j]

        mymin *= factor
        mymax *= factor
        #self.set_box(mymin[0],mymax[0],mymin[1],mymax[1],mymin[2],mymax[2])
        self.set_box(mymin,mymax,mymin,mymax,mymin,mymax)

    def set_com(self,pos):
        ''' shift all positions so that COM is at "pos" '''
        com = self.get_com()
        shift = com + pos
        for i in range(len(self.ellipsoids)):
            self.ellipsoids[i].pos -= shift


    def get_com(self):
        com = np.zeros((3,))
        masstotal = 0;
        for i in range(len(self.ellipsoids)):
            mytype = self.ellipsoids[i].mytype
            mass = self.atom_types[mytype-1].mass
            com += (self.ellipsoids[i].pos * mass)
            masstotal += mass
        com /= masstotal
        return com


    def write_psf(self,fnme):
        f = open(fnme,'w')
        for i in range(5):
          f.write("*\n")
        f.write("\n")
        f.write("\t%d !NATOMS\n" % len(self.ellipsoids))
        for i in range(len(self.ellipsoids)):
            #f.write("%8d %s %-4s %s %-5s %3s %13s %7s 0 \n" % (i+1, "MOL1", self.residue_index[i], self.residue_name3[i], self.site_name1[i], self.site_index[i], self.charge[i], self.mass[i]))
            if self.ellipsoids[i].mytype == 1:
              f.write("%8d %s %-4s %s %-5s %3s %13s %7s 0 \n" % (i+1, "MOL1", "RESI", "RESN", "H",1,self.ellipsoids[i].charge, 100.0))
            elif self.ellipsoids[i].mytype == 2:
              f.write("%8d %s %-4s %s %-5s %3s %13s %7s 0 \n" % (i+1, "MOL1", "RESI", "RESN", "He",1,self.ellipsoids[i].charge, 100.0))
            # right now i am just setting a basic type for the linker histone
            elif self.ellipsoids[i].mytype > 2:
              f.write("%8d %s %-4s %s %-5s %3s %13s %7s 0 \n" % (i+1, "MOL1", "RESI", "RESN", "He",1,self.ellipsoids[i].charge, 100.0))

        nbond = len(self.bonds)
        f.write("\n%8d !NBONDS\n" % (nbond))
        count = 0
        while count < nbond:
            for i in range(4):
                f.write("%8d%8d" % (int(self.bonds[count].atom1+1), int(self.bonds[count].atom2+1)))
                count += 1
                if count >= nbond:
                    break
            f.write("\n")

        f.close()

        pass

    def read_dump(self,fnme):
        nlines = sum(1 for line in open(fnme))
        f = open(fnme)
        natoms = -1
        x = []
        q = []
        mytype = []
        timestep = -1
        xlo = xhi = ylo = yhi = zlo = zhi = -1

        for i in range(nlines):
            line = f.readline()
            if "NUMBER OF ATOMS" in line:
                line=f.readline()
                l = line.split()
                n = int(l[0])
                natoms = n
                x = np.zeros((n,3))
                q = np.zeros((n,4))
                mytype = np.zeros((n,))
            elif "TIMESTEP" in line:
                l = f.readline().split()
                timestep = int(l[0])

            elif "ITEM: BOX BOUNDS" in line:
                l = f.readline().split()
                xlo = float(l[0])
                xhi = float(l[1])
                l = f.readline().split()
                ylo = float(l[0])
                yhi = float(l[1])
                l = f.readline().split()
                zlo = float(l[0])
                zhi = float(l[1])

            elif "ITEM: ATOMS" in line:
                if natoms == 0:
                    print "Error! no natoms!"
                    exit(1)
                for i in range(natoms):
                    line=f.readline()
                    l = line.split()
                    idx = int(l[0])-1
                    mytype[idx] = int(l[1])
                    x[idx][0] = float(l[2])
                    x[idx][1] = float(l[3])
                    x[idx][2] = float(l[4])
                    q[idx][0] = float(l[5])
                    q[idx][1] = float(l[6])
                    q[idx][2] = float(l[7])
                    q[idx][3] = float(l[8])
        f.close()

        #now put into molecule
        for i in range(natoms):
            self.ellipsoids.append(Ellipsoid(i+1,int(mytype[i]),x[i],q[i],0,[0,0,0],0))

    def write_dump(self,fnme):
        f=open(fnme,"w")
        f.write("ITEM: TIMESTEP\n")
        f.write("0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write("%d\n" % len(self.ellipsoids))
        f.write("ITEM: BOX BOUNDS\n")
        f.write("%f %f\n" % (self.boxl[0],self.boxl[1]))
        f.write("%f %f\n" % (self.boxl[0],self.boxl[1]))
        f.write("%f %f\n" % (self.boxl[0],self.boxl[1]))
        f.write("ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4]\n")
        for i in range(len(self.ellipsoids)):
            index = self.ellipsoids[i].index
            mytype = self.ellipsoids[i].mytype
            x = self.ellipsoids[i].pos[0]
            y = self.ellipsoids[i].pos[1]
            z = self.ellipsoids[i].pos[2]
            q0 = self.ellipsoids[i].quat[0]
            q1 = self.ellipsoids[i].quat[1]
            q2 = self.ellipsoids[i].quat[2]
            q3 = self.ellipsoids[i].quat[3]
            f.write("%d %d %f %f %f %f %f %f %f\n" % (index,mytype,x,y,z,q0,q1,q2,q3))
        f.close()
    def write_xyz(self,fnme):
        shift = 5
        n = len(self.ellipsoids)
        f=open(fnme,"w")
        f.write("%d\n\n" % (4.0*n))
        for i in range(n):
            #write position
            pos0 = self.ellipsoids[i].pos
            pos = pos0
            mytype = self.ellipsoids[i].mytype
            #name = inv_typemap[mytype]
            #f.write("%s\t%f %f %f\n" % (name,pos[0],pos[1],pos[2]))
            f.write("%d\t%f %f %f\n" % (mytype,pos[0],pos[1],pos[2]))

            #write orientation
            q = self.ellipsoids[i].quat
            fvu = quat_fvu_rot(np.eye(3),q)
            mymap = {0:'F', 1:'V', 2:'U'}
            for j in range(3):
                pos = np.add(pos0,np.multiply(shift,fvu[j]))
                f.write("%s\t%f %f %f\n" % (mymap[j],pos[0],pos[1],pos[2]))
        f.close()
    def check_types(self):
        # check that all atom,bond,angle types are in range
        natom = len(self.ellipsoids)
        nbond = len(self.bonds)
        nangle = len(self.angles)
        ndihedral = len(self.dihedrals)

        for i in range(natom):
            mytype = self.ellipsoids[i].mytype
            if (mytype > self.natom_type):
                print "Error! Atom of type %d > natomtypes %d" % (mytype,self.natom_type);
                sys.exit(1)

        for i in range(nbond):
            mytype = self.bonds[i].mytype
            if (mytype > self.nbond_type):
                print "Error! bond of type %d > nbondtypes %d" % (mytype,self.nbond_type);
                sys.exit(1)

        for i in range(nangle):
            mytype = self.angles[i].mytype
            if (mytype > self.nangle_type):
                print "Error! angle of type %d > nangletypes %d" % (mytype,self.nangle_type);
                sys.exit(1)

        for i in range(ndihedral):
            mytype = self.dihedrals[i].mytype
            if (mytype > self.ndihedral_type):
                print "Error! dihedral of type %d > ndihedraltypes %d" % (mytype,self.ndihedral_type);
                sys.exit(1)


    def write_lammps_input(self, fnme):
        self.check_types();
        natom = len(self.ellipsoids)
        nbond = len(self.bonds)
        nangle = len(self.angles)
        ndihedral = len(self.dihedrals)
        #do this more elegantly eventually (want molecule.py to be stand alone)

        f=open(fnme,"w")
        f.write("LAMMPS data file\n\n")
        f.write("%d atoms\n" % natom)
        f.write("%d ellipsoids\n" % natom)
        f.write("%d bonds\n" % nbond)
        f.write("%d angles\n" % nangle)
        f.write("%d dihedrals\n" % ndihedral)
        f.write("\n")

        f.write("%d atom types\n" % self.natom_type)
        f.write("%d bond types\n" % self.nbond_type)
        f.write("%d angle types\n" % self.nangle_type)
        f.write("%d dihedral types\n" % self.ndihedral_type)
        f.write("\n")

        f.write("%f %f xlo xhi\n" % (self.boxl[0],self.boxl[1]))
        f.write("%f %f ylo yhi\n" % (self.boxl[2],self.boxl[3]))
        f.write("%f %f zlo zhi\n" % (self.boxl[4],self.boxl[5]))
        f.write("\n")

        f.write("Masses\n\n")
        for i in range(self.natom_type):
            mytype = self.atom_types[i].mytype
            mass = self.atom_types[i].mass
            f.write("%d %f\n" % (mytype, mass))
        f.write("\n")

        f.write("Atoms\n\n")
        for i in range(natom):
            mytype = self.ellipsoids[i].mytype
            pos = self.ellipsoids[i].pos
            #sha = self.ellipsoids[i].shape
            ellipsoidflag = 1
            density = 1
            molid = self.ellipsoids[i].molid + 1
            cha = self.ellipsoids[i].charge
            f.write("%d %d %f %f %f %d %f %d %f\n" % (i+1,mytype,pos[0],pos[1],pos[2],ellipsoidflag,density,molid,cha))
        f.write("\n")

        f.write("Ellipsoids\n\n")
        for i in range(natom):
            q = self.ellipsoids[i].quat
            s = self.ellipsoids[i].shape
            f.write("%d %f %f %f %f %f %f %f\n" % (i+1,s[0],s[1],s[2],q[0],q[1],q[2],q[3]))
        f.write("\n")

        if (nbond != 0):
          f.write("Bonds\n\n")
        for i in range(nbond):
          mytype=self.bonds[i].mytype
          a1 = self.bonds[i].atom1
          a2 = self.bonds[i].atom2
          f.write("%d %d %d %d\n" % (i+1,mytype,a1+1,a2+1))
        f.write("\n")

        if (nangle != 0):
          f.write("Angles\n\n")
        for i in range(nangle):
          mytype=self.angles[i].mytype
          a1 = self.angles[i].atom1
          a2 = self.angles[i].atom2
          a3 = self.angles[i].atom3
          f.write("%d %d %d %d %d\n" % (i+1,mytype,a1+1,a2+1,a3+1))
        f.write("\n")

        if (ndihedral != 0):
          f.write("Dihedrals\n\n")
        for i in range(ndihedral):
          mytype=self.dihedrals[i].mytype
          a1 = self.dihedrals[i].atom1
          a2 = self.dihedrals[i].atom2
          a3 = self.dihedrals[i].atom3
          a4 = self.dihedrals[i].atom4
          f.write("%d %d %d %d %d %d\n" % (i+1,mytype,a1+1,a2+1,a3+1,a4+1))
        f.write("\n")

        f.close()

class Ellipsoid:
    def __init__(self,index,mytype,pos,quat,charge,shape,molid):
        self.index = index
        self.mytype = mytype
        self.pos = pos

        #check quaternions to make sure they have correct length
        norm = quat_norm(quat)
        if (abs(norm - 1.0) > 1e-5):
            print "Error! Norm of ellipse not equal to 1! You messed something up with your quat math!"
            sys.exit(1)

        self.quat = quat
        self.charge = charge
        self.shape = shape
        self.molid = molid

class Bond:
    def __init__(self,mytype,atom1,atom2):
        #self.index = index
        self.mytype = mytype
        self.atom1 = atom1
        self.atom2 = atom2

class Angle:
    def __init__(self,mytype,atom1,atom2,atom3):
        #self.index = index
        self.mytype = mytype
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
class Dihedral:
    def __init__(self,mytype,atom1,atom2,atom3,atom4):
        #self.index = index
        self.mytype = mytype
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
class AtomType:
    def __init__(self,mytype,mass):
       self.mytype = mytype
       self.mass = mass


