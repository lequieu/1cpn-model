#!/usr/bin/env python 

import numpy as np
import os, glob, sys
import subprocess
import math
import copy

def ener_conserve():

    if not os.path.isfile("energy.dat"):
        print("No energy file cannot test system")
        return True

    d=np.loadtxt("energy.dat")

    erot_threshold = 1e-6
    d_low_erot = d[d[:,1] < erot_threshold]
    n = d_low_erot.shape[0]

    tol = 1e-1
    etot0 =d[0,2]
    for i in range(n):
        etot = d_low_erot[i,2]
        if (np.abs(etot - etot0) > tol):
            print "TEST FAILED! etot at step %d is %f!" %(d_low_erot[i,0], etot)
            return False
    return True

def mom_conserve():
    if not os.path.isfile("traj-angmom.dump"):
        print("No momentum file found. Cannot test system")
        return True

    f = open("traj-angmom.dump",'r')
    lines = f.readlines()

    momentum = 0.0
    momentum_0 = 0.0
    flag = 0
    for line in lines:
        q = line.split()
        if len(q) == 12:
            #Calculate system and individual angular momentum
            #L = R x P + sum(I_i * w_i)
            mass = float(q[8])
            if flag == 0:
                momentum_0 = momentum_0 + mass*(float(q[9]) + float(q[10]) + float(q[11]))
            else:
                momentum = momentum + mass*(float(q[9]) + float(q[10]) + float(q[11]))
        # Check after each frame for the system
        elif 'TIMESTEP' in line:

            if abs(momentum-momentum_0) > 1E-3:
                print("Warning system is not conserving momentum (%f)\n" % mom)
                return False
            #Reinitialize everything
            momentum = 0.0

            if flag == 0:
                flag = 1
    return True

def amom_conserve():
    if not os.path.isfile("traj-angmom.dump"):
        print("No momentum file found. Cannot test system")
        return True

    f = open("traj-angmom.dump",'r')
    lines = f.readlines()

    ang_mom = 0.0
    mtot = 0
    natoms = int(lines[3])
    com = np.zeros(3)
    v_com = np.zeros(3)
    v = []
    r = []
    flag = 0
    for line in lines:
        q = line.split()
        if len(q) == 12:
            #Calculate system and individual angular momentum
            #L = R x P + sum(I_i * w_i)
            ang_mom = ang_mom + float(q[5]) + float(q[6]) + float(q[7])
            mass = float(q[8])
            v.append(np.array([float(q[9]),float(q[10]),float(q[11])]))
            r.append(np.array([float(q[2]),float(q[3]),float(q[4])]))
            com = com + r[-1]*mass
            v_com = v_com + v[-1]*mass
            mtot = mtot + mass
        # Check after each frame for the system
        elif 'TIMESTEP' in line and mtot != 0:
            com = com*1.0/mtot
            v_com = v_com*1.0/mtot
            #Calculate initial angular momentum
            if flag == 0:
                ang_mom_0 = ang_mom
                for i in range(0,natoms):
                    ang_mom_0 = ang_mom_0 + np.sum(mass*np.cross(r[i]-com,v[i]-v_com))
                ang_mom = ang_mom_0
                flag = 1
            else:
                for i in range(0,natoms):
                    ang_mom = ang_mom + np.sum(mass*np.cross(r[i]-com,v[i]-v_com))
            if abs(ang_mom-ang_mom_0) > 1E-3:
                print("Warning system is not conserving angular momentum (%f)\n" % ang_mom)
                return False
            #Reinitialize everything
            mtot = 0.0
            com = 0.0*com
            v_com = 0.0*v_com
            ang_mom = 0.0
            v = []
            r = []
    return True

def onecpn_test(indir,lammpsfile):
    if not os.path.isfile("%s" % (lammpsfile)):
        print("Executable file (%s) does not exist" % lammpsfile)
        sys.exit(1)

    #Run lammps
    stdout = open('stdout.txt','wb')
    stderr = open('stderr.txt','wb')
    subprocess.call(["lmp_mpi", "-in", lammpsfile],stdout=stdout,stderr=stderr)

    ener_test = ener_conserve()
    mom_test = mom_conserve()
    amom_test = amom_conserve()

    if ener_test and mom_test and amom_test:
        print("%s conserves necessary observable quantities. Continuing..." % indir)
    else:
        if not ener_test:
            print("%s does not conserve energy" % indir)
            sys.exit(2)
        if not mom_test:
            print("Warning: %s does not conserve momentum" % indir)
            sys.exit(3)
        if not amom_test:
            print("Warning: %s does not conserve angular momentum" % indir)
            #sys.exit(4)

print("Running 1CPN integration tests")

#Currently only providing this option for the mpi version of lammps
#Will probably throw an option to change later
lammps_exec = "lmp_mpi"
infile = "./ffld_list.dat"

##Check for the force field list
##Taken out for now because a walk might be cleaner
#if not os.path.isfile(infile):
#    print("Force field list file not found.")
#    sys.exit(1)

#Run all tests in the directory
owd = os.getcwd()
for (dirpath,dirnames,filenames) in os.walk("."):
    for file in filenames:
        if 'in.zewdie' in file or 'in.wlc' in file or 'in.1cpn' in file:
            os.chdir(dirpath)
            test = onecpn_test(dirpath,file)
            os.chdir(owd)

print("Tests have now finished")
