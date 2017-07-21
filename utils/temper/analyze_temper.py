#!/usr/bin/env python

import sys
import numpy as np
import pdb

if len(sys.argv) != 2:
    print "Usage: %s <lammps temper file>" % sys.argv[0]
    sys.exit(1)

filename = sys.argv[1]

nlines = sum(1 for line in open(filename))
f = open(filename,'r')

#parse up until tempering starts
found = False
for line in f:
    if "Step T0" in line:
        found = True
        break

if not found:
  print "Error! This doesn't look like a tempering file!"
  sys.exit(1)

# get number of replicas 
l = line.split()
nreplica = len(l)-1
print "nreplica: %d" % nreplica

first = True
T = np.zeros((nreplica,))
nswap = np.zeros((nreplica,))
nattempt = 0
nroundtrips = 0
direction_Tlow = 1   
direction_Thigh = -1
print "need to compute nroundtrips!"
for line in f:
    l = line.split()
    if (len(l) == nreplica +1) and (l[0] != "Step"):
      timestep = int(l[0])
      for r in range(nreplica):
         T[r] = int(l[r+1])

      if (not first) :
         nswap += (T != Tprev)

      Tprev = np.copy(T)
      
      # compute round trips!
      if (nreplica-1) == T[0]:  #if TMax reached boxMin
        if direction_Tlow == -1: # if moving left
            direction_Tlow = 1   # shift to moving right
        elif direction_Tlow == 1: # if moving right
            direction_Tlow = 1    # move left
            nroundtrips += 1      # and was roundtrip!

      if 0 == T[nreplica-1]:  #if TMin reached boxMax
        if direction_Tlow == 1: # if moving right
            direction_Tlow = -1   # shift to moving left
        elif direction_Tlow == -1: # if moving left
            direction_Tlow = 1    # move right
            nroundtrips += 1      # and was roundtrip!


      nattempt += 1
      first = False
    else:
      break 
print "Finished replica exchange section of file after %d attempted exchanges over %d timesteps" % (nattempt, timestep)
f.close()

print "Swap Probabilities: "
swapprob = nswap / nattempt
print swapprob

print "Number of Roundtrips: "
print nroundtrips
