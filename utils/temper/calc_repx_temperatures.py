#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
    print "Usage: %s <n> <T_0> <T_f>" % sys.argv[0]
    sys.exit(1)

n = int(sys.argv[1])
T0 = float(sys.argv[2])
Tf = float(sys.argv[3])
R = (Tf/T0)**(1./(n-1.))
print "variable T world",
for i in range(n):
    print "%0.2f" % (T0*R**i),
print ''
print "variable rep world",
for i in range(n):
    print "%d" % (i),
print ''
    

