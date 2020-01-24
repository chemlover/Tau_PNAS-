#!/usr/bin/python
# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian
# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/
# Modified by Weihua Zheng 03/2014 
# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import networkx as nx
import sys
import Zheng_func
from VectorAlgebra import *
import math
def IsConnected(a,b,len,ca_atom,boxsize,cutoff):
	contact=0;
	for i in range(len*a-len,len*a):
		for j in range(len*b-len,len*b):
			xi=ca_atom[i][0];yi=ca_atom[i][1];zi=ca_atom[i][2];
			xj=ca_atom[j][0];yj=ca_atom[j][1];zj=ca_atom[j][2];
			#print xi;
			distij=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj);
			#print math.sqrt(distij)
			if distij <= (cutoff)*(cutoff):
				contact=contact+1;
			distij=0;
	if contact >= 10:
		return 1;
	else:
		return 0;


def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return 180*dihedral_angle(v1, v2, v3)/3.14159265358979

if len(sys.argv)!=5 :
	print "\n.py dump_file output ChainNo ChainLength\n"
	print "Take dump_file, output psi for all residues in each snapshot.\n"
	print
	sys.exit()

dump_file = sys.argv[1]
out_file  = sys.argv[2]
FRETNo  = int (sys.argv[3])
ChainLength =int (sys.argv[4])
file_out = open(out_file,'w');

an = 0.4831806                                                                                                        
bn = 0.7032820                                                                                                        
cn = -0.1864262                                                                                                       
ap = 0.4436538                                                                                                        
bp = 0.2352006                                                                                                        
cp = 0.3211455

#Get n_snapshot
file_len       = Zheng_func.get_file_len(dump_file)
nline_snapshot = Zheng_func.get_nline_snapshot(dump_file)
n_snapshot     = file_len / nline_snapshot
print "n_snapshot of the dump file is = ", n_snapshot


#o=open(out_file,'w')
for i in range(n_snapshot):
	## i_snapshot of the dump_file
	i_dump  = Zheng_func.get_dump_i(dump_file, i) 
	## get atom coordinates
	ca_atoms, cb_atoms, o_atoms= Zheng_func.get_atoms_dump(i_dump)
	n_res = len(ca_atoms)
        n_chain = n_res // ChainLength
        for i in range(n_chain-1):
            for j in range(i+1,n_chain):
                res_i = i*ChainLength + FRETNo - 1
                res_j = j*ChainLength + FRETNo - 1
                distance = vabs(vector(cb_atoms[res_i],cb_atoms[res_j]))
                file_out.write(str(distance)+'\n')
file_out.close()
