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
ChainNo  = int (sys.argv[3])
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


o=open(out_file,'w')
file_contact=open("contact_number",'w');
for i in range(n_snapshot):
	## i_snapshot of the dump_file
	i_dump  = Zheng_func.get_dump_i(dump_file, i) 
	## get atom coordinates
	ca_atoms, cb_atoms, o_atoms= Zheng_func.get_atoms_dump(i_dump)
	G=nx.Graph();
	G.add_nodes_from([0,1,2,3,4,5]); 
	for j in range(ChainNo):
		for k in range(j, ChainNo):
			if IsConnected(j,k,ChainLength,ca_atoms,500,10)==1:
				G.add_edge(j,k);	
	
	kkk=sorted(nx.connected_components(G) , key = len , reverse = True);
	#out_file.write(kkk);
	setlen = len(kkk[0]);
	a=kkk[0];### this is also a set object, should be further dealt with.
	#index=[list(a)[0],list(a)[0],list(a)[0],list(a)[]]
	index=[];new_cb_atom=[];
	for jj in range(setlen):
		index.append(list(a)[jj]);
		for kk in range(20):
			new_cb_atom.append(cb_atoms[ChainLength*jj+kk])
	ContactNumber=Zheng_func.compute_N_contacts(new_cb_atom, 6.5, 2)
	#print ContactNumber;
	file_contact.write(str(ContactNumber)+'\n');
	file_out.write(str(len(kkk[0]))+'\n')

o.close()
