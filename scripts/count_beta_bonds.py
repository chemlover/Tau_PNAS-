#!/usr/bin/python
# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian
#
# Mod by Zheng 09/2014####
# Function:
# Designed for Abeta oligomers
# Count the number of intra & inter-monomer beta hydrogen bonds
# Output: Intra_count, Inter_count
# ----------------------------------------------------------------------

import sys
import Zheng_func

if len(sys.argv)!=7:
	print "\n.py mono.pdb dump_file Output_file Nmer len_linker HB_CA_cutoff(recommend:5.36)\n"
	print
	sys.exit()

pdb_file=Zheng_func.pdbname(sys.argv[1])
dump_file = sys.argv[2]
out_file = sys.argv[3]
Nmer=     int(sys.argv[4])
len_linker =int(sys.argv[5]) 
cutoff=float(sys.argv[6])

n_chains=Zheng_func.get_n_chains(pdb_file)
if n_chains > 1 :
	print "mono.pdb should have only one chain!"
	print " Exit with Error."
	sys.exit()

ca_atoms_pdb=Zheng_func.get_ca(pdb_file)
file_len=Zheng_func.get_file_len(dump_file)
print file_len
line = Zheng_func.getline_range(dump_file, 4, 4); line=line[0].split() ; #get number of atoms
n_atoms=int(line[0]) ; nline=n_atoms + 9
print n_atoms

n_snapshot=file_len / nline
out = open(out_file, 'w')
for i in range(n_snapshot):
	line_start = 1 + nline*i ; line_end = nline*(i+1)
	dump_part=Zheng_func.getline_range(dump_file, line_start, line_end)
	ca_atoms=Zheng_func.get_ca_dump(dump_part)
   	assert(len(ca_atoms) == n_atoms/3)
	#native_count, nonnative_count = Zheng_func.fused_count_contact(ca_atoms_pdb, ca_atoms, Nmer, len_linker, cutoff)
	#OUTPUB [intra_count, inter_count] = [(anti-para HB, para HB), (anti-para HB, para HB)]
	intra_count, inter_count = Zheng_func.count_beta_contact(ca_atoms_pdb, ca_atoms, Nmer, len_linker, cutoff) #recommend cutoff=5.36 A

	str1=' '.join(map(str,intra_count))
	str2=' '.join(map(str,inter_count))
	#for item in inter_count :
	#	str0=' '.join(map(str,item))
	#	str2+=' '+str0
        #print intra_count
        #print inter_count
	#print str1+str2
	out.write(str1+' '+str2)
	out.write('\n')

out.close()
