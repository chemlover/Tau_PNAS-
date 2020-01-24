import os
import re
import sys
import numpy as np
print "*py inputfile_id output_id length_single_chain repeat_already num_add"

PdbID = sys.argv[1]
outputfile = sys.argv[2]
length= int(sys.argv[3])
repeat = int(sys.argv[4])
num_add = int(sys.argv[5])
atom_type = ['N','CA','C','O','CB','H']
attach_type = ['N','C','C','O','C','H']
#Pdb_type = ['ALA','VAL','LEU','ILE','PRO','PHF','TRP','MET','GLY','SER','THR','CYS','TYR','ASN','GLN','ASP','GLU','LYS','ARG','HIS']
chain_type = ['X','Y','Z',]
#first = np.zeros((5*(num_add+repeat)*length),7)1
Pdb_seq = []
with open(PdbID+".pdb",'r') as fopen:
  lines = fopen.readlines()
  for line in lines:
 #      N = 0 
#       print line
       if line.split()[0] == "TER":
 #         N += 1
         break
       if line.split()[0] != "TER" and line.split()[0] != "END":
           if line.split()[2] == 'CA':
              Pdb_seq.append(line[17:20])
 #         
first = np.zeros((5*length,3))
end = np.zeros((5*length,3))
dist = np.zeros((5*length,3))
N_lin = 0
for line in lines:      
       if line.split()[0] == "TER":
          Finish =line
          break
       else:
          N_seq = int(line[22:26])
          atomtype = line.split()[2]
       for i in range(5): 
               
        if atomtype == atom_type[i] or (atomtype == 'H' and line[17:20] == 'GLY'):
    #     if line[17:20] == 'GLY' or  (atomtype == 'H' and line[17:28] == 'GLY'):
          
    #        if atomtype == 'CB':
    #           atomtype = 'H'   
          number = N_seq*5-5+i 
          if atomtype == 'H':
           number = N_seq*5-1 
#           print line
          first[number][0] = float(line[30:38])
          first[number][1] = float(line[38:46])
          first[number][2] = float(line[46:54]) 
N_chain = 0
N_lin = 0
for line in lines:
       if line.split()[0] == "TER":
          N_chain += 1
       if N_chain == repeat - 1:  
          if line.split()[0] != 'END':
             N_seq = int(line[22:26])
             atomtype = line.split()[2]
#             print atomtype
             if atomtype == 'H' and line[17:20] == 'GLY':
                print line
             for i in range(5):
              if atomtype == atom_type[i] or  (atomtype == 'H' and line[17:20] == 'GLY'): 
  #             if line[17:20] == 'GLY':
  #               if atomtype == 'CB':
  #                  atomtype = 'H'       
               number = N_seq*5-5+i
               if atomtype == 'H':
                number = N_seq*5-1
  
               end[number][0] = float(line[30:38])
               end[number][1] = float(line[38:46])
               end[number][2] = float(line[46:54])

for i in range(length*5):
     print i
#     print dist[i] 
     dist[i][0] = (end[i][0]-first[i][0])/(repeat-1.0)
     dist[i][1] = (end[i][1]-first[i][1])/(repeat-1.0)
     dist[i][2] = (end[i][2]-first[i][2])/(repeat-1.0)
     print dist[i]
data = ''
#print dist
for line in lines:
    if line.split()[0] != 'END':
       data += line
#       print line
       N_new_line = int(line.split()[1])
       TER_chain_id = line[21:22]
N_new_line += 1
data +=  'TER ' + str(N_new_line).rjust(7,' ') + '      ' + Pdb_seq[length-1]+ ' ' + TER_chain_id + str(length).rjust(4,' ' ) + '\n'
for i in range(num_add):
   for j in range(length):
     for k in range(5):
      N_new_line += 1
      squeue = j*5+ k
      x = format(end[squeue][0] + (i+1)*dist[squeue][0],'.3f')
      y = format(end[squeue][1] + (i+1)*dist[squeue][1],'.3f')
      z = format(end[squeue][2] + (i+1)*dist[squeue][2],'.3f')
      if Pdb_seq[j] == 'GLY' and atom_type[k] == 'CB':
        k = k + 1
      data += 'ATOM' + str(N_new_line).rjust(7,' ') + '  ' + atom_type[k].ljust(4,' ') +  Pdb_seq[j] + ' ' + chain_type[i] + str(j+1).rjust(4,' ') + str(x).rjust(12,' ')  + str(y).rjust(8,' ') + str(z).rjust(8,' ') + '  1.00  0.00' + attach_type[k].rjust(12,' ') + '\n'
   N_new_line += 1
   if i != num_add -1 :
     data += 'TER ' + str(N_new_line).rjust(7,' ') + '      ' + Pdb_seq[length-1]+ ' ' + chain_type[i] + str(length).rjust(4,' ' ) + '\n'
   else:
     data += 'END'
with open(outputfile+'.pdb','w')  as fwrite:
     fwrite.writelines(data)    
