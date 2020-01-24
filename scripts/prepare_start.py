import os
import math
import sys
print "*.py monomer_name outputname num"
monomer = sys.argv[1] + '.pdb'
outputfile = sys.argv[2] + '.pdb'
num = int(sys.argv[3])
chain = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
with open(monomer,'r') as fopen:
   lines = fopen.readlines()
distz = 0
distx = 0
disty = 0
for line in lines:
  if line.split()[0] != 'END':
    if abs(float(line[30:38])) > distx:
       distx = abs(float(line[30:37]))
    if abs(float(line[38:46])) > disty:
       disty = abs(float(line[38:45]))
    if abs(float(line[46:54])) > distz:
       distz = abs(float(line[46:53]))
       resi_num = int(line.split()[1])
       relic = line[11:26]

#distz = 100.000/(num+1)
#distx = 0
#disty = 0
distx += 100.0
disty += 100.0
distz += 100.0
data = ''
for i in range(num):
    for line in lines:
     if  line.split()[0] != 'END':
         dx = (i//3 + i%3//1 )*distx
         dy = (i//3 + i%3//2 )*disty
         dz = (i//3)*distz
         #print dx,dy,dz   
         x= format(float(line[30:37]) + (i//3 + i%3//1 )*distx,'.3f')
         y= format(float(line[38:45]) + (i//3 + i%3//2 )*disty,'.3f')
         z= format(float(line[46:53]) + (i//3)*distz,'.3f')
    #     x= format(float(line[30:37]) ,'.3f')
    #     y= format(float(line[38:45]) ,'.3f')
    #     z= format(float(line[46:53]) + i*distz,'.3f')
         single = line[0:21] + chain[i] + line[22:30] + str(x).rjust(8,' ') + str(y).rjust(8,' ') + str(z).rjust(8,' ') + line[54:80] + '\n'
     else:
         single = 'TER    ' + str(resi_num + 1)   + relic[0:10] + chain[i] + relic[11:15] + '\n'
     data += single
    print i
    print i//3
    print chain[i]
    print dx,dy,dz
with open(outputfile,'w') as fwrite:
     fwrite.writelines(data)
