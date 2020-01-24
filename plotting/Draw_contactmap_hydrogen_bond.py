import numpy as np
import pandas
import prody
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import cm as cm
import sys
import math
import os
np.set_printoptions(threshold=np.nan)


def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    return math.sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))
def checkIfNative(xyz_CAi, xyz_CAj):
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r< 8: return True
    else: return False

def isNative(r):
        if r<12.0: return True
        else: return False

def get_ca_atoms(pdb_part):
    ca_atoms = []
    o_atoms = []
    n_atoms = []
    for line in pdb_part: 
        if line.split()[2] == 'N':
           x=float(line[30:38])
           y=float(line[38:46])
           z=float(line[46:54])
           #atom = [x,y,z]
           #n_atoms.append(atom)
        else:
           atom = [0,0,0]
           n_atoms.append(atom)
        break
    print line
    print n_atoms          
    for line in pdb_part:
       if len(line.split()) > 5:
        #if line.split()[0] != 'END':
           if (line.split()[2] == 'CA') and line.split()[0] == "ATOM":
           #print line
              x=float(line[30:38])
              y=float(line[38:46])
              z=float(line[46:54])
              atom = [x,y,z]
              ca_atoms.append(atom)
           if line.split()[2] == 'N':
              x=float(line[30:38])
              y=float(line[38:46])
              z=float(line[46:54])
              atom = [x,y,z]
              n_atoms.append(atom)
           if line.split()[2] == 'O':
              x=float(line[30:38])
              y=float(line[38:46])
              z=float(line[46:54])
              atom = [x,y,z]
              o_atoms.append(atom)
    h_atoms = [[100,100,100]]
    for i in range(1,len(n_atoms)):
        #print i
        x = 0.84100 * ca_atoms[i-1][0] + 0.89296 * ca_atoms[i][0] - 0.73389 * o_atoms[i-1][0]
        y = 0.84100 * ca_atoms[i-1][1] + 0.89296 * ca_atoms[i][1] - 0.73389 * o_atoms[i-1][1]
        z = 0.84100 * ca_atoms[i-1][2] + 0.89296 * ca_atoms[i][2] - 0.73389 * o_atoms[i-1][2]
        atom = [x,y,z]
        h_atoms.append(atom)
        #print h_atoms[i],ca_atoms[i],o_atoms[i]
    return ca_atoms,o_atoms,n_atoms,h_atoms

def compute_contactmap(ca_atoms,o_atoms,n_atoms,h_atoms):
    width = len(ca_atoms)
    contactmap = np.zeros((width,width))
    cutoff = 0.8
    for i in range(width):
        for j in range(i+4,width):
            #if ca_atoms[i] != [0,0,0] and ca_atoms[j] != [0,0,0]: 
               #if checkIfNative(ca_atoms[i], ca_atoms[j]):
                ron = vabs(vector(o_atoms[j],n_atoms[i]))
                roh = vabs(vector(o_atoms[j],h_atoms[i]))
                #print ron,roh
                stone = math.exp(-(ron-2.98)*(ron-2.98)/2/0.68/0.68 - (roh-2.06)*(roh-2.06)/2/0.76/0.76)
                #print stone
                if stone > cutoff:
                  contactmap[i][j] = 1.0
                  contactmap[j][i] = contactmap[i][j]
    #print  np.max(contactmap)
    return contactmap

import glob
def convert_hb_contactmap(contactmap):
    [m,n] = np.shape(contactmap)
    new_contactmap = np.zeros((m,n))
    for i in range(1,m-1):
        for j in range(0,n-1):
            if contactmap[i][j] ==  1 and (contactmap[i-1][j+1] == 1 or contactmap[i+1][j+1] == 1): 
               new_contactmap[i][j] = 1
               #print new_contactmap[i][j]
               #print i,j
    #print new_contactmap
    return new_contactmap 
def get_average_contactmap(n_res):
    path = './'
    filelist = glob.glob(path+'*.pdb')
    contactmaps = np.zeros((n_res,n_res))
    s = 0
    for a in filelist:
        with open(a,"r") as fopen:
             pdb_part = fopen.readlines()
        ca_atoms,o_atoms,n_atoms,h_atoms = get_ca_atoms(pdb_part)
        if len(ca_atoms) == n_res:
           contactmap =  compute_contactmap(ca_atoms,o_atoms,n_atoms,h_atoms)
           #contactmap = convert_hb_contactmap(contactmap)
           contactmaps = contactmaps + contactmap
           s += 1
    #print s
    contactmap = np.true_divide(contactmaps,s)
    return contactmap

def separate_contactmap_inter_intra(contactmap,os):
    [m,n] = np.shape(contactmap)
    n_res = m
    n_res_chain = n_res / os
    contactmap_inter = np.zeros((n_res_chain,n_res_chain))
    contactmap_intra = np.zeros((n_res_chain,n_res_chain))
    n_intra = os
    for i in range(os):
        delta = i*n_res_chain
        for l in range(n_res_chain):
            for m in range(n_res_chain):
                contactmap_intra[l][m] += contactmap[l+delta][m+delta]/os
    n_inter = 0
    for i in range(os-1):
        for j in range(i+1,os):
            n_inter += 1
    for i in range(os-1):
        for j in range(i+1,os):
            delta_i = i*n_res_chain
            delta_j = j*n_res_chain
            for l in range(n_res_chain):
                for m in range(n_res_chain):
                    contactmap_inter[l][m] += (contactmap[l+delta_i][m+delta_j]+contactmap[l+delta_j][m+delta_i])/os/2
    return contactmap_intra,contactmap_inter
     

def draw_contactmap(matrix,xname,yname,title,maplim):
    (n_atomsx,n_atomsy) = np.shape(matrix)
    if n_atomsx == n_atomsy:
       n_atoms = n_atomsx
    else:
       print  "matrix is not square"
       sys.exit()
    plt.figure(figsize=(10,10))
    X=np.arange(1,n_atoms+1);
    cmap = cm.get_cmap('viridis')
    cmap = cm.get_cmap('jet')
    Y=np.arange(1,n_atoms+1)
    #n = 30
    #x = 0.0125
    #white = plt.cm.seismic(np.ones(80-79)*x)
    #upper = plt.cm.seismic(np.linspace(1-x, 1, n))
    #colors = np.vstack((white, upper))
    #tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)
    #cmap = plt.cm.jet  # define the colormap
    # extract all colors from the .jet map
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    ## force the first color entry to be grey
    #cmaplist[0] = (1, 1, 1, 1.0)
    # create the new map
    #cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    #print matrix
    #cmax = np.max(matrix)*1.2
    #for i in range(int(n_atoms)):
    #  for j in range(int(n_atoms)):
      #   print matrix
      #   print i,j,matrix[i][j]
         #if matrix[int(i)][int(j)] == 0.0:
         #    matrix[i][j]=float("inf")
    plt.pcolormesh(X,Y,matrix,vmin=0,vmax=1,cmap=cmap)
    #,edgecolors='k')
#colorbar;
    plt.colorbar()
    maplim = np.max(matrix)
    plt.clim(0,maplim)
    #plt.xlabel(xname,fontsize=30)
    hfont = {'fontname':'helvetica'}
    plt.ylabel(yname,fontsize=30,**hfont)
    plt.title(title,fontsize=30,**hfont)
    plt.xlim([1,n_atoms])
    plt.ylim([1,n_atoms])
    major_ticks = np.arange(1, n_atoms+1,8)
    minor_ticks = np.arange(1, n_atoms+1,8)
    ax = plt.axes()
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    #ax.set_yticks(major_ticks)
    #ax.set_yticks(minor_ticks, minor=True)
#plt.grid()
    #ax.grid(which='both')
    #ax.grid(which='minor', alpha=0.2)
    #ax.grid(which='major', alpha=0.5)
    #ax.set_yticks([1,709,1007,1176,1441], minor=False)
    #x.set_yticks([1,709,1007,1176,1441], minor=True)
    #ax.yaxis.grid(True, which='major')
    #ax.yaxis.grid(True, which='minor')
    #ax.set_xticks([1,709,1007,1176,1441], minor=False)
    #ax.set_xticks([1,709,1007,1176,1441], minor=True)
    #ax.xaxis.grid(True, which='major')
    #ax.xaxis.grid(True, which='minor')
    #plt.axis([0,n_atoms+1, 0, n_atoms+1])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    np.savetxt(yname+".txt",matrix,fmt="%.6f",delimiter=" ")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('%s.png'%(yname))
    #plt.show()

def main():
    #if len(sys.argv) != 5:
    #   print "*py pdbfile1 pdbfile xname yname title"
    #   sys.exit()
    #pdbfile1 = sys.argv[1]
    #pdbfile2 = sys.argv[2]
    title = sys.argv[1]
    n_res_chain = int(sys.argv[2])
    n_chain = int(sys.argv[3])
    maplim = float(sys.argv[4])
    titlename = sys.argv[5]
    n_res = n_res_chain*n_chain
    contactmap = get_average_contactmap(n_res)
    contactmap_intra,contactmap_inter = separate_contactmap_inter_intra(contactmap,n_chain)
    #intername = title + "-interchain"
    #intraname = title + "-intrachain"
    xname = "residue"
    #yname = "residue"
    draw_contactmap(contactmap_intra,xname,"intrachain","",maplim)
    draw_contactmap(contactmap_inter,xname,"interchain",titlename,maplim)
if __name__ == '__main__':
    main()

