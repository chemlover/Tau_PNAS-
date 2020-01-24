import os
import sys
import numpy as np 
pmffile = sys.argv[1]
perturbfile = sys.argv[2]
pmf = np.loadtxt(pmffile)
perturb = np.loadtxt(perturbfile)
for i in range(len(perturb)):
    pmf[i][1] += perturb[i]
np.savetxt(pmffile,pmf,fmt="%.6f",delimiter=" ")
