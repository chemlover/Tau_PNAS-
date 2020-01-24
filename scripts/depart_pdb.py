import os
import sys
import numpy as np
qfile = sys.argv[1]
q=np.loadtxt(qfile)
qmax =np.max(q)
qmin =np.min(q)
dq=(qmax-qmin)/2
os.system("mkdir fibril prefibril")
for i in range(len(q)):
    if q[i]-qmin <= dq:
       order = i+1
       os.system("cp %d.pdb prefibril/"%(order))
    else:
       order = i+1
       os.system("cp %d.pdb fibril/"%(order))
