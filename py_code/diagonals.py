#The purpose of this code is to just get the diagonals
#of the covariance matrix
#Don't want to load a 1650x1650 array every time
#That thing is 60 MB, not efficient for web interface to have to wait
#for that thing to load

import numpy as np
import pdb

cov_data = np.loadtxt('/Users/penafiel/JPL/nora_data/cov_3x2pt_4.500000e+01_1.800000e+04_WFIRST_Ncl15_Ntomo10_2pt_inv')
x = cov_data[:,0]
y = cov_data[:,1]
cov = cov_data[:,2]

aux = (x==y)
cov_diag = cov[aux]
index = range(1650)

np.savetxt('cov_diag.dat', np.transpose([index,cov_diag]))
pdb.set_trace()


