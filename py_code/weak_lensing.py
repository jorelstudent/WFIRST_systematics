#The purpose of this code is to get the cov matrix
#For just a specific subset


import numpy as np
import pdb

#Load the data
#inv_cov_data = np.loadtxt('/Users/penafiel/JPL/nora_data/cov_3x2pt_4.500000e+01_1.800000e+04_WFIRST_Ncl15_Ntomo10_2pt_inv')
#x = inv_cov_data[:,0]
#y = inv_cov_data[:,1]
#inv_cov = inv_cov_data[:,2]


f = open('/Users/penafiel/JPL/nora_data/cov_3x2pt_4.500000e+01_1.800000e+04_WFIRST_Ncl15_Ntomo10_2pt_inv')
#Fill our inverse covariance matrix
inv_cov = np.zeros((1650,1650))
for line in f:
	i,j, val = line.split()
	i = int(i)
	j = int(j)
	inv_cov[i][j] = val

#Invert to get covariance
#Strip to only get the weak lensing ones
#Then re-invert

cov = np.linalg.inv(inv_cov)
cov_825 = cov[:825,:825]
inv_cov_825 = np.linalg.inv(cov_825)

#This array should be (1,2,3, ...., 1,2,3....)
y = range(len(inv_cov_825)) * len(inv_cov_825)

#This array should be (0,0,0,0,....,1,1,1,....)
x = []
for i in range(len(inv_cov_825)):
	x = np.append(x, [i] * len(inv_cov_825))

#Reshape our inverse covariance matrix
inv_cov_825_col = np.reshape(inv_cov_825, 825*825)

np.savetxt('inv_cov_wl.dat', np.transpose([x,y,inv_cov_825_col]))


pdb.set_trace()

	
