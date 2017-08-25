import numpy as np
from collections import OrderedDict
import pdb

parameters = np.genfromtxt('/Users/penafiel/JPL/nora_data/parameter_samples_2')

#Create a dictionary for the parameters, since it is easier to call
param_dict = {}

param_list = ['Z01', 'Z02', 'Z03', 'Z04', 'Z05', 'Z06', 'Z07', 'Z08', 'Z09', 'Z10', 'SZ1',
			  'M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10',
			  'AIA', 'BIA', 'EIA', 'EIH']
param_ind = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
			 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
			 50, 51, 52, 53]

#The parameters are varied 2 at a time
#x = 1st parameter varied
#y = 2nd parameter varied
max_par = len(param_list) - 1 # Max number of plots on the x/y axis
ind_1 = 0 #beginning index where the parameters are located
ind_2 = 100 #ending index where the parameters are located
x_ax = 0 # for the x parameter

for i in range(max_par):
	y_ax = np.copy(i) + 1 # For the y parameter
	while y_ax <= max_par:#(x_ax < y_ax) & (x_ax<max_par):
		x_string = param_list[x_ax]
		y_string = param_list[y_ax]
		
		print x_string, param_ind[x_ax], y_string, param_ind[y_ax]
		print ind_1, ind_2
		param_dict[x_string + '_' + y_string + '_x'] = parameters[ind_1:ind_2, param_ind[x_ax]]
		param_dict[x_string + '_' + y_string + '_y'] = parameters[ind_1:ind_2, param_ind[y_ax]]
		ind_1 += 100
		ind_2 += 100
		y_ax += 1
	x_ax += 1


#Now we have our dictionary!
pdb.set_trace()

