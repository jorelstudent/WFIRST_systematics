from flask import Flask, render_template, request
import pandas as pd
from bokeh.embed import components

from bokeh.palettes import Spectral6
from bokeh.layouts import column, widgetbox, WidgetBox, layout
from bokeh.models import CustomJS, Button, HoverTool, ColumnDataSource, LinearColorMapper, BasicTicker, PrintfTickFormatter, ColorBar, OpenURL, TapTool#For the button and a hovertool
from bokeh.models.widgets import Paragraph, PreText, CheckboxGroup, Slider, Dropdown, Select, RangeSlider #For the sliders and dropdown
from bokeh.plotting import figure, curdoc, show
from bokeh.io import gridplot, output_file, show #allows you to make gridplots
from bokeh.charts import HeatMap, bins, output_file, show #Allows you to craete heatmaps
from bokeh.models import Rect

import numpy as np
import pdb


app = Flask(__name__)

#Function to only add points to a list if it is there or not
#Keeps order unlike sets
def point_to_index(points, point):
    try: 
        return points.index(point)
    except ValueError:
        points.append(point)
        return len(points)-1


#Create a set of names for what to be plotted by the interface
choices = ['Cross Section (5)',
           'Galaxy Bias, All (10)',
           '- Galaxy Bias, Low-z (5)',
           '- Galaxy Bias, Hi-z (5)',
           'Delta Source Photo-z And Sigma, All(10+1)',
           '- Delta Source Photo-z And Sigma, Low-z (5+1)',
           '- Delta Source Photo-z And Sigma, Hi-z (5+1)',
           'Delta Lens Photo-z and Sigma, All(10+1)',
           '- Delta Lens Photo-z And Sigma, Low-z (5+1)',
           '- Delta Lens Photo-z And Sigma, Hi-z (5+1)',
           'Shear Calibration, All (10)',
           '- Shear Calibration, Low-z (5)',
           '- Shear Calibratino, Hi-z (5)',
           'Intrinsic Alignment (4)',
           'All (46) (Not for the faint of heart/machine)']

@app.route('/', methods=['GET', 'POST'])
def home():
    parameters = np.loadtxt('/Users/penafiel/JPL/nora_data/parameter_samples_nuisance')
    #index values
    #index = parameters[:,0]
    #Galaxy bias parameters
    b_1 = parameters[:1000,8]
    b_2 = parameters[:1000,9]
    b_3 = parameters[:1000,10]
    b_4 = parameters[:1000,11]
    b_5 = parameters[:1000:,12]
    b_6 = parameters[:1000,13]
    b_7 = parameters[:1000,14]
    b_8 = parameters[:1000,15]
    b_9 = parameters[:1000,16]
    b_10 = parameters[:1000,17]

    #Source photo-z parameters
    delta_source_1 = parameters[:1000,18]
    delta_source_2 = parameters[:1000,19]
    delta_source_3 = parameters[:1000,20]
    delta_source_4 = parameters[:1000,21]
    delta_source_5 = parameters[:1000,22]
    delta_source_6 = parameters[:1000,23]
    delta_source_7 = parameters[:1000,24]
    delta_source_8 = parameters[:1000,25]
    delta_source_9 = parameters[:1000,26]
    delta_source_10 = parameters[:1000,27]
    sigma_souce = parameters[:1000,28]

    #Lens photo-z parameters
    delta_lens_1 = parameters[:1000,29]
    delta_lens_2 = parameters[:1000,30]
    delta_lens_3 = parameters[:1000,31]
    delta_lens_4 = parameters[:1000,32]
    delta_lens_5 = parameters[:1000,33]
    delta_lens_6 = parameters[:1000,34]
    delta_lens_7 = parameters[:1000,35]
    delta_lens_8 = parameters[:1000,36]
    delta_lens_9 = parameters[:1000,37]
    delta_lens_10 = parameters[:1000,38]
    sigma_lens = parameters[:1000,39]

    #Shear calibration parameters
    m_1 = parameters[:1000,40]
    m_2 = parameters[:1000,41]
    m_3 = parameters[:1000,42]
    m_4 = parameters[:1000,43]
    m_5 = parameters[:1000,44]
    m_6 = parameters[:1000,45]
    m_7 = parameters[:1000,46]
    m_8 = parameters[:1000,47]
    m_9 = parameters[:1000,48]
    m_10 = parameters[:1000,49]

    #Intrinsic alignment parameters
    A_0 = parameters[:1000,50]
    beta = parameters[:1000,51]
    eta = parameters[:1000,52]
    eta_z = parameters[:1000,53]

    index = np.arange(1000) + 1

    #Load the fiducial data vector
    fiduc_data = np.loadtxt('/Users/penafiel/JPL/nora_data/FIDUCIAL_WFIRST_all_2pt_WFIRST_LSST_test')
    fiduc_i = fiduc_data[:825,0]
    fiduc_c_l = fiduc_data[:825,1]

    #Load the inverse of the covariance matrix 
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

    #This array should be (1,2,3, ...., 1,2,3....)
    y = range(len(cov_825)) * len(cov_825)

    #This array should be (0,0,0,0,....,1,1,1,....)
    x = []
    for i in range(len(cov_825)):
        x = np.append(x, [i] * len(cov_825))
    aux_diag = (x==y)
    cov_825_col = np.reshape(cov_825, 825*825)
    cov_825_diag = cov_825_col[aux_diag]
    z_val = []

    if request.method == 'POST':
        if request.form['submit'] == 'submit_z':
            print request.form.getlist('z bin 1')
            #Get the z values not chosen
            #We will then make these 0
            for i in (np.arange(10) + 1):
                
                if request.form.get('z bin %i' %i) == None:
                    z_val = np.append(z_val, i)
            print z_val

            #script = ""
            #div = ""
            #return render_template('homepage.html', script=script, div=div)

    #Default behavior for z_val
    if z_val == []:
        z_val = [-1]
    #Arrays of z_value cross correlation
    z_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 10]
    z_j = [1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,3,4,5,6,7,8,9,10,4,5,6,7,8,9,10,5,6,7,8,9,10,6,7,8,9,10,7,8,9,10,8,9,10,9,10,10]

    #Make them length 825
    z_i = z_i*15
    z_j = z_j*15

    print z_i
    #Makes the z values not chosen = 0
    for i in z_val:
        aux_i = np.asarray(z_i) == int(i)
        aux_j = np.asarray(z_j) == int(i)

        for j in range(len(cov_825_diag)):
            if (aux_i[j] == True) | (aux_j[j]== True):
                cov_825[j] = 0
                cov_825[:,j] = 0
                cov_825[j,j] = 1#Make diagonals 1
        print cov_825 
    #Get the inverse covariance matrix
    inv_cov_825 = np.linalg.inv(cov_825)
    #Get the l values
    l = [2.404149e+01,3.473953e+01,5.019803e+01,7.253529e+01,1.048122e+02,1.514519e+02,2.188454e+02,   3.162278e+02,4.569437e+02,6.602757e+02,9.540871e+02,1.378639e+03,1.992110e+03,2.878565e+03,4.159476e+03]
    l_full = l*55
    l_full = np.asarray(l_full)
    aux_lin = (l_full < 600)
    aux_quasi = (l_full >=600) & (l_full < 1500)
    aux_non = (l_full >=1500)

    #Get delta_chi_2 for different l ranges
    #1 corresponds to linear scales
    #2 corresponds to quasilinear
    #3 corresponds to nonlinear

    #Redefine our covariance matrix
    cov_825_1_2 = np.copy(cov_825)
    cov_825_1_1 = np.copy(cov_825)
    cov_825_2_3 = np.copy(cov_825)
    cov_825_2_2 = np.copy(cov_825)
    cov_825_3_3 = np.copy(cov_825)

    aux_diag = (x==y)
    #Go through the matrix and 0 out the rows and columns for the
    #ranges that don't belong in the range, e.g. 1_2, we 0 out 3
    #Reshape the array and get the diagonals
    #Then calculate the inverse matrix
    #Do this for 1_2
    cov_825_1_2_col = np.reshape(cov_825_1_2, 825*825)
    cov_825_1_2_diag = cov_825_1_2_col[aux_diag]
    for i in range(len(cov_825_1_2_diag)):
        if aux_non[i] == True:
            cov_825_1_2[i] = 0
            cov_825_1_2[:,i] = 0
            cov_825_1_2[i,i] = 1 #make the diagonals 1
    inv_cov_825_1_2 = np.linalg.inv(cov_825_1_2)


    #Do this for 1_1
    cov_825_1_1_col = np.reshape(cov_825_1_1, 825*825)
    cov_825_1_1_diag = cov_825_1_2_col[aux_diag]
    for i in range(len(cov_825_1_1_diag)):
        if (aux_non[i] == True) | (aux_quasi[i] == True):
            cov_825_1_1[i] = 0
            cov_825_1_1[:,i] = 0
            cov_825_1_1[i,i] = 1 # make the diagonals 1
    inv_cov_825_1_1 = np.linalg.inv(cov_825_1_1)

    #Do this for 2_3
    cov_825_2_3_col = np.reshape(cov_825_2_3, 825*825)
    cov_825_2_3_diag = cov_825_2_3_col[aux_diag]
    for i in range(len(cov_825_2_3_diag)):
        if aux_lin[i] == True:
            cov_825_2_3[i] = 0
            cov_825_2_3[:,i] = 0
            cov_825_2_3[i,i] = 1 #Make the diagonals 1
    inv_cov_825_2_3 = np.linalg.inv(cov_825_2_3)

    #Do this for 2_2
    cov_825_2_2_col = np.reshape(cov_825_2_2, 825*825)
    cov_825_2_2_diag = cov_825_2_2_col[aux_diag]
    for i in range(len(cov_825_2_2_diag)):
        if (aux_lin[i] == True) | (aux_non[i] == True):
            cov_825_2_2[i] = 0
            cov_825_2_2[:,i] = 0
            cov_825_2_2[i,i] = 1 #Make the diagonals 1
    inv_cov_825_2_2 = np.linalg.inv(cov_825_2_2)

    #Do this for 3_3
    cov_825_3_3_col = np.reshape(cov_825_3_3, 825*825)
    cov_825_3_3_diag = cov_825_3_3_col[aux_diag]
    for i in range(len(cov_825_3_3_diag)):
        if (aux_lin[i] == True) | (aux_quasi[i] == True):
            cov_825_3_3[i] = 0
            cov_825_3_3[:,i] = 0
            cov_825_3_3[i,i] = 1 #Make the diagonals 1
    inv_cov_825_3_3 = np.linalg.inv(cov_825_3_3)

    delta_chi_2_1_3 = []
    delta_chi_2_1_2 = []
    delta_chi_2_1_1 = []
    delta_chi_2_2_3 = []
    delta_chi_2_2_2 = []
    delta_chi_2_3_3 = []
    
    delta_chi_2 = []
    #err_ggl = []
    #err_gc = []
    for i in index:
        #Load the data of your vector
        c_l_data = np.loadtxt('/Users/penafiel/JPL/nora_data/PCA_nuisance_only/WFIRST_all_2pt_%i' %i)
        c_l_i = c_l_data[:825,0]
        c_l = c_l_data[:825,1]
        
        #Get the difference matrix
        delta_mat = c_l - fiduc_c_l

        #Make sure this becomes 0 for the proper z values
        for i in z_val:
            aux_i = np.asarray(z_i) == int(i)
            aux_j = np.asarray(z_j) == int(i)
            delta_mat[aux_i | aux_j] = 0

        delta_1_2 = np.copy(delta_mat)
        delta_1_1 = np.copy(delta_mat)
        delta_2_3 = np.copy(delta_mat)
        delta_2_2 = np.copy(delta_mat)
        delta_3_3 = np.copy(delta_mat)

        #Make the ones corresponding to certain l s 0
        delta_1_2[aux_non] = 0
        delta_1_1[aux_quasi | aux_non] = 0
        delta_2_3[aux_lin] = 0
        delta_2_2[aux_lin | aux_non] = 0
        delta_3_3[aux_lin | aux_quasi] = 0
        
        #From here calculate the value of the delta chi_squared
        delta_chi_2_i = np.dot(np.dot(delta_mat, inv_cov_825), delta_mat)
        delta_chi_2_1_2_i = np.dot(np.dot(delta_1_2, inv_cov_825_1_2), delta_1_2)
        delta_chi_2_1_1_i = np.dot(np.dot(delta_1_1, inv_cov_825_1_1), delta_1_1)
        delta_chi_2_2_3_i = np.dot(np.dot(delta_2_3, inv_cov_825_2_3), delta_2_3)
        delta_chi_2_2_2_i = np.dot(np.dot(delta_2_2, inv_cov_825_2_2), delta_2_2)
        delta_chi_2_3_3_i = np.dot(np.dot(delta_3_3, inv_cov_825_3_3), delta_3_3)

        delta_chi_2 = np.append(delta_chi_2, delta_chi_2_i)
        delta_chi_2_1_2 = np.append(delta_chi_2_1_2, delta_chi_2_1_2_i)
        delta_chi_2_1_1 = np.append(delta_chi_2_1_1, delta_chi_2_1_1_i)
        delta_chi_2_2_3 = np.append(delta_chi_2_2_3, delta_chi_2_2_3_i)
        delta_chi_2_2_2 = np.append(delta_chi_2_2_2, delta_chi_2_2_2_i)
        delta_chi_2_3_3 = np.append(delta_chi_2_3_3, delta_chi_2_3_3_i)

        #We can separaete them into components
        #0 - 824 is weak lensing
        #825-1499 is galaxy galaxy lensing
        #1500 - 1649 is galaxy clustering
        #aux_wl = (c_l_i >= 0) and (c_l_i <= 824)
        #aux_ggl = (c_l_i >= 825) and (c_l_i <= 1499)
        #aux_gc = (c_l_i >= 1500) and (c_l_i <= 1649)

        #Get the absolute value of the difference
        #diff_wl_i = np.abs(fiduc_c_l[aux_wl] - c_l[aux_wl])
        #diff_ggl_i = np.abs(fiduc_c_l[aux_ggl] - c_l[aux_ggl])
        #diff_gc_i = np.abs(fiduc_c_l[aux_gc] - c_l[aux_gc])

        #Now we get the mean of these values and do a log base 10 
        #And that serves as our metric
        #err_wl = np.append(err_wl, np.log10(np.mean(diff_wl_i)))
        #err_ggl = np.append(err_ggl, np.log10(np.mean(diff_ggl_i)))
        #err_gc = np.append(err_gc, np.log10(np.mean(diff_gc_i)))

    #Now that we have the difference
    #We'll prototype this thing for now
    #We'll generate two dictionaries just in case
    data_dict = {'delta_chi_2':delta_chi_2,
                 'b_1':b_1,
                 'b_2':b_2,
                 'b_3':b_3,
                 'b_4':b_4,
                 'b_5':b_5,
                 'b_6':b_6,
                 'b_7':b_7,
                 'b_8':b_8,
                 'b_9':b_9,
                 'b_10':b_10,
                 'delta_source_1':delta_source_1,
                 'delta_source_2':delta_source_2,
                 'delta_source_3':delta_source_3,
                 'delta_source_4':delta_source_4,
                 'delta_source_5':delta_source_5,
                 'delta_source_6':delta_source_6,
                 'delta_source_7':delta_source_7,
                 'delta_source_8':delta_source_8,
                 'delta_source_9':delta_source_9,
                 'delta_source_10':delta_source_10,
                 'sigma_source':sigma_souce,
                 'delta_lens_1':delta_lens_1,
                 'delta_lens_2':delta_lens_2,
                 'delta_lens_3':delta_lens_3,
                 'delta_lens_4':delta_lens_4,
                 'delta_lens_5':delta_lens_5,
                 'delta_lens_6':delta_lens_6,
                 'delta_lens_7':delta_lens_7,
                 'delta_lens_8':delta_lens_8,
                 'delta_lens_9':delta_lens_9,
                 'delta_lens_10':delta_lens_10,
                 'sigma_lens': sigma_lens,
                 'm_1':m_1,
                 'm_2':m_2,
                 'm_3':m_3,
                 'm_4':m_4,
                 'm_5':m_5,
                 'm_6':m_6,
                 'm_7':m_7,
                 'm_8':m_8,
                 'm_9':m_9,
                 'm_10':m_10,
                 'A_0':A_0,
                 'beta':beta,
                 'eta':eta,
                 'eta_z':eta_z,
                 'index':index}

    #err_dict = {'err_wl': err_wl,
    #            'err_ggl':err_ggl,
    #            'err_gc':err_gc}
    delta_dict = {'delta_chi_2_1_3': delta_chi_2,
                  'delta_chi_2_1_2': delta_chi_2_1_2,
                  'delta_chi_2_1_1': delta_chi_2_1_1,
                  'delta_chi_2_2_3': delta_chi_2_2_3,
                  'delta_chi_2_2_2': delta_chi_2_2_2,
                  'delta_chi_2_3_3': delta_chi_2_3_3}
    
    #Create a source file in the context of bokeh
    source_data = ColumnDataSource(data=data_dict)

    #Also do it for the other data arrays
    source_delta = ColumnDataSource(data=delta_dict)

    #As a test, let's do 10 plots, since that's what my code previous code has
    #This one goes from some really light yellow thing (fff5ee) to Red
    colors = ['#fff5ee', '#ffe4e1', '#ffc1c1', '#eeb4b4', '#f08080', '#ee6363', '#d44942', '#cd0000', '#ff0000']
    #Mapper corresponding to the tot_tot_data
    mapper = LinearColorMapper(palette=colors, low=0, high=750)

    #Get the values chosen for x and y
    #Determine the selected choices
    current_x = request.args.get('choice_for_x')
    current_y = request.args.get('choice_for_y')

    if current_x == None:
        current_x = '- Galaxy Bias, Low-z (5)'

    if current_y == None:
        current_y = '- Galaxy Bias, Low-z (5)'

    #Fill our parameter values that we will be plotting
    if current_x == 'Cross Section (5)':
        param = ['b_1', 'delta_source_1', 'delta_lens_1', 'm_1', 'A_0']
        lengths = [0.08001, 0.0002001, 0.00005001, 0.0005001, 0.4001]
    #    widths = set([
    if current_x == 'Galaxy Bias, All (10)':
        param = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10']
        lengths = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005, 0.08006, 0.08007, 0.08008, 0.08009, 0.0801] 
    if current_x == '- Galaxy Bias, Low-z (5)':
        param = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5']
        lengths = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005]
    if current_x == '- Galaxy Bias, Hi-z (5)':
        param = ['b_6', 'b_7', 'b_8', 'b_9', 'b_10']
        lengths = [0.08006, 0.08007, 0.08008, 0.08009, 0.0801] 
    if current_x == 'Delta Source Photo-z And Sigma, All(10+1)':
        param = ['delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'delta_source_6', 'delta_source_7', 'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source']
        lengths = [0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004]
    if current_x == '- Delta Source Photo-z And Sigma, Low-z (5+1)':
        param = ['delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'sigma_source']
        lengths = [0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0004]
    if current_x == '- Delta Source Photo-z And Sigma, Hi-z (5+1)':
        param = [ 'delta_source_6', 'delta_source_7', 'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source']
        lengths = [0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004]
    if current_x == 'Delta Lens Photo-z and Sigma, All(10+1)':
        param = ['delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5', 'delta_lens_6', 'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens']
        lengths = [0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008]
    if current_x ==  '- Delta Lens Photo-z And Sigma, Low-z (5+1)':
        param = ['delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5','sigma_lens']
        lengths = [0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00008]
    if current_x == '- Delta Source Photo-z And Sigma, Hi-z (5+1)':
        param = ['delta_lens_6', 'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens']
        lengths = [ 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008]
    if current_x == 'Shear Calibration, All (10)':
        param = ['m_1', 'm_2', 'm_3', 'm_4', 'm_5', 'm_6', 'm_7', 'm_8', 'm_9', 'm_10']
        lengths = [0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005, 0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501]
    if current_x == '- Shear Calibration, Low-z (5)':
        param = ['m_1', 'm_2', 'm_3', 'm_4', 'm_5']
        lengths = [0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005]
    if current_x == '- Shear Calibration, Hi-z (5)':
        param = ['m_6', 'm_7', 'm_8', 'm_9', 'm_10']
        lengths = [0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501]
    if current_x == 'Intrinsic Alignment (4)':
        param = ['A_0', 'beta', 'eta', 'eta_z']
        lengths = [0.4001, 0.0752, 0.135, 0.3] 
    if current_x == 'All (46) (Not for the faint of heart/machine)':
        param = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10','delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'delta_source_6', 'delta_source_7',
                     'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source','delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5', 'delta_lens_6', 
                     'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens', 'm_1', 'm_2', 'm_3', 'm_4', 'm_5', 'm_6', 'm_7', 'm_8', 'm_9', 'm_10','A_0', 'beta', 'eta', 'eta_z']
        lengths = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005, 0.08006, 0.08007, 0.08008, 0.08009, 0.0801,
                       0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004,
                       0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008,
                       0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005, 0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501,
                       0.4001, 0.0752, 0.136, 0.3]
    if current_y == 'Cross Section (5)':
        param_y = ['b_1', 'delta_source_1', 'delta_lens_1', 'm_1', 'A_0']
        lengths_y = [0.08001, 0.0002001, 0.00005001, 0.0005001, 0.4001]
    if current_y == 'Galaxy Bias, All (10)':
        param_y = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10']
        lengths_y = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005, 0.08006, 0.08007, 0.08008, 0.08009, 0.0801] 
    if current_y == '- Galaxy Bias, Low-z (5)':
        param_y = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5']
        lengths_y = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005]
    if current_y == '- Galaxy Bias, Hi-z (5)':
        param_y = ['b_6', 'b_7', 'b_8', 'b_9', 'b_10']
        lengths_y = [0.08006, 0.08007, 0.08008, 0.08009, 0.0801] 
    if current_y == 'Delta Source Photo-z And Sigma, All(10+1)':
        param_y = ['delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'delta_source_6', 'delta_source_7', 'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source']
        lengths_y = [0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004]
    if current_y == '- Delta Source Photo-z And Sigma, Low-z (5+1)':
        param_y = ['delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'sigma_source']
        lengths_y = [0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0004]
    if current_y == '- Delta Source Photo-z And Sigma, Hi-z (5+1)':
        param_y = [ 'delta_source_6', 'delta_source_7', 'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source']
        lengths_y = [0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004]
    if current_y == 'Delta Lens Photo-z and Sigma, All(10+1)':
        param_y = ['delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5', 'delta_lens_6', 'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens']
        lengths_y = [0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008]
    if current_y ==  '- Delta Lens Photo-z And Sigma, Low-z (5+1)':
        param_y = ['delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5','sigma_lens']
        lengths_y = [0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00008]
    if current_y == '- Delta Source Photo-z And Sigma, Hi-z (5+1)':
        param_y = ['delta_lens_6', 'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens']
        lengths_y = [ 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008]
    if current_y == 'Shear Calibration, All (10)':
        param_y = ['m_1', 'm_2', 'm_3', 'm_4', 'm_5', 'm_6', 'm_7', 'm_8', 'm_9', 'm_10']
        lengths_y = [0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005, 0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501]
    if current_y == '- Shear Calibration, Low-z (5)':
        param_y = ['m_1', 'm_2', 'm_3', 'm_4', 'm_5']
        lengths_y = [0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005]
    if current_y == '- Shear Calibration, Hi-z (5)':
        param_y = ['m_6', 'm_7', 'm_8', 'm_9', 'm_10']
        lengths_y = [0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501]
    if current_y == 'Intrinsic Alignment (4)':
        param_y = ['A_0', 'beta', 'eta', 'eta_z']
        lengths_y = [0.4001, 0.0752, 0.135, 0.3] 
    if current_y == 'All (46) (Not for the faint of heart/machine)':
        param_y = ['b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10','delta_source_1', 'delta_source_2', 'delta_source_3', 'delta_source_4', 'delta_source_5', 'delta_source_6', 'delta_source_7',
                     'delta_source_8', 'delta_source_9', 'delta_source_10', 'sigma_source','delta_lens_1', 'delta_lens_2', 'delta_lens_3', 'delta_lens_4', 'delta_lens_5', 'delta_lens_6', 
                     'delta_lens_7', 'delta_lens_8', 'delta_lens_9', 'delta_lens_10', 'sigma_lens', 'm_1', 'm_2', 'm_3', 'm_4', 'm_5', 'm_6', 'm_7', 'm_8', 'm_9', 'm_10','A_0', 'beta', 'eta', 'eta_z']
        lengths_y = [0.08001, 0.08002, 0.08003, 0.08004, 0.08005, 0.08006, 0.08007, 0.08008, 0.08009, 0.0801,
                       0.0002001, 0.0002002, 0.0002003, 0.0002004, 0.0002005, 0.0002006, 0.0002007, 0.0002008, 0.0002009, 0.000201, 0.0004,
                       0.00005001, 0.00005002, 0.00005003, 0.00005004, 0.00005005, 0.00005006, 0.00005007, 0.00005008, 0.00005009, 0.0000501, 0.00008,
                       0.0005001, 0.0005002, 0.0005003, 0.0005004, 0.0005005, 0.0005006, 0.0005007, 0.0005008, 0.0005009, 0.000501,
                       0.4001, 0.0752, 0.135, 0.3]

    for i, j in zip(param_y, lengths_y):
        point_to_index(param, i)
        point_to_index(lengths, j)
    #Append the 2 parameter sets
    #param |= param_y
    #param_list = list(param)
    #lengths |= lengths_y
    #lengths_list = list(lengths)

    print 'param_list', param
    print 'lengths_lust', lengths
    #Get the number of figures
    fig_num = 0
    i = len(param) - 1
    while i > 0:
        fig_num += i
        i -= 1
    hover_list = []
    fig_list = []

    #Create a program that creates a corner plot
    #Can't use corner so had to make it myself
    #WHYYYYYYYYYYY

    max_par = len(param) - 1 #Get the max number of plots on the x/y axis

    for i in range(max_par):
        y_ax = np.copy(i)+1
        x_ax = 0
        while (x_ax < y_ax) & (x_ax < max_par):
            hover_i = HoverTool(tooltips=[('index', '$index'),
                                      ('(x,y,)', '($x, $y)'),
                                      ('Delta_chi_2', '@delta_chi_2')])
            fig_i = figure(plot_width=300, plot_height=300,tools=[hover_i, TapTool()])

            fig_i.grid.grid_line_color = None
            #Plots the rectangles
            fig_i_rect = fig_i.rect(param[x_ax], param[y_ax],width=lengths[x_ax], height=lengths[y_ax], alpha=0.8, source=source_data,fill_color={'field':'delta_chi_2', 'transform':mapper}, line_color=None)
            fig_i.yaxis.axis_label = param[y_ax]
            fig_i.xaxis.axis_label = param[x_ax]

            #Create glyphs for the highlighting portion, so that when it is tapped
            #the colors don't change

            selected = Rect(fill_color={'field':'delta_chi_2', 'transform':mapper}, fill_alpha=0.8, line_color=None)
            nonselected = Rect(fill_color={'field':'delta_chi_2', 'transform':mapper}, fill_alpha=0.8, line_color=None)

            fig_i_rect.selection_glyph = selected
            fig_i_rect.nonselection_glyph = nonselected
            fig_list.append(fig_i)
            x_ax += 1

    #Generate our matrix for the grid plot
    none_val = np.full(max_par * max_par, None) #Fill all of these values with None, much easier than making values become None
                                                #No error thi way
    fig_matrix = none_val.reshape(max_par, max_par) #Reshape it to a matrix
    fig_lower = np.tril_indices(max_par) #Call the lower triangular matrix
    fig_matrix[fig_lower] = fig_list #Make those values the ones in our figure

    fig_matrix = fig_matrix.tolist()
    #Creates the gridplot to be reminscient of a corner plot
    plot = gridplot(fig_matrix)


    #Creates the color bar and adds it to the right side of the big plot
    
    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size='12pt',
                    ticker=BasicTicker(desired_num_ticks=len(colors)),
                    label_standoff=6, border_line_color=None, location=(0,0))

    s_color = figure()
    #Since this basically creates another plot, we want to remove it
    #That's what the next couple of lines does
    s_color.grid.grid_line_color = None
    s_color.axis.axis_line_color = None
    s_color.add_layout(color_bar, 'left')
    s_color.toolbar.logo = None
    s_color.toolbar_location = None


    #Code to be utilized for the checkboxes in the interface
    #Note this code is in JavaScript
    code_l = """

    //All of our data
    var delta_chi_2 = source_data.data;
    var delta_chi_sub = source_delta.data;

    var sum = 0;
    // Get the value from our checkbox group
    var l_check = l_checkbox.active;

    // Create an if statement to make sure it doesn't give
    // an error if no checkbox is chosen
    if (l_check != []) {
        var l_start = String(l_checkbox.active[0] + 1);
        var l_end = String(l_checkbox.active[l_check.length - 1] + 1);
    
        //Change the value being plotted depending on the values
        // of the checkbox
        var delta_string = 'delta_chi_2_' + l_start + '_' + l_end;
        delta_chi_2['delta_chi_2'] = delta_chi_sub[delta_string];
    } //l_bin

    //if (l_check == []) {
    //   for (var i = 0; i<delta_chi_2['delta_chi_2'].length; i++) {
    //        delta_chi_2['delta_chi_2'][i] = sum;
    //    }//delta_chi_2
    //} //l_bin
    console.log(delta_chi_2['delta_chi_2']);
    source_data.trigger('change');
    """

    callback_l = CustomJS(args=dict(source_data=source_data, source_delta=source_delta), code=code_l)

    #Create a checkbox for the l ranges
    l_checkbox = CheckboxGroup(labels=['Linear Scales, l < 600', 'Quasi-Linear Scales, 600 <= l < 1500', 'Nonlinear Scales, l >= 1500'], active=[0,1,2], callback=callback_l)
    
    callback_l.args['l_checkbox'] = l_checkbox
    

    l = layout([[WidgetBox(l_checkbox)],[plot,s_color]])
    script, div = components(l)
    return render_template('homepage.html', script=script, div=div,
           choice_x=choices, current_x=current_x,
           choice_y=choices, current_y=current_y)

def create_z(z_range):
    #Get the z_values
    z_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 10]
    z_j = [1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,3,4,5,6,7,8,9,10,4,5,6,7,8,9,10,5,6,7,8,9,10,6,7,8,9,10,7,8,9,10,8,9,10,9,10,10]



#With debug=True, Flask Render will auto-reload when there are code changes
if __name__ == '__main__':
    #set debug to False in a production environment
    app.run(port=5001, debug=True)





