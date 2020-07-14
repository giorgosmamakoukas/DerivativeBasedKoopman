# This file tests the accuracy of trained Koopman

######## 0. IMPORT PYTHON FUNCTIONS ########

from functions import * # import all functions in functions.py file
from TrainingKoopman import * # import all functions in functions.py file
import matplotlib.pyplot as plt
from numpy import arange, insert, linspace

######## 1. IMPORT EXPERIMENTAL DATA ########

mat = io.loadmat('InterpolatedData_200Hz.mat', squeeze_me=True)
positions =  mat['Lengths'] - 1 # subtract 1 to convert MATLAB indices to python
# positions includes indices with last measurement of each experiment
x = mat['x_int_list']
y = mat['y_int_list']
psi = mat['psi_int_list']
v_x = mat['v1_int_list']
v_y = mat['v2_int_list']
omega = mat['omega_int_list']
u1 = mat['u1_list']
u2 = mat['u2_list']
positions = insert(positions, 0, -1) # insert -1 as index that precedes the 1st experiment

######## 2. PREDICT DATA USING TRAINED KOOPMAN ########

Kd = TrainKoopman() # Train Koopman
for exp_i in range(0, positions.size -2): # for each experiment
    indx = positions[exp_i]+1 # beginning index of each trial

    Psi_predicted = empty((positions[exp_i+1]-(indx), 62))
    s0 = [x[indx], y[indx], psi[indx], v_x[indx], v_y[indx], omega[indx]]
    Psi_predicted[0,:] = Psi_k(s0, [u1[indx], u2[indx]]).transpose() # Initialize with same initial conditions as experiment

    for j in range(0, positions[exp_i+1]-1-(indx)):
        Psi_predicted[j+1, :] = dot(Kd,Psi_predicted[j, :])

    ######## 3. PLOT EXPERIMENTAL VS PREDICTED DATA ########

    ylabels = ['x (m)', 'y (m)', 'ψ (rad)', r'$\mathregular{v_x (m/s)}$', r'$\mathregular{v_x (m/s)}$', 'ω (rad/s)']
    exp_data = [x, y, psi, v_x, v_y, omega]

    time = linspace(0, 1./200*(j+2), j+2) # create time vector
    fig = plt.figure()
    for states_i in range(6):
        plt.subplot('23'+str(states_i+1)) # 2 rows # 3 columns
        plt.plot(time, Psi_predicted[:, states_i])
        plt.plot(time, exp_data[states_i][indx:positions[exp_i+1]])
        plt.ylabel(ylabels[states_i])
    plt.gca().legend(('Predicted','Experimental'))

    Amp_values = [15, 20, 25, 30]
    Bias_values = [-20, -30, -40, -50, 0, 20, 30, 40, 50]
    titles = 'Amp: ' + str(Amp_values[(exp_i)//18]) + ' Bias: ' + str(Bias_values[(exp_i % 18) //2])
    fig.suptitle(titles)
    plt.show(block=False)
    plt.pause(2)
    plt.close()
