# This code uses experimental data to train an approximate Koopman operator for the tail-actuated robotic fish

######## 0. IMPORT PYTHON FUNCTIONS ########

from functions import * # import all functions in functions.py file
from numpy import zeros
from scipy import io, linalg

def TrainKoopman():

    ######## 1. IMPORT DATA ########

    mat = io.loadmat('InterpolatedData_200Hz.mat', squeeze_me=True)
    positions =  mat['Lengths'] - 1 # subtract 1 to convert MATLAB indices to python
    x = mat['x_int_list']
    y = mat['y_int_list']
    psi = mat['psi_int_list']
    v_x = mat['v1_int_list']
    v_y = mat['v2_int_list']
    omega = mat['omega_int_list']
    u1 = mat['u1_list']
    u2 = mat['u2_list']

    ######## 2. INITIALIZE A and G matrices
    A = zeros((62, 62)) # 62 is the size of the Ψ basis functions
    G = zeros((62, 62))

    ######## 3. TRAINING KOOPMAN ########
    for i in range(x.size-1):

        # print('{:.2f} % completed'.format(i/x.size*100))
        if i in positions:
            i += 1 # jump to next trial at the end of each trial

        # Create pair of state measurements
        s0 = [x[i], y[i], psi[i], v_x[i], v_y[i], omega[i]]
        sn = [x[i+1], y[i+1], psi[i+1], v_x[i+1], v_y[i+1], omega[i+1]]

        Atemp, Gtemp = A_and_G(s0,sn,[u1[i],u2[i]])
        A = A+Atemp;
        G = G+Gtemp;

    Koopman_d = dot(A,linalg.pinv2(G)) # more accurate than numpy
    # Koopman_d = dot(A,numpy.linalg.pinv(G))

    # io.savemat('SavedData.mat', {'A' : A, 'G': G, 'Kd': Koopman_d}) # save variables to Matlab file

    return Koopman_d
