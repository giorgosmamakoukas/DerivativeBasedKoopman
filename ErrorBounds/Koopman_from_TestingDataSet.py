from scipy.io import loadmat, savemat
from numpy import zeros, empty
from settings import *
from PredictionError import *
from scipy import linalg

mat = loadmat('Data - Errors.mat', squeeze_me=True)

# Train Koopman from testingDataSet
A = zeros((NKoopman, NKoopman))
G = zeros((NKoopman, NKoopman)) # 7 is the number of elements in the basis functions

Ps0_list = empty((Samples*timeSteps-1, NKoopman))
Psi_list = empty((Samples*timeSteps-1, NKoopman))

for samples_i in range(Samples):
    for samples_j in range(timeSteps-1):
        s0 = TestingDataSet[samples_j,   0:Nstates, samples_i] # NEED TO CHANGE --> 1:Nstates
        sn = TestingDataSet[samples_j+1, 0:Nstates, samples_i] # NEED TO CHANGE --> 1:Nstates

        u10 = TestingDataSet[samples_j+1, Nstates, samples_i]

        [Atemp, Gtemp] = A_and_G(s0,sn,u10)

        Ps0_list[samples_j+(samples_i-1)*timeSteps-1,:] = Psi_k(s0, u10).T
        Psi_list[samples_j+(samples_i-1)*timeSteps-1,:] = Psi_k(sn, u10).T
        A = A+Atemp
        G = G+Gtemp

Kd = dot(A,linalg.pinv2(G)) # more accurate than numpy

savemat('Data - Koopman_TestingDataSet.mat', {'Kd' : Kd, 'Ps0_list': Ps0_list, 'Psi_list': Psi_list}) # save variables to Matlab file
