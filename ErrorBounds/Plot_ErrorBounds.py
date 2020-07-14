#!/usr/bin/env python
# -*- coding: utf-8 -*-
from settings import *

####### 1. Calculate error bounds
mat1 = loadmat('Data - Koopman_and_BasisFunctions.mat', squeeze_me=True)
max_local_errors = mat1['max_local_errors']

nn = arange(NKoopman-1, NKoopman-3, -1) # excludes NKoopman-3
coeff = ts**(nn+1) / factorial(nn+1)
MaxDer = max_local_errors[0:Nstates]/coeff # max magnitude of derivatives

E_bound = empty((timeSteps-1, Nstates))

# % Mvalues = [51.8451; 51.8865];
t = ts
for t_k in range(timeSteps-1):
    for states in range(Nstates):
        E_bound[t_k, states] = MaxDer[states] * new_function_error(nn[states],t)
    t += ts

E_bound= concatenate((zeros((1,Nstates)), E_bound), axis = 0) # add 0 error for t = 0

data = loadmat('Data - Errors.mat', squeeze_me=True)
E = data['error_history']

####### 2. Compare via plotting calculated error bounds with actual errors
time = arange(0, tFinal+ts, ts) # create time vector
ylabels = ['θ error (rad)', '$dθ$ error (rad/s)']
for states in range(Nstates):
    subplot('12'+str(states+1))
    plot(time, E_bound[:,states], 'g') # Verification of data bound from data
    plot(time, E[:,states], 'b.'); # actual error
    ylabel(ylabels[states])

tight_layout()
gca().legend(('Data-driven bound', 'Actual error'))
show(block=False)
pause(2)
close()
