#!/usr/bin/env python
# -*- coding: utf-8 -*-
from settings import *

# Load learned Koopman operator
KoopmanOperators = loadmat('Data - Koopman_and_BasisFunctions.mat', squeeze_me=True)
Kd = KoopmanOperators['Kd']

# Initialize variable to store error in states over time horizon for each sample
error_history = empty((timeSteps, Samples, Nstates))

for random_i in range(0, Samples):

    # Generate random initial conditions for states and controls
    th0 = uniform(-2*pi, 2*pi)
    dth0 = uniform(-5, 5)
    u10 = uniform(-5, 5)

    s0 = [th0, dth0]
    Psi = Psi_k(s0, u10) # evaluate basis functions at t = 0

    s_Koopman = empty((timeSteps, Nstates)) # store state trajectories
    s_Koopman[0,:] = Psi[0:Nstates].T

    tt = linspace(0,tFinal, timeSteps)

    for i in range(1, timeSteps):
        Psi = dot(Kd,Psi) # Propagate with Koopman: Ψ_{k+1} = K_d * Ψ_k
        s_Koopman[i,:] = Psi[0:Nstates].T

    # Real dynamics evolution
    tt = linspace(0,tFinal, timeSteps)
    s_real = odeint(single_pendulum, s0, tt, args=(u10,))

    # Errors
    for i in range(Nstates):
        error_history[:,random_i, i] = abs(s_real[:,i] - s_Koopman[:,i])

error_history = amax(error_history, axis = 1)
savemat('Data - Errors.mat', {'error_history' : error_history}) # save variables to Matlab file
