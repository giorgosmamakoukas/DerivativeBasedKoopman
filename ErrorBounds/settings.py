from numpy import zeros, sin, cos, asarray, dot, empty, amax, ones, concatenate, linspace, arange
from random import seed, uniform
from scipy import linalg, io
from math import pi
from scipy.integrate import odeint
from scipy.io import loadmat, savemat
# import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from scipy.special import factorial

Samples = 5000 # Number of random initial conditions for both training and testing data
ts = 0.02 # time spacing between training state measurements
tFinal = 0.3 # time horizon --- used in measuring error

NKoopman = 2 # Number of basis functions
Nstates = 2 # Number of system states
Ncontrol = 1 # Number of system inputs

timeSteps = round(tFinal/ts)+1

####### System parameters
g = 9.81 # gravitational constant
l = 1 # pendulum length

def single_pendulum(y,t, u):
    theta, dtheta = y
    f = asarray([dtheta, g/l * sin(theta) +  u ])
    return f

def Psi_k(s, u): # Evaluates basis functions Ψ(s(t_k))
    theta, dtheta = s
    psi = zeros([NKoopman,1])
    psi[0:2, 0] = s
    if NKoopman >= 3:
        psi[2, 0] = g/l * sin(theta) + u
    if NKoopman >= 4:
        psi[3, 0] = g/l * cos(theta) * dtheta
    return psi

def A_and_G(s_1, s_2, u): # Uses measurements s(t_k) & s(t_{k+1}) to calculate A and G
    A = dot(Psi_k(s_2, u), Psi_k(s_1, u).T)
    G = dot(Psi_k(s_1, u), Psi_k(s_1, u).T)
    return A, G

def new_function_error(n,t): # error bound of an arbitrary function f using its n derivatives at time t 
    total_error = t**(n+1)/factorial(n+1)
    return total_error
