from settings import *

seed(1)

# Samples = Number_of_Samples
A = zeros((NKoopman, NKoopman))
G = zeros((NKoopman, NKoopman))

Ps0_list = empty((Samples,NKoopman))
Psi_list = empty((Samples,NKoopman))

for i in range (Samples):

    # Sample states
    th0 = uniform(-2*pi, 2*pi)
    dth0 = uniform(-5, 5)
    s0 = [th0, dth0]
    u10 = uniform(-5, 5)

    # Simulate system forward
    sn = odeint(single_pendulum, s0, [0, ts], args=(u10,))
    sn = sn[-1,:]

    # Evaluate basis functions at t = 0 and t = ts
    Ps0_list[i,:] = Psi_k(s0, u10).T
    Psi_list[i,:] = Psi_k(sn, u10).T

    [Atemp, Gtemp] = A_and_G(s0,sn,u10);
    A = A+Atemp
    G = G+Gtemp

Kd = dot(A,linalg.pinv2(G)) # more accurate than numpy
print(Kd)

## Measure maximum local (across one time step) errors in Ψ(s_{k+1}) - Kd*Ψ(s_k)
local_errors = empty([Samples, NKoopman])
for i in range(Samples):
    local_errors[i,:] = abs(Psi_list[i,:]- dot(Kd,Ps0_list[i,:]))
max_local_errors = amax(local_errors, axis = 0)
print('Max local errors in theta: %.5f and dtheta: %.5f ' % tuple(max_local_errors[0:2]))

# Save trained Koopman and basis functions measurements used to obtain it
io.savemat('Data - Koopman_and_BasisFunctions.mat', {'Kd' : Kd, 'max_local_errors' : max_local_errors}) # save variables to Matlab file
