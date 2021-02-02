from settings import *

# load Koopman Model
KoopmanOperators = loadmat('Data - Koopman_and_BasisFunctions.mat', squeeze_me=True)
Kd = KoopmanOperators['Kd']

# define A, B, Q, R
a = Kd[0:3, 0:3]
b = Kd[0:3, 3, None]
q = np.diag([3, 1, 0.1])
r = np.diag([0.1])

# Calculate LQR gains
kLQR = dlqr(a,b,q,r)

# Sample states
th0 = uniform(-pi, pi)
dth0 = uniform(-2, 2)
s0 = [th0, dth0]

# sn = np.asarray(s0[:,:, None])
sn = np.asarray(s0)
sn = sn.reshape(2, 1)

print(s0)
# print(kLQR.shape)
# print(Psi_k(s0, 0)[0:3])
# print(dot(kLQR, Psi_k(s0, 0)[0:3]))
for i in np.arange(0.0, 10, ts):
    u10 = -dot(kLQR, Psi_k(s0, 0)[0:3])[0,0]
    stemp = odeint(single_pendulum, s0, [0, ts], args=(u10,))
    # s0 = stemp[-1,:, None]
    s0 = stemp[-1,:]
    sn = np.hstack((sn, stemp[-1,:, None]))

print(sn[:,-10:])
#
# # Simulate system forward
# sn = odeint(single_pendulum, s0, [0, ts], args=(u10,))
# sn = sn[-1,:]
#
# # Evaluate basis functions at t = 0 and t = ts
# Ps0_list[i,:] = Psi_k(s0, u10).T
# Psi_list[i,:] = Psi_k(sn, u10).T
#
# [Atemp, Gtemp] = A_and_G(s0,sn,u10);
# A = A+Atemp
# G = G+Gtemp
