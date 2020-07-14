# import numpy as np
# def pend(y, t, b, c):
#     # theta, omega = y
#     theta, dtheta = y
#     dydt = [dtheta, 9.81 * np.sin(theta) + b]
#     # dydt = [omega, -b*omega - c*np.sin(theta)]
#     return dydt
# # ...
# # We assume the constants are b = 0.25 and c = 5.0:
#
# b = 10
# c = 5.0
# # For initial conditions, we assume the pendulum is nearly vertical with theta(0) = pi - 0.1, and it initially at rest, so omega(0) = 0. Then the vector of initial conditions is
#
#
# y0 = [3.1 - 0.1, 0.0]
# #We generate a solution 101 evenly spaced samples in the interval 0 <= t <= 10. So our array of times is:
#
# t = np.linspace(0, 10, 101)
#
# from scipy.integrate import odeint
# sol = odeint(pend, y0, t, args=(b, c))
#
# import matplotlib.pyplot as plt
# plt.plot(t, sol[:, 0], 'b', label='theta(t)')
# plt.plot(t, sol[:, 1], 'g', label='omega(t)')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.grid()
# plt.show()
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt

def single_pendulum(y,t, u):
    g = 9.81
    l = 1
    theta, dtheta = y
    f = [dtheta, g/l * np.sin(theta) +  u ]
    print(u)
    return f
# def pend(y, t, b, c):
#     # theta, omega = y
#     theta, dtheta = y
#     dydt = [dtheta, 9.81 * np.sin(theta) + b]
#     # dydt = [omega, -b*omega - c*np.sin(theta)]
#     return dydt
# initial condition
y0 = [5, 5]
u = 0.1
# time points
t = np.linspace(0,20)

# solve ODE
sol = odeint(single_pendulum, y0, t, args=(u,))
print(sol.shape)

import matplotlib.pyplot as plt
plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
