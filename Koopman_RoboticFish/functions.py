# This file includes all user-defined functions

from math import atan, sqrt, sin, cos
from numpy import empty, sign, dot

def Psi_k(s,u): #Creates a vector of basis functions using states s and control u
    x, y, psi, v_x, v_y, omega = s # Store states in local variables

    if (v_y == 0) and (v_x == 0):
        atanvXvY = 0; # 0/0 gives NaN
        psi37 = 0;
        psi40 = 0;
        psi52 = 0;
        psi56 = 0;
    else:
        atanvXvY = atan(v_y/v_x);
        psi37 = v_x * pow(v_y,2) * omega/sqrt(pow(v_x,2)+pow(v_y,2));
        psi40 = pow(v_x,2) * v_y * omega / sqrt(pow(v_x,2) + pow(v_y,2)) * atanvXvY;
        psi52 = pow(v_x,2) * v_y * omega / sqrt(pow(v_x,2) + pow(v_y,2));
        psi56 = v_x * pow(v_y,2) * omega * atanvXvY / sqrt(pow(v_x,2) + pow(v_y,2));

    Psi = empty([62,1]); # declare memory to store psi vector

    # System States
    Psi[0,0] = x;
    Psi[1,0] = y;
    Psi[2,0] = psi;
    Psi[3,0] = v_x;
    Psi[4,0] = v_y;
    Psi[5,0] = omega;

    # f(t): terms that appear in dynamics
    Psi[6,0] = v_x * cos(psi) - v_y * sin(psi);
    Psi[7,0] = v_x * sin(psi) + v_y * cos(psi);
    Psi[8,0] =  v_y * omega;
    Psi[9,0] = pow(v_x,2);
    Psi[10,0] = pow(v_y,2);
    Psi[11,0] = v_x * omega;
    Psi[12,0] = v_x * v_y;
    Psi[13,0] = sign(omega) * pow(omega,2);

    # df(t)/dt: terms that appear in derivative of dynamics
    Psi[14,0] = v_y * omega * cos(psi);
    Psi[15,0] = pow(v_x,2) * cos(psi);
    Psi[16,0] = pow(v_y,2) * cos(psi);
    Psi[17,0] = v_x * omega * sin(psi);
    Psi[18,0] = v_x * v_y * sin(psi);

    Psi[19,0] = v_y * omega * sin(psi);
    Psi[20,0] = pow(v_x,2) * sin(psi);
    Psi[21,0] = pow(v_y,2) * sin(psi);
    Psi[22,0] = v_x * omega * cos(psi);
    Psi[23,0] = v_x * v_y * cos(psi);

    Psi[24,0] = v_x * pow(omega,2);
    Psi[25,0] = v_x * v_y * omega;
    Psi[26,0] = v_x * pow(v_y,2);
    Psi[27,0] = v_y * sign(omega) * pow(omega,2);
    Psi[28,0] = pow(v_x,3);

    Psi[29,0] = v_y * pow(omega,2);
    Psi[30,0] = v_x * omega * sqrt(pow(v_x,2) + pow(v_y,2));
    Psi[31,0] = v_y * omega * sqrt(pow(v_x,2) + pow(v_y,2)) * atanvXvY;
    Psi[32,0] = pow(v_x,2) * v_y;
    Psi[33,0] = v_x * sign(omega) * pow(omega,2);
    Psi[34,0] = pow(v_y,3);
    Psi[35,0] = pow(v_x,3) * atanvXvY;
    Psi[36,0] = v_x * pow(v_y,2) * atanvXvY;
    Psi[37,0] = psi37;
    Psi[38,0] = pow(v_x,2) * v_y * pow(atanvXvY,2);
    Psi[39,0] = pow(v_y,3) * pow(atanvXvY,2);
    Psi[40,0] = psi40;

    Psi[41,0] = pow(v_y,2) * omega;
    Psi[42,0] = v_x * v_y * sqrt(pow(v_x,2) + pow(v_y,2));
    Psi[43,0] = pow(v_y,2) * sqrt(pow(v_x,2) + pow(v_y,2)) * atanvXvY;
    Psi[44,0] = pow(v_x,2) * omega;
    Psi[45,0] = pow(v_x,2) * sqrt(pow(v_x,2) + pow(v_y,2)) * atanvXvY;
    Psi[46,0] = v_x * v_y * sign(omega) * omega;
    Psi[47,0] = pow(omega, 3);

    Psi[48,0] = v_y * omega * sqrt(pow(v_x,2) + pow(v_y,2));
    Psi[49,0] = pow(v_x,3);
    Psi[50,0] = v_x * pow(v_y,2);
    Psi[51,0] = pow(v_x,2) * v_y * atanvXvY;
    Psi[52,0] = psi52;

    Psi[53,0] = v_x * omega * sqrt(pow(v_x,2) + pow(v_y,2)) * atanvXvY;
    Psi[54,0] = pow(v_x,3) * pow(atanvXvY,2);
    Psi[55,0] = v_x * pow(v_y,2) * pow(atanvXvY,2);
    Psi[56,0] = psi56;
    Psi[57,0] = pow(v_y, 3) * atanvXvY;

    Psi[58,0] = v_x * pow(omega,2);
    Psi[59,0] = v_y * sign(omega) * pow(omega,2);

    # add control inputs
    Psi[60,0] = u[0];
    Psi[61,0] = u[1];

    return Psi

def A_and_G(s_1, s_2, u): # Uses measurements s(t_k) & s(t_{k+1}) to calculate A and G 
    A = dot(Psi_k(s_2, u), Psi_k(s_1, u).transpose());
    G = dot(Psi_k(s_1, u), Psi_k(s_1, u).transpose());
    return A, G
