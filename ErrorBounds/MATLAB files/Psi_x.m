function psi = Psi_x(s, u)

theta = s(1); dtheta = s(2); 
g = 9.81; l = 1;
% States
Psi(1) = theta; 
Psi(2) = dtheta; 
% Psi(3) = g/l * sin(theta) + u;

Psi = Psi';
psi = Psi;
end