function [ f ] = sys_dynam( t, s, u )

g = 9.81; l = 1;
theta = s(1); 
dtheta = s(2);

% System dynamics
f(1) = dtheta;
f(2) = g / l * sin(theta) + u;
f = f';

end