% Parameters 
Number_of_Samples = 5000; % Number of random initial conditions for both training and testing data
ts = 0.01; % time spacing between state measurements
tFinal = 0.3; % time horizon --- used in measuring error

NKoopman = 2; % Number of basis functions, including control terms
Nstates = 2; % Number of system states
Ncontrol = 1; % Number of system inputs

m = tFinal/ts; % Number of prediction steps
n = 1; % Number of derivatives used in basis functions

