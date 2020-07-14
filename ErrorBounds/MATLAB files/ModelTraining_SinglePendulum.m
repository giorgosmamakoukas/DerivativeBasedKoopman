
addpath('C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\Sagital Glider');
Samples = Number_of_Samples;    
A = zeros(NKoopman, NKoopman);
G = zeros(length(A),length(A)); % 7 is the number of elements in the basis functions

rng(1)
Ps0_list = nan(length(Samples),NKoopman);
Psi_list = nan(length(Samples),NKoopman);

 for i = 1:1: Samples
    clc
    i
    pause(0.005)

    % States
    th0 = 4*pi*rand(1)-2*pi;
    dth0 = 10*rand(1)-5;
    s0 = [th0; dth0];
    u10 = 10 * rand(1) -5; 

    [~, sn] = ode45(@(t, s) sys_dynam(t,s, u10), [0, ts], s0);  
    sn = sn(end,:)';

    Ps0_list(i,:) = Psi_x(s0, u10)';
    Psi_list(i,:) = Psi_x(sn, u10)';
       
    [Atemp, Gtemp] = Koop_K(s0,sn,u10);
    A = A+Atemp;
    G = G+Gtemp;
end

save('Data - Observation Measurements', 'Ps0_list', 'Psi_list');
Kd = A * pinv(G);
K = real(logm(Kd)/ts);

save('Data - Learned Koopman Operators', 'K', 'Kd' );