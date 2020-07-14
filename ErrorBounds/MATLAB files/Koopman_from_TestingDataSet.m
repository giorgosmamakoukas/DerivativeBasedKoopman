% Train Koopman from testingDataSet
load('Data - Errors.mat');

A = zeros(NKoopman, NKoopman);
G = zeros(NKoopman, NKoopman); % 7 is the number of elements in the basis functions
Ps0_list = nan(length(TestingDataSet)*length(TestingDataSet(:,:,1))-1, NKoopman);
Psi_list = nan(length(TestingDataSet)*length(TestingDataSet(:,:,1))-1, NKoopman);
for samples_i = 1 : Number_of_Samples
    for samples_j = 1 : tFinal/ts
        s0 = TestingDataSet(samples_j,   1:Nstates, samples_i);
        sn = TestingDataSet(samples_j+1, 1:Nstates, samples_i);
        
        u10 = TestingDataSet(samples_j+1, Nstates+1, samples_i);
        
        [Atemp, Gtemp] = Koop_K(s0,sn,u10);
        
        Ps0_list(samples_j+(samples_i-1)*(length(TestingDataSet(:,:,1))-1),:) = Psi_x(s0, u10)';
        Psi_list(samples_j+(samples_i-1)*(length(TestingDataSet(:,:,1))-1),:) = Psi_x(sn, u10)';
        A = A+Atemp;
        G = G+Gtemp;
    end
end

Kd = A * pinv(G);

save('Data - Koopman_TestingDataSet', 'Kd', 'Ps0_list', 'Psi_list');
