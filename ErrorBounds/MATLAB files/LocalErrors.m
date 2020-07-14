%% load Psi(x_k) and Psi(x_{k+1}), and Kd
if isequal(set,'training')
    load('Data - Observation Measurements.mat')
    load('Data - Learned Koopman Operators.mat')
elseif isequal(set,'testing')
    load('Data - Koopman_TestingDataSet.mat');
end

%% Measure maximum local (across one time step) errors in Psi(x_{k+1}) - Kd * Psi(x_k)
local_errors = nan(Number_of_Samples, NKoopman);
for i = 1 : Number_of_Samples
    local_errors(i,:) = abs((Psi_list(i,:)' - Kd*Ps0_list(i,:)')');
end
max_local_errors = max(local_errors)';
load('Data - Errors.mat');  % to save error_history

% Calculate maximum actual error as a function of time for each of the
% states
max_real_error_progression =  max(error_history, [], 2);
max_real_error_progression = reshape(max_real_error_progression, [round(tFinal/ts)+1,Nstates]);
% for i = 1 : Nstates
% max_real_error_progression(i,:) = max(error_history(:,:,i)');
% end

% delete('*.mat');
max_global_errors = reshape(max(max(error_history)), [Nstates,1]);
save('Data - MaxErrors', 'max_local_errors', 'max_global_errors', 'max_real_error_progression');

fprintf('Max local errors\n')
fprintf(' theta: % .4f \n d\theta: % .4f \n', max_local_errors(1:Nstates)) 

% Range of training data
max_range = max(max(Psi_list), max(Ps0_list));
min_range = min(min(Psi_list), min(Ps0_list));
