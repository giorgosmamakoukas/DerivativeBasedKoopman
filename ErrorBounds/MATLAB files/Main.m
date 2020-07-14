clear; clc; close all;
settings;

% 1. Train Koopman operators K and Kd from simulation data
    rng(1);
    ModelTraining_SinglePendulum;
    % Saves [K, Kd] in "Data - Learned Koopman Operators"
    % Saves [Psi_k, Psi_{k+1}] in "Data - Observation Measurements"

clearvars -except m n Ncontrol NKoopman Nstates Number_of_Samples tFinal ts
clc;

% 2. Measure prediction error
    rng(1);
    KoopmanPredictionError_SinglePendulum;
    % Saves the error as a function of time (error_history) and the domain of the states from the true dynamics (TestingDataSet) in "Data - Errors"


clearvars -except m n Ncontrol NKoopman Nstates Number_of_Samples tFinal ts
clc;

% 3. Train Koopman and [Psi_k, Psi_{k+1}] from testing data set (Step 2.)
    Koopman_from_TestingDataSet;
    % Saves [Kd, Psi_k, Psi_{k+1}] in "Data - Koopman_TestingDataSet"

clearvars -except m n Ncontrol NKoopman Nstates Number_of_Samples tFinal ts
clc;

%%
% 4. Measure local errors in testing data set using measurements from Step 3.
    set = 'training'; % Use Kd and [Psi_k, Psi_{k+1}] from training dataset to measure error
%     set = 'testing'; % Use Kd and [Psi_k, Psi_{k+1}] from testing dataset to measure error
    LocalErrors;

clearvars -except m n Ncontrol NKoopman Nstates Number_of_Samples tFinal ts error_history max_local_errors max_global_errors

% 5. Calculate
    FunctionErrorFormula;
    % Saves max real local and global errors for each state in "'Data - MaxErrors'"

clearvars -except m n Ncontrol NKoopman Nstates Number_of_Samples tFinal ts error_history max_local_errors max_global_errors Mvalues

% 6. Plot error bounds as a function of time
Plot_ErrorBounds_SinglePendulum

% 7. Delete data points
% delete('*.mat');
