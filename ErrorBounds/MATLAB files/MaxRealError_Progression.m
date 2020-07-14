
settings;
load('Data - MaxErrors.mat');

labels = {'x (m)', 'y (m)', '\psi (rad)', 'v_x (m/s)', 'v_y (m/s)', '\omega (rad/s)'};
addpath('C:\Users\mamak\Dropbox\Important Documents\Research\My Papers\Journal Papers\TRO (2019)')
figure()
for states = 1 : Nstates
    subplot(1,6,states); 
    for i = 1 : length(n)
        plot(linspace(0,tFinal,round(tFinal/ts)+1),max_real_error_progression(:,states)); hold on;
    end
    ylabel(labels(states));
    xlabel('T (s)')
    TROplots;
end
l = legend('n = 0', 'n = 1', 'n = 2');
l.FontSize = 8;
% suptitle('Error bound dependency on prediction horizon and number of derivatives')

function[total_error] = new_function_error(n,t)

total_error = (t)^(n+1)/factorial(n+1);

end