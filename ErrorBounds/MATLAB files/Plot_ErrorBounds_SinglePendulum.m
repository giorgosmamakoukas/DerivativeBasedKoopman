nn = n:-1:n-1;
errors_n = nan(tFinal/ts, Nstates);
Mvalues = Mvalues';
% Mvalues = [51.8451; 51.8865];
t_k = 1;
for t = ts : ts : tFinal
    for states = 1 : Nstates
%         for i = 1 : length(n)
            M = Mvalues(states,1);
            errors_n(t_k, states) = M * new_function_error(nn(states),t);
%         end
    end
    t_k = t_k+1;    
end

errors_n = [zeros(1,Nstates); errors_n];

labels = {'\theta (rad)', 'd\theta (rad/s)'};
addpath('C:\Users\mamak\Dropbox\Important Documents\Research\My Papers\Journal Papers\tFinalRO (2019)')
figure(1)
for states = 1 : Nstates
    subplot(1,Nstates,states)
    plot( linspace(0, tFinal, tFinal/ts+1),errors_n(:,states), 'g'); hold on;
    ylabel(labels(states));
    xlabel('T (s)')
    TROplots;
end
% l = legend('nn = 0', 'nn = 1', 'nn = 2');
% l.FontSize = 8;

% Analytical derivatives
Mvalues_analytical = [14.81];
errors_n_analytical = nan(tFinal/ts, Nstates);

t_k = 1;
for t = ts : ts : tFinal
    for states = 1 : Nstates
            M = Mvalues_analytical(1);
            errors_n_analytical(t_k, states) = M * new_function_error(nn(states),t);
%         end
    end
    t_k = t_k+1;    
end
errors_n_analytical = [zeros(1,Nstates); errors_n_analytical];

figure(1)
for states = 1 : Nstates
    subplot(1,Nstates,states)
    plot( linspace(0, tFinal, tFinal/ts+1),errors_n_analytical(:,states), 'r'); hold on;
    ylabel(labels(states));
    xlabel('T (s)')
    TROplots;
end



% Plot measured/actual errors
load('Data - Errors.mat');
E = reshape(max(error_history,[], Nstates), tFinal/ts+1,Nstates);
figure(1)
for states = 1 : Nstates
    subplot(1,Nstates,states)
    plot( linspace(0, tFinal, tFinal/ts+1),E(:,states), 'b.'); hold on;
    ylabel(labels(states));
    xlabel('T (s)')
    TROplots;
    box on;

end

h= legend('Data bound ', 'Taylor bound', 'Actual Error')
set(h,'FontSize',10);
% suptitle(' n = 1 ')
function[total_error] = new_function_error(n,t)

total_error = (t)^(n+1)/factorial(n+1);

end

