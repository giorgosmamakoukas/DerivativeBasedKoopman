% Load learned Koopman operator for Electric Fish dynamics
KoopmanOperators = load('Data - Learned Koopman Operators.mat');
Kd = real(KoopmanOperators.Kd);


states = [1,2]; labels = {'\theta (rad)', 'd\theta (rad/s)'};
settings;

Samples = Number_of_Samples;
errors_allRandomSamples = nan(Samples,length(states));
error_history = nan(round(tFinal/ts+1),Samples, Nstates);

% Keep range of state values to predict error
max_values = zeros(Nstates,1);
min_values = zeros(Nstates,1);
max_control_values = zeros(Ncontrol,1);
min_control_values = zeros(Ncontrol,1);

TestingDataSet = nan(tFinal/ts+1, Nstates+Ncontrol, Samples);
for random_i = 1 : 1: Samples
%     rng(random_i)
    clc
    random_i
    % Generate random initial conditions
            th0 = 4*pi*rand(1)-2*pi;
            dth0 = 5*rand(1)-5;
            s0 = [th0; dth0];
            u10 = 10 * rand(1) -5;

            Psi_0 = Psi_x(s0, u10);



        % Koopman prediction
        s = nan(round(tFinal/ts), Nstates);
        i = 0;
        for t = 0 : ts : tFinal
            i = i+1;
            temp = Kd^(t/ts) * Psi_0;
            s(i,:) = real(temp(1:Nstates)');
        end

    % Real dynamics evolution
    t = 0:ts:tFinal;
    [~, s_real] = ode45(@(t, s) sys_dynam(t,s, u10), t, s0);
    TestingDataSet(:,:, random_i) = [s_real, u10*ones(length(s_real),1)];

    temp_control_values = u10;
    max_control_values = max(max_control_values, temp_control_values);
    min_control_values = min(min_control_values, temp_control_values);

   max_values_temp = max(s_real);
   min_values_temp = min(s_real);

   for max_i = 1 : length(max_values_temp)
       if max_values(max_i) < max_values_temp(max_i)
           max_values(max_i) = max_values_temp(max_i);
       end
      if min_values(max_i) > min_values_temp(max_i)
           min_values(max_i) = min_values_temp(max_i);
       end
   end

    % Errors
    for i = 1 : Nstates
%         error(i) = norm(s_real(:,states(i)) - s(:,states(i)));
        error_history(:,random_i, i) = abs(s_real(:,states(i)) - s(:,states(i)));
%         average_error(random_i, i) = mean(abs(s_real(:,states(i)) - s(:,states(i))));

    end
%     errors_allRandomSamples(random_i,:) = error;
end


figure()
for i  = 1 : length(states)
%     subplot(5,3,i); stdshade(error_history(:,:, i), 0.9, 'b'); hold on;
subplot(1,2,i); plot(max(error_history(:,:,i)), 'b');
%     subplot(3,3,i); plot(mean(error_history(:,:,i))-mean(data.error_history(:,:,i)), 'k', 'linewidth', 2);
ylabel(labels(i));
% mean(mean(error_history(:,:,i))-mean(data.error_history(:,:,i)))

end
% fprintf('Errors \n x: %.3f \n y: %.3f \n z: %.3f \n vx : %.3f \n vy : %.3f \n vz : %.3f \n p : %.3f \n q : %.3f \n r: %.3f \n', error)
save('Data - Errors', 'error_history', 'TestingDataSet');
