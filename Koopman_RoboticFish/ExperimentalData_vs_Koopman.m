%% load all data from Maria's experiments
clear
close all;
clc 

ts = 0.005;
% Measure accuracy of trained Koopman
load('InterpolatedData.mat');
% positions = Lengths;

clearvars -except positions Lengths ts
load('TailFish_LQR.mat');
load('Data.mat');

addpath('C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences')
% addpath('C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\Sagital Glider')
global K
% tailfish_parameters();
% 
% Not sure these are helping
K(end-1:end,:) = 0; %% Do not predict control change (keep control constant)
% K(1:3,1:3) = 0; % dynamics don't depend on x and y for tail fish

Kd = expm(K*ts);
% positions = positions - 1;
% positions(end) = [];
% error = NaN(length(Ps0_list), 6);
% for i = 1 : 6 % scans all basis functions
%     for j = 1 : length(Ps0_list) % scans all trials
%         if intersect(positions,j) == j % skip crossover between trials
%             j = j+1;
%         end
%         error(j,i) = Psi_list(j,i) - Kd(i,:) * Ps0_list(j,:)';
%     end
% end
% %remove NaN between switching trials
% for i = length(positions) : -1 : 1
%     error(positions(i), :) = 0;
% end
%% Koopman simulated data
tic 

clear K_Psi
positions = Lengths; % this is the last indeces for a given trial --- I should predict the next at the -1 spot. 
K_Psi = NaN(positions(end), length(Kd)); % each column for a separate state
positions(end) = [];
K_Psi(1,:) = Ps0_list(1,:);
% K_Psi(j+1, :) = Kd * K_Psi(j, :)';
for j = 1 : length(Ps0_list)
    if intersect(positions,j) == j % skip crossover between trials
%         K_Psi(j+1,:) = Ps0_list(j,:); % last one, then switch
        j = j+1;
        K_Psi(j,:) = Ps0_list(j,:);
    end
    K_Psi(j+1, :) = Kd * K_Psi(j, :)';
end

% save('KoopmanTrajectory_matchingExperiments', 'K_Psi')
toc

plot(K_Psi(:,1));
%%
%% load all data from Maria's experiments
tic 

workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;

% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox\images\imdemos');
% Ask user to confirm or change.
% topLevelFolder = uigetdir(start_path);
topLevelFolder = 'C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences';
if topLevelFolder == 0
	return;
end
% Get list of all subfolders
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

clear singleSubFolder remain start_path topLevelFolder allSubFolders

    a0_exp = [];
    aa_exp = [];
    t_exp = [];
    v1_exp = [];
    v2_exp = [];
    omega_exp = [];
    x_exp = [];
    y_exp = [];
    psi_exp = [];
for i = 2 : numberOfFolders
    name = listOfFolderNames{i};
    cd(name);

%     files = dir('*.txt');
    files = dir('STATES*.txt');
    for ii=1:length(files)
        eval(['load ' files(ii).name ' -ascii']);
    end
    % control parameters
    omega_a = 2*pi;
    numbers_in_folder_name = regexp(name,'-?\d+\.?\d*|-?\d*\.?\d+','match');
    a0 = str2num(numbers_in_folder_name{end});
    aa = str2num(numbers_in_folder_name{end-1});
    
   
    for trial = 1 : 2
        no = num2str(trial);
        filename = eval((['STATES', no]));
        t_exp = [t_exp;filename(:,1)];
        v1_exp = [v1_exp;filename(:,2)];
        v2_exp = [v2_exp;filename(:,3)];
        omega_exp = [omega_exp;filename(:,4)];
        x_exp = [x_exp;filename(:,6)/1];
        y_exp = [y_exp;filename(:,7)/1];
%         psi_exp = [psi_exp;filename(:,8)];
        psi_exp = [psi_exp;filename(:,5)];

%         omega_exp = [omega_exp;omega(:)];
%         x_exp = [x_exp;x(:)];
%         y_exp = [y_exp;y(:)];
%         psi_exp = [psi_exp;psi(:)];
            
        a0_exp = [a0_exp;a0];
        aa_exp = [aa_exp; aa];
    end
end
switching_inx_exp = find(diff(t_exp)<0); % last indx before new trial
switching_inx_exp = [0;switching_inx_exp];
% t_exp(switching_inx_exp+1) = 0; % so that all start from zero seconds
format short;
clc
toc
%%
clc 
close all
cd 'C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences';
addpath('C:\Users\mamak\Dropbox\Important Documents\Research\Simulations\Inverted_Cart_Pendulum\altmany-export_fig-76bd7fa')
load('InterpolatedData.mat');
FinalTimePlot = 5; % Don't plot over the entire trajectory, just over 5 seconds

positions = Lengths;
positions = [0; positions];
BlueDotSize = 22;
TextSize = 18;


%Animate result

% writerObj = VideoWriter('KoopmanFitness_200ms.mp4');
% writerObj.Quality = 95;
% writerObj.FrameRate = 2;
% open(writerObj);
figure('Renderer', 'painters', 'Position', [0 200 1400, 800])
for trials = 1 : length(Lengths)-1 %% 22 trials
    clf;
%     close all
    t = 0 : ts : (positions(trials+1)-positions(trials)-1)*ts;
    t(1)
    t(end)
    t(find(t>FinalTimePlot)) = []; % remove points after 5 seconds
    x = x_int_list(positions(trials)+1:positions(trials+1));            x(length(t)+1:end) = [];
    subplot(2,3,1)
    plot(t, x,'color', 'g','linewidth', 3);
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 1);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    texp = t_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1))-t_exp(switching_inx_exp(trials)+1);  texp(find(texp>FinalTimePlot)) = [];
    xexp = x_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); xexp(length(texp)+1:end) = [];
    hold on; scatter(texp,xexp, BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)'); ylabel('x (m)');
    set(gca,'FontSize',TextSize)

    subplot(2,3,2)
    hold on;
    y = y_int_list(positions(trials)+1:positions(trials+1));            y(length(t)+1:end) = [];
    plot(t, y,'color', 'g','linewidth', 3)
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 2);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    yexp = y_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); yexp(length(texp)+1:end) = [];
    hold on; scatter(texp,yexp, BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)'); ylabel('y (m)');set(gca,'FontSize',TextSize)
    box on;
%         hold on; scatter(t,y,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('y');
    % 
    subplot(2,3,3)
    psi = psi_int_list(positions(trials)+1:positions(trials+1));            psi(length(t)+1:end) = [];
    plot(t, psi,'color', 'g','linewidth', 3)
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 3);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    psiexp = psi_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); psiexp(length(texp)+1:end) = [];
    hold on; scatter(texp,psiexp,BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)');ylabel('\psi (rad)'); set(gca,'FontSize',TextSize)
%         hold on; scatter(t,psi,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('psi');
    % 
    subplot(2,3,4)
    v1 = v1_int_list(positions(trials)+1:positions(trials+1));            v1(length(t)+1:end) = [];
    plot(t, v1,'color', 'g','linewidth', 3)
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 4);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    v1exp = v1_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); v1exp(length(texp)+1:end) = [];
    hold on; scatter(texp,v1exp, BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)'); ylabel('v_x (m/s)'); set(gca,'FontSize',TextSize)
%         hold on; scatter(t,v1,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('v1');

    subplot(2,3,5)
    v2 = v2_int_list(positions(trials)+1:positions(trials+1));            v2(length(t)+1:end) = [];
    plot(t, v2,'color', 'g','linewidth', 3)
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 5);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    v2exp = v2_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); v2exp(length(texp)+1:end) = [];
    hold on; scatter(texp,v2exp, BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)'); ylabel('v_y (m/s)'); set(gca,'FontSize',TextSize)
%         hold on; scatter(t,v2,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('v2');

    subplot(2,3,6)
    omega = omega_int_list(positions(trials)+1:positions(trials+1));            omega(length(t)+1:end) = [];
    plot(t, omega,'color', 'g','linewidth', 3)
    Ksim = K_Psi(positions(trials)+1:positions(trials+1), 6);           Ksim(length(t)+1:end) = [];
    hold on; plot(t, Ksim, 'color', 'r', 'linewidth', 3);
    omegaexp = omega_exp(switching_inx_exp(trials)+1:switching_inx_exp(trials+1)); omegaexp(length(texp)+1:end) = [];
    hold on; scatter(texp,omegaexp,BlueDotSize, 'filled', 'MarkerFaceColor', 'b'); 
    xlabel('t (sec)'); ylabel('\omega (rad/s)'); set(gca,'FontSize',TextSize)
    
%     str = sprintf('Trial',trial)
    suptitle(['Trial ', num2str(1+(mod(trials,2)==0)), ': Bias = ', num2str(a0_exp(trials)), ', Amp = ', num2str(aa_exp(trials))])
    set(gcf, 'color', 'w')
    pdfname = sprintf(['Koopman_vs_Experiment_2piAmp', num2str(aa_exp(trials)),'Bias',num2str(a0_exp(trials)), '_Trial ', num2str(1+(mod(trials,2)==0)),'.pdf']);
%     pdfname
%     export_fig(pdfname);
%     close all;
        pause(0.5);
%     writeVideo(writerObj,getframe(gcf));
%         pause(0.1);
end
%  close(writerObj);
