%% load all data from Maria's experiments
clear 
close all;
clc 


workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;

addpath 'C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences'
% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox\images\imdemos');
% Ask user to confirm or change.
% topLevelFolder = uigetdir(start_path);
topLevelFolder = 'C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences';
if topLevelFolder == 0
	return;
end
% Get list of all subfolders.
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

    x_int_list = [];
    y_int_list = [];
    psi_int_list = [];
    v1_int_list = [];
    v2_int_list = [];
    omega_int_list = [];
    omega_a_list = [];
    u1_list = [];
    u2_list = [];
    
    addpath('C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\Sagital Glider');
    tailfish_parameters();
    Lengths = [0]; % length of each data set

for i = 2 : numberOfFolders
    name = listOfFolderNames{i}
    cd(name)

%     files = dir('*.txt');
    files = dir('STATES*.txt');
    for j=1:length(files)
        eval(['load ' files(j).name ' -ascii']);
    end
    % control parameters
    omega_a = 2*pi;
    numbers_in_folder_name = regexp(name,'-?\d+\.?\d*|-?\d*\.?\d+','match');
    a0 = str2num(numbers_in_folder_name{end});
    
    aa = str2num(numbers_in_folder_name{end-1});
%     a0 = 40;  
    a0 = a0 * pi/180; % [-40, + 40] degrees convert to radians        
    aa = aa * pi/180;
%             close all;

    % In each subfolder, store both trials
    for trial = 1 : 3
        no = num2str(trial);
        currentFolder = cd;
        
        if isfile(['STATES', no, '.txt'])
        
            filename = eval((['STATES', no]));
            t = filename(:,1);
            v1 = filename(:,2);
            v2 = filename(:,3);
            omega = filename(:,4);
            x = filename(:,6);
            y = filename(:,7);
            psi = filename(:,5);
% 
            %% Calculate body-frame velocities 
            % The TVplusLPfilter filters data, but discards 2 first and
            % last values. 
%             clear v1 v2 omega
%             dt = 0.3333;
%             [x, vx_w] = TVplusLPfilter(x, dt);
%             [y, vy_w] = TVplusLPfilter(y, dt);
%             [psi, omega] = TVplusLPfilter(psi, dt);
%             t(end-1:end) = []; t(1:2) = []; 
            % Convert world-frame linear velocities to body-frame linear
            % velocities
%             
%             for jj = 1 : length(x)
%                 th = psi(jj);
%                 R = [cos(th), -sin(th); sin(th), cos(th)];
%                 temp = R\[vx_w(jj); vy_w(jj)];
%                 v1(jj) = temp(1);
%                 v2(jj) = temp(2);
%             end
            
            %%
            tF = t(end) - t(1); t = t - t(1); 
            NN = tF/0.005+1;
%             tint = 0: 0.33: tF; 
%             tint = 0: 0.005: tF; 

%             NN = length(tint);


            figure('Renderer', 'painters', 'Position', [0 400 1000 600])
            subplot(2,3,1)
            x_int=interp1(t, x,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), x_int,'color', 'g','linewidth', 2)
            hold on; scatter(t,x, 10, 'filled', 'MarkerFaceColor', 'b'); ylabel('x');

            % 
    %         figure(2)
            subplot(2,3,2)
            hold on;
            y_int=interp1(t, y,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), y_int,'color', 'g',  'linewidth', 2)
            hold on; scatter(t,y,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('y');
%             % 
            subplot(2,3,3)
    %         figure(3)
            psi_int=interp1(t, psi,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), psi_int,'color', 'g',  'linewidth', 2)
            hold on; scatter(t,psi,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('psi');
            % 
    %         figure(4)
            subplot(2,3,4)
            v1_int=interp1(t, v1,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), v1_int,'color', 'g',  'linewidth', 2)
            hold on; scatter(t,v1,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('v1');

    %         figure(5)
            subplot(2,3,5)
            v2_int=interp1(t, v2,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), v2_int,'color', 'g',  'linewidth', 2)
            hold on; scatter(t,v2,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('v2');

    %         figure(6)
            subplot(2,3,6)
            omega_int=interp1(t, omega,linspace(0,tF,NN), 'spline');
            plot(linspace(0,tF,NN), omega_int,'color', 'g',  'linewidth', 2)
            hold on; scatter(t,omega,  10, 'filled', 'MarkerFaceColor', 'b'); ylabel('omega');

            title_name = sprintf('Amp %d Bias %.1f, Trial: %d', aa*180/pi, a0*180/pi, trial);
            suptitle(title_name);

            u1 = aa^2 * (3 - 3/2 * a0^2 - 3/8 * aa^2);
            u2 = aa^2 * a0; 
            x_int_list = [x_int_list; x_int'];
            y_int_list = [y_int_list; y_int'];
            psi_int_list = [psi_int_list; psi_int'];
            v1_int_list = [v1_int_list; v1_int'];
            v2_int_list = [v2_int_list; v2_int'];
            omega_int_list = [omega_int_list; omega_int'];
            u1_list = [u1_list; u1 * ones(length(x_int),1)];
            u2_list = [u2_list; u2 * ones(length(x_int),1)];
            omega_a_list = [omega_a_list; omega_a * ones(length(x_int),1)];
            Lengths = [Lengths; Lengths(end) + length(x_int)];

            clear x_int y_int psi_int v1_int v2_int omega_int ASTATES1 ASTATES2
            
        end
        % Save
%         if trial == 1
%             save('InterpolatedData', 'x_int', 'y_int', 'psi_int', 'v1_int', 'v2_int', 'omega_int', 'omega_a', 'u1', 'u2');
%         else
%             save('InterpolatedData_v2', 'x_int', 'y_int', 'psi_int', 'v1_int', 'v2_int', 'omega_int', 'omega_a', 'u1', 'u2');
%         end
    end
%     close all;
end

Lengths(1) = [];
cd('C:\Users\mamak\Dropbox\Important Documents\Research\Project 3 - Koopman\TailFish\Fish Data - April 2019 (30 Hz  Frame Rate)\2pi - SmoothDataAndFiniteDifferences');
save('InterpolatedData', 'x_int_list', 'y_int_list', 'psi_int_list', 'v1_int_list', 'v2_int_list', 'omega_int_list', 'omega_a_list', 'u1_list', 'u2_list', 'Lengths');
format short