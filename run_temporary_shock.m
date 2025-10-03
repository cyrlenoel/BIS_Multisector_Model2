%% run_temporary_shock.m
% This script runs the simulations for a temporary shock
close all; clear; clc; tic; 
fprintf('-------------------\n\n');
fprintf('Script execution has started. Estimated runtime: approximately 520 seconds.\n');
fprintf('Tested on: CPU - 13th Gen Intel(R) Core(TM) i7-1370P (1.90 GHz), RAM - 64.0 GB, OS - Windows 11 Enterprise.\n');
fprintf('Performance may vary depending on your system configuration.\n\n');
fprintf('-------------------\n');
pause(5);

% Folder settings
[folder_path, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename());
cd(folder_path); addpath(genpath(folder_path));
allPaths = strsplit(path, pathsep);
isDynarePath = false(1, numel(allPaths));
for f = 1:numel(allPaths)
    isDynarePath(f) = exist(fullfile(allPaths{f}, 'dynare.m'), 'file') == 2 && exist(fullfile(allPaths{f}, 'dynare_version.m'), 'file') == 2;
end
d_path = allPaths(isDynarePath);
if ~isempty(d_path) % Check if Dynare path is already added
    fprintf('-------------------\n\n');
    disp(d_path');
    fprintf("The path(s) listed above containing 'Dynare' have been found in the MATLAB search path. Proceeding to the next step...\n");
    pause(5);
else % if not, ask user
    d_path = input( ...
        "-------------------\n\n" + ...
        "\nThe model requires Dynare 5.5 to run. Please ensure it is installed before proceeding.\n" + ...
        "Enter the full Dynare path on your local drive (eg, C:\\Program Files\\Dynare\\5.5\\matlab) and press `Enter`.\n > ", 's');
    d_path = strrep(strrep(d_path, '''', ''), '"', '');
    if exist(d_path, "dir") && ischar(d_path)
        addpath(d_path); fprintf("The following path was successfully added to the MATLAB search path:\n%s\n\n", d_path);
    else
        error("Error: The specified Dynare path does not exist or is invalid. Check the path and try again. Execution has been stopped.");
    end
end

% Select the country and source
[iso2, source] = f_select_country_source(); % function saved under `input\progs'
% `iso2' and `source` must be entered by the user when prompted. Available options and instructions are provided in the prompts. 
source = ['_' source];

% Option to use the algorithm by Callum Jones
bypass_Callum = 1; % bypass the dyn_to_str algorithm by Callum Jones; not with permanent shocks

% Set the horizon of the IRFs and the number of shocks
horizon = 20;
% If numshk is larger than 1, selected_shock and selected_target must be of
% the same size as the number of shocks
numshk = 1;

% Select the shocks
selected_shocks = {'eps_mm'};
exovars1 = selected_shocks;

% Select the target
selected_targets = {'gam_2'};
targetvar1 = selected_targets;

% Specify the number of periods for which the shock is active
number_of_periods = 4;

% Specify the size of the change in the target variable 
% 1 stands for doubling the steady state value
target_change_size = repmat(.25,length(targetvar1),number_of_periods);

% Specify the number of quarters for AIT rule
% n_ait_parameter = 1 >>> IT
n_ait_parameter = 8;

run 'multiple_shocks.m'
% Here we back out the path of the shock to make sure that we apply the
% same shock in each of the following exercises
scenario_path_baseline = scenario_path;
SET_baseline = SET;

% Check if the relevant folder exist
if ~exist(strcat("output\temporary_shock\",iso2), 'dir')
    mkdir(strcat("output\temporary_shock\",iso2)); % Create the folder
end
save(['output\temporary_shock\' iso2 '\baseline_' extractAfter(selected_shocks{1},'eps_')],"scenario_path_baseline","SET_baseline");
targetpath1 = [];
targetnum1 = [];
if  contains(selected_shocks,'eps_mm')
    eval(['targetpath1 = scenario_path(SET.variable.eps_m,1:' num2str(number_of_periods) ');']);
    eval(['targetnum1 = SET.variable.eps_m;']);
    save(['output\temporary_shock\' iso2 '\targetpath1_eps_m' ],"targetpath1","targetnum1");
else
    for targp = 1:length(selected_shocks)
        eval(['targetpath1 = [targetpath1; scenario_path(SET.variable.' extractAfter(selected_shocks{targp},'eps_') ',1:' num2str(number_of_periods) ')];']); % Check that it becomes a matrix once num_of_periods>1
        eval(['targetnum1 = [targetnum1; SET.variable.' extractAfter(selected_shocks{targp},'eps_') '];']);
    end
    save(['output\temporary_shock\' iso2 '\targetpath1_' extractAfter(selected_shocks{1},'eps_')],"targetpath1","targetnum1");
end

% Save the original parameter values
params_orig = M_.params;

% Define the exercises
% Core
exercises = struct();
exercises.id_2 = struct();
exercises.id_2.id = 2; % Must be at least 2 because baseline exercise is 1
% Set the names of the parameters to change
exercises.id_2.names =        {'phi_pie';'phi_pie_core'}; 
exercises.id_2.names_legend = {'\phi_{\pi}';'\phi_{core}'}; % LaTeX format
% Set the new values of the parameters to change
exercises.id_2.values = [0;1.5];
% Set the target variable
exercises.id_2.targetvar1 = extractAfter(selected_shocks{1},'eps_');

% AIT (arithmetic average)
exercises.id_3 = struct();
exercises.id_3.id = 3;
exercises.id_3.names =        {'phi_pie';'phi_pie_bar'}; 
exercises.id_3.names_legend = {'\phi_{\pi}';'\phi_{\overline{\pi}}'}; % LaTeX format
exercises.id_3.values = [0;1.5];
exercises.id_3.targetvar1 = extractAfter(selected_shocks{1},'eps_');

exercises_fields = fieldnames(exercises);

% Change parameters
for exercise_ = 1:length(exercises_fields)
    SET = SET_baseline;
    i = str2double(extractAfter(char(exercises_fields(exercise_)),'id_'));
    
    % Input the target variable
    eval(['targetvar1 = exercises.', char(exercises_fields(exercise_)), '.targetvar1;']);

    % Input the new parameter values
    eval(['param_new_names = exercises.', char(exercises_fields(exercise_)), '.names;']);
    eval(['param_new_value = exercises.', char(exercises_fields(exercise_)), '.values;']);
    
    % Changes the parameter and steady-state values and resolves the model
    run 'change_param_ss_resolve.m'
    impose_target_path_sim(i, exercise_, exercises, exercises_fields, iso2);
    % Having set the scenario, now simulate it
    scenario_path = scenario_multiple(oo_, M_, out_new, horizon, exonum1, M_.exo_nbr, targetpath1, targetnum1, bypass_Callum);
    
    eval(['scenario_path_',  char(exercises_fields(exercise_)), '= scenario_path;']);
    eval(['SET_',  char(exercises_fields(exercise_)), '= SET;']);
    % eval(['save output\temporary_shock\',iso2,'\', char(exercises_fields(exercise_)),'_ait_',num2str(n_ait_parameter), ' scenario_path_',  char(exercises_fields(exercise_)), ' SET_',  char(exercises_fields(exercise_))]);
    eval(['save output\temporary_shock\',iso2,'\', char(exercises_fields(exercise_)), '_', extractAfter(selected_shocks{1},'eps_'), ' scenario_path_',  char(exercises_fields(exercise_)), ' SET_',  char(exercises_fields(exercise_))]);
end



%% Variable to plot IRFs for
vcell = {'pie','invest','r','y_va','c','rr'}; 
annualisation_cell = {'400','100','400','100','100','400'};
titlecell = {'Inflation','Investment','Policy rate','Real output',...
    'Real consumption','Real interest rate'};
obs_units = {'pp','pp','pp','%','%','pp'};
fontSize = 15;
fontName = 'Times New Roman';

%% Plot the IRFs comparing the different parameters assigned under each
% exercise and the same monetary policy rule

% Check if the relevant folder exists
if ~exist(strcat("graphs\temporary_shock\",iso2), 'dir')
    mkdir(strcat("graphs\temporary_shock\",iso2)); % Create the folder
end

% Import the results 
clearvars -except exercises exercises_fields source iso2 n_ait_parameter vcell annualisation_cell titlecell selected_shocks folder_path fontSize fontName obs_units
load(['output\temporary_shock\' iso2 '\baseline_' extractAfter(selected_shocks{1},'eps_') '.mat'])
for i = 1:length(exercises_fields)
    load(['output\temporary_shock\' iso2 '\' exercises_fields{i} '_' extractAfter(selected_shocks{1},'eps_') '.mat'])
end

% Figure
figure_name = strcat(folder_path,'\graphs\temporary_shock\',iso2,'\',iso2,'_temp_shock_',selected_shocks{1});
fprintf('Saving figures, please wait...\n %s.jpg\n %s.eps\n', figure_name, figure_name);
close all;
jj=0;figure;
set(gcf, 'Position', [100, 100, 1400, 1400]);  % [left, bottom, width, height]
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', fontSize);
set(gca, 'FontName', fontName);
for n = 1:numel(vcell)
    v = vcell{n};
    annualisation = annualisation_cell{n};
    jj = jj +1;
    subplot(2,3,jj)
    grid on; grid minor; box on; hold on;

    eval(['plot(scenario_path_baseline(SET_baseline.variable.', v,',:)*', annualisation,'); hold on;']);
    for i = 1:length(exercises_fields)
        eval(['plot(scenario_path_', exercises_fields{i}, '(SET_', exercises_fields{i} ,'.variable.', v,',:)*', annualisation,')']);
    end
    title([iso2  ' ' lower(titlecell{n})], 'FontSize', fontSize, 'FontName', fontName);
    yline(0, '-black', 'LineWidth', 1.5); 
    xline(0, '-black', 'LineWidth', 1.5);
    ylabel([obs_units{n}], 'FontSize', fontSize, 'FontName', fontName);
    xlabel('Quarters', 'FontSize', fontSize, 'FontName', fontName)
    ax = gca;
    ax.XAxis.FontSize = fontSize-2; 
    ax.YAxis.FontSize = fontSize-2; 
end
legend_handle = legend(['Headline inflation targeting'], ...
       ['Core inflation targeting'], ...
       ['Average inflation targeting'], ...
       'Location', 'northeast','FontSize',fontSize-3,'FontName', fontName,'Orientation','Horizontal', 'Interpreter', 'latex');
set(legend_handle,'Position',[0.39, 0.01, 0.25, 0.04],'Units','normalized')
print(figure_name,'-djpeg','-r600');
print(figure_name,'-depsc','-r300');
fprintf('Script execution completed successfully.\n');

toc
