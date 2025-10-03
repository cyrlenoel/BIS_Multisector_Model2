% Calculates a temporary shock
cd(folder_path);

% i = 1 defines the baseline exerecise
i = 1;

% Estimate the model
run 'model_run.m';

y_bar = out.mats.y_bar;

% Store the endogenous variable names-indices
for j_  =1 : M_.endo_nbr
    eval(['SET.variable.', char(M_.endo_names(j_,:)), ' = ', num2str(j_), ';']);
end

for jj = 1 : numshk
    eval(['targetnum1(' num2str(jj) ',1)= strmatch(targetvar1(' num2str(jj) '),M_.endo_names,''exact'');']);
end

for kk = 1 : length(exovars1)
    eval(['exonum1(' num2str(kk) ',1) = strmatch(exovars1(' num2str(kk) '),M_.exo_names,''exact'');']);
end


% This vector determines:
% The size of the change in the target variable -- i.e. 1. stands for doubling the steady state value
% The number of periods for which the shock is active -- i.e. 4 stands for 4 periods
if i == 1
    % eval(['targetpath1    = +' num2str(target_change_size) '.*ones(1,' num2str(number_of_periods) ')']);%- log(y_bar(targetnum1));        % Desired change in commodity prices (in deviation from steady state)
    for tp_ = 1 : length(targetvar1)
        eval(strcat("targetpath1(",num2str(tp_),",1:",num2str(number_of_periods),")    = +", num2str(target_change_size(tp_)), ".*ones(1,", num2str(number_of_periods), ");"));%- log(y_bar(targetnum1));        % Desired change in commodity prices (in deviation from steady state)
    end
else
    load(['Simulations\output\' iso2 '\targetpath1_' extractBefore(extractAfter(selected_shocks{1},'eps_'),'_') '.mat'])
end

% Having set the scenario, now simulate it
scenario_path = scenario_multiple(oo_, M_, out, horizon, exonum1, M_.exo_nbr, ...
    targetpath1, targetnum1, bypass_Callum);


