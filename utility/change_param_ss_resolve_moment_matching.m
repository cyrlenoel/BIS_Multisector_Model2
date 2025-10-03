% Change the parameter values
try
    for j_ = 1 : size(param_new_value,1);

        eval(['SET.param.', char(param_new_names(j_,:)), ' = ', num2str(param_new_value(j_,:)), ';']);
        warning(['Changed ' char(param_new_names(j_,:)) ' to ' num2str(param_new_value(j_,:))])

    end
catch
    warning('Nothing changed');
end

% Updated the steady state
disp('Calculating initial steady state')

ss_params_update = calc_ss(SET);
ss_update = update_steady_state(SET,ss_params_update);
SET.ss = ss_update;

% Store the new parameter values
SET_temp = struct;
SET_fields = fieldnames(SET.param);
for i_ = 1: length(SET_fields)
    if size(eval(['SET.param.', char(SET_fields(i_))]),1) == 1
        eval(['SET_temp.param.', char(SET_fields(i_)), ' = SET.param.', char(SET_fields(i_)), ';']);
    elseif (size(eval(['SET.param.', char(SET_fields(i_))]),1) > 1) & (size(eval(['SET.param.', char(SET_fields(i_))]),2) == 1)
        for ii_ = 1:size(eval(['SET.param.', char(SET_fields(i_))]),1)
            eval(['SET_temp.param.', char(SET_fields(i_)), num2str(ii_), ' = SET.param.', char(SET_fields(i_)), '(', num2str(ii_),')',';']);
        end
    else
        for ii_ = 1:size(eval(['SET.param.', char(SET_fields(i_))]),1)
            for iii_ = 1:size(eval(['SET.param.', char(SET_fields(i_))]),2)
                eval(['SET_temp.param.', char(SET_fields(i_)), num2str(ii_), '_', num2str(iii_), ' = SET.param.', char(SET_fields(i_)), '(', num2str(ii_),',', num2str(iii_),')',';']);
            end
        end
    end
end

% Store the new steady state values
SET_fields = fieldnames(SET.ss);
for i_ = 1: length(SET_fields)
    if size(eval(['SET.ss.', char(SET_fields(i_))]),1) == 1
        eval(['SET_temp.param.ss_', char(SET_fields(i_)), ' = SET.ss.', char(SET_fields(i_)), ';']);
    elseif (size(eval(['SET.ss.', char(SET_fields(i_))]),1) > 1) & (size(eval(['SET.ss.', char(SET_fields(i_))]),2) == 1)
        for ii_ = 1:size(eval(['SET.ss.', char(SET_fields(i_))]),1)
            eval(['SET_temp.param.ss_', char(SET_fields(i_)), num2str(ii_), ' = SET.ss.', char(SET_fields(i_)), '(', num2str(ii_),')',';']);
        end
    else
        for ii_ = 1:size(eval(['SET.ss.', char(SET_fields(i_))]),1)
            for iii_ = 1:size(eval(['SET.ss.', char(SET_fields(i_))]),2)
                eval(['SET_temp.param.ss_', char(SET_fields(i_)), num2str(ii_), '_', num2str(iii_), ' = SET.ss.', char(SET_fields(i_)), '(', num2str(ii_),',', num2str(iii_),')',';']);
            end
        end
    end
end

% Change the values in M_
for j_ = 1 : size(M_.params,1);
    if isequal(char(M_.param_names(j_,:)), 'ss_n')
        M_.params(j_) = M_.params(j_);
    else
        eval(['M_.params(j_) = SET_temp.param.', char(M_.param_names(j_,:)), ';']);
    end

end

params_new = M_.params;
params_comparison = [params_orig, params_new];

% Resolve the model (resol is a Dynare function)
% [oo_.dr, info, M_, oo_] = resol(0, M_, options_,oo_);
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);