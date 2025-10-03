% Change the parameter values
try
    size_vector_changes = 0;
    for j_ = 1 : size(param_new_names,1)
        if eval(['size(SET.param.' char(param_new_names(j_,:)) ',1)']) == 1
            size_vector_changes = size_vector_changes + eval(['size(SET.param.' char(param_new_names(j_,:)) ',1)']);
            eval(['SET.param.', char(param_new_names(j_,:)), ' = ', num2str(param_new_value(size_vector_changes,:)), ';']);
            warning(['Changed ' char(param_new_names(j_,:)) ' to ' num2str(param_new_value(size_vector_changes,:))])
        elseif eval(['size(SET.param.' char(param_new_names(j_,:)) ',1)']) > 1
            for k_ = 1 : eval(['size(SET.param.' char(param_new_names(j_,:)) ',1)'])
                eval(['SET.param.', char(param_new_names(j_,:)), '(' num2str(k_) ',1) = ', num2str(param_new_value(k_+size_vector_changes,1)), ';']);
                warning(['Changed ' char(param_new_names(j_,:)) ' to ' num2str(param_new_value(k_+size_vector_changes,1))])
            end
            size_vector_changes = size_vector_changes + eval(['size(SET.param.' char(param_new_names(j_,:)) ',1)']);
        end
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
[oo_.dr, info, M_, oo_] = resol(0, M_, options_,oo_);

% Convert Dynare matrices to structural form
dyn_in.M_  = M_;
dyn_in.oo_ = oo_;
dyn_in.solve = 1;
dyn_in.options_ = options_;

% Calculate steady state for new parameter values
varlist_temp = erase(oo_.var_list, '_flex');
y_bar = [];
for r = 1:length(oo_.var_list)
    try y_bar(r,1) = M_.params(strmatch(['ss_',varlist_temp{r,1}],M_.param_names,'exact'));, catch y_bar(r,1) = 1;end
end

% format longg;
% [[1:length(out.mats.y_bar)]' out.mats.y_bar out_new.mats.y_bar]
% clear varlist_temp;
% pause
%
%
% % Calculate new Q_adj and G_adj matrices (including the constant)
% Q    = out_new.mats.Q(1:n_bp,1:n_bp);
% C    = (eye(n_bp) - Q) * log(out_new.mats.y_bar);
%
% out_new.mats.C = C;
%
% out_new.mats.Q_adj = [out_new.mats.Q(1:n_bp,1:n_bp), C; zeros(1,n_bp),1];
% out_new.mats.G_adj = [out_new.mats.G];

dyn_in.linearize_around_diff_y = 1;
dyn_in.dyn_linearize_point = log(y_bar);

if bypass_Callum == 1
    G2=dyn_in.oo_.dr.ghu*sqrt(M_.Sigma_e);
    out_new.mats.G = G2(dyn_in.oo_.dr.inv_order_var,:);
    Dynare_Q = oo_.dr.ghx(dyn_in.oo_.dr.inv_order_var,:);
    for ii = 1:size(Dynare_Q,1)
        if ismember(ii,dyn_in.oo_.dr.state_var)
            out_new.mats.Q(:,ii) = Dynare_Q(:,find(dyn_in.oo_.dr.state_var==ii));
        else
            out_new.mats.Q(:,ii) = zeros(size(Dynare_Q,1),1);
        end
    end



else
    out_new = dyn_to_str(dyn_in);
end





out_new.mats.y_bar = y_bar;
