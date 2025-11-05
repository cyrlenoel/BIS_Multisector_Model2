SET.param = param ;

ss_params_update   = calc_ss(SET) ; 
ss_update          = update_steady_state(SET, ss_params_update) ;
SET.ss             = ss_update ;

% Store the new parameter values
SET_temp = struct ;
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
for j_ = 1 : size(M_.params,1) ;
    if isequal(char(M_.param_names(j_,:)), 'ss_n')
        M_.params(j_) = M_.params(j_) ;                   % Need to be fixed
    else
        eval(['M_.params(j_) = SET_temp.param.', char(M_.param_names(j_,:)), ';']) ;
    end

end

params_new = M_.params;
params_comparison = [params_orig, params_new]; 

% Resolve the model for new parameter values
[oo_.dr, info, M_, oo_] = resol(0, M_, options_,oo_) ;

disp('Converting matrices to structural form')
 
% Convert Dynare matrices to structural form
dyn_in.M_  = M_ ;
dyn_in.oo_ = oo_ ;
dyn_in.solve = 1 ;
dyn_in.options_ = options_ ;
 
update_all_steady_states_for_Dynare ;


% Extract the vector of steady-state-values and create constant matrix
save('output/test_ss.mat', '-regexp','^ss_')
test_ss = load('output/test_ss.mat');
y_bar = [] ;
for r = 1:length(oo_.var_list)
    try y_bar(r,1) = test_ss.(['ss_' strrep(oo_.var_list{r,1},'_flex','')]);, catch y_bar(r,1) = 1;, end
end
 
dyn_in.linearize_around_diff_y = 1 ;
dyn_in.dyn_linearize_point     = log(y_bar) ;
 
out             = dyn_to_str(dyn_in) ;
 
out.mats.y_bar  = y_bar ;

Q_new = out.mats.Q ;
G_new = out.mats.G ;
A0_new = out.mats.A0 ;
A1_new = out.mats.A1 ;
B0_new = out.mats.B0 ;
D0_new = out.mats.D0 ;

A_new = A0_new \ A1_new ; B_new = A0_new \ B0_new ;
D_new = A0_new \ D0_new ;