%% model_run_auto.m
% This script loads the inputs needed and runs the model 

% Load io and employment matrices
[io, empl] = f_main_prepare_data(folder_path,iso2,source);
% folder_path, d_path, iso2, source are set in the `run_*_shock.m' script

% Add Dynare path again (In f_main_prepare_data.m, we restore the MATLAB default path 
% to avoid conflicts with commands that might mistakenly use functions or files from the Dynare path.
if ischar(d_path)
    addpath(d_path);
elseif iscell(d_path)
    for p = 1:length(d_path)
        addpath(d_path{p});
    end
end

% Check if the relevant folder exist
if ~exist(strcat(folder_path,'\output'), 'dir')
    mkdir(strcat(folder_path,'\output')); % Create the folder
end

% Settings about the sectors
SET.nsectors    = io.nsectors ; % Number of sectors in the economy

flexible_price_sectors      = {'Agriculture'; 'Mining'} ; 
idx_flex                    = find(contains(io.industry_names,flexible_price_sectors)) ;
semi_flexible_price_sectors = { 'Manufacturing', 'Utilities', 'Transport', 'Retail', 'Wholesale','WholesaleAndRetail'} ; 
idx_semi_flex               = find(contains(io.industry_names,semi_flexible_price_sectors)) ;

% Set the names of various other files
parameter_file = {@model_2_parameters} ;

% File to call parameters
parameter_file{1}() ;

if bypass_Callum == 0
    param.rho_a_(:,1) = 0;
end

SET.param       = param ;
SET.param_init  = param ;       % Save initial set of parameter values

% Used for simulations
% try
%     SET.param.vphi      = vphi_new_value ;      
%     warning(['Changed vphi to ' num2str(vphi_new_value)])
% catch
%     warning('Nothing changed');
% end

disp('Calculating initial steady state')

% Calculate initial steady state
ss_params   = calc_ss_initial(SET) ;
SET.x0      = ss_params ;    % Update steady state parameter guess

% Based on those results, update parameters also targeting income and VA shares
ss_params = calc_ss_initial_extended(SET) ;

SET.x0    = log(ss_params(1 : end-nsectors*2-(nsectors-1))) ;
SET.x0    = [SET.x0; zeros(nsectors-1,1)] ;
SET.x0    = SET.x0(2:end) ;

save output\steady_state_results ss_params ;

ss = update_steady_state_initial(SET,ss_params) ;   % Fill out the steady state

% Reset parameters with new steady state values
param.A_N   = ss.A_N ;                      % Constant in labour supply equation
param.G_ss  = ss.g ;                        % Government spending

% Weights in value-added and gross output bundles
param.omega_nd_ = ss.omega_nd_ ; 
param.omega_y_  = ss.omega_y_ ;
param.omega_ns_ = ss.omega_ns_ ;
param.a_ss_     = ss.a_ ;

% Update SET
SET.param = param ;
SET.ss    = ss ;

%==========================================================================
% From here call dynare

SET.param.n_ait = n_ait_parameter;

%dynare input_parameter_test.mod noclearall
eval(['dynare bis_multisector_model.mod noclearall -Dn_ait_loop=' num2str(SET.param.n_ait) ' -Dnsector=' num2str(SET.nsectors)]); 

%table(M_.param_names, M_.params) 

disp('Dynare step finished')

% =========================================================================

disp('Converting matrices to structural form')

% Convert Dynare matrices to structural form
dyn_in.M_  = M_ ;
dyn_in.oo_ = oo_ ;
dyn_in.solve = 1 ;
dyn_in.options_ = options_ ;

% Extract the vector of steady-state-values and create constant matrix
save('output\test_ss.mat', '-regexp','^ss_')
test_ss = load('output\test_ss.mat');
y_bar = [] ;
for r = 1:length(oo_.var_list)
    try y_bar(r,1) = test_ss.(['ss_' strrep(oo_.var_list{r,1},'_flex','')]);, catch y_bar(r,1) = 1;, end
end

dyn_in.linearize_around_diff_y = 1 ;
dyn_in.dyn_linearize_point     = log(y_bar) ;

if bypass_Callum == 1
    G2=dyn_in.oo_.dr.ghu*sqrt(M_.Sigma_e) ;
    out.mats.G = G2(dyn_in.oo_.dr.inv_order_var,:);
    Dynare_Q = oo_.dr.ghx(dyn_in.oo_.dr.inv_order_var,:);
    for ii = 1:size(Dynare_Q,1)
        if ismember(ii,dyn_in.oo_.dr.state_var) 
            out.mats.Q(:,ii) = Dynare_Q(:,find(dyn_in.oo_.dr.state_var==ii));
        else 
            out.mats.Q(:,ii) = zeros(size(Dynare_Q,1),1);
        end
    end



else
    out             = dyn_to_str(dyn_in) ;
end







out.mats.y_bar  = y_bar ;

% The order of variables is in M_.endo_names

save output\test_dynare out ;


