% Set parameters for the model 

%% Initialise parameters

% *** Economic model parameters ***

% Aggregate parameters
param.eta       = 0.9;             % Elasticity of substitution of goods in consumption and investment
param.eta_g     = 0.9;             % Elasticity of substitution in government consumption bundle
param.zeta      = 0.95;            % Elasticity of substitution between capital and labour
param.vphi      = 0.25;            % Elasticity of substitution between value-added and intermediates
param.psi       = 0.4;             % Elasticity of substitution between intermediates
param.xi        = 2;               % Elasticity of substitution between labour types
param.bet       = 0.995;           % Discount rate
param.h         = 0.7;             % Habits
param.capup     = 0.5;             % Capital utilisation cost parameter
param.eps_w     = 5;               % Substitution between labour types within sectors
param.nu        = 2;               % Frisch
param.delt      = 0.02;            % Depreciation rate
param.Spp       = 3;               % Investment adjustment cost
param.chi_w     = 0.2;             % Backward looking component of wages PC 
param.rho_r     = 0.5;             % Lagged rates in Taylor rule
param.phi_pie   = 1.5;             % Inflation in Taylor rule
param.phi_pie_core = 0;            % Core inflation
param.phi_y     = 0;               % Output in Taylor rule
param.phi_gap   = 0.125;           % Output gap in Taylor rule
param.sig_r     = 1;               % Stdev of MP shock
param.rho_g     = 0.5;             % Persistence of government expenditure shock
param.sig_g     = 1;               % Standard deviation of public expenditure shock
param.n_ss      = 1/3;             % Steady state hours worked
param.A_N       = 7;               % Constant weight on disutility of labour
param.thet_w    = 0.75;            % Adjustment cost in wage PC
param.G_ss      = 1;               % Parameter to control the steady state of government spending
param.rk_wedge  = 1.8;             % Wedge on return on capital
param.cshift_ss = 1;               % Steady state of consumption shifter
param.rho_cshift= 0.5;             % Persistence of consumption shifter
param.sig_cshift= 1;               % Standard deviation of shock to consumption shifter
param.omega_r   = 0.99;             % Share of Ricardian consumers
param.rho_trans = 0.99;            % Persistence of shocks to transfers
param.sig_trans = 1;               % Standard deviation of shock to transfers
param.rho_tfp     = 0.5;           % Persistence of aggregate TFP shock
param.sig_tfp     = 1;             % Standard deviation of aggregate TFP shock
param.sig_mm      = 1;             % Standard deviation of aggregate correlation shock
param.n_ait = n_ait_parameter;      % AIT window (this is only relevant if phi_pie_bar > 0)
param.phi_pie_bar = 0;              % AIT parameter (default = 0 --> IT)
param.phi_pie_NYFed = 0;            % NY Fed AIT parameter (default = 0 --> IT)
param.phi_pie_forecast = 0;         % Targeting inflation forecast
param.phi_pie_headmin = 0;           % Targeting 50% inflation and 50% inflaiton in the mining sector

% ======= Sectoral parameters ==========

nsectors = SET.nsectors;

% *** Consumption weights ***
% Load matrix of sectoral consumption shares
for k_ = 1 : nsectors
     param.omega_c_(k_,1) = io.consumption_weights(k_);
     param.rho_om_c_(k_,1)   = 0.5;
     param.sig_om_c_(k_,1)   = 1;     
end

% *** Investment weights ***
for k_ = 1 : nsectors
    param.omega_i_(k_,1) = io.investment_weights(k_);
end

% *** Government expenditure ***
for k_ = 1 : nsectors
    param.omega_g_(k_,1) = io.government_weights(k_);
end

% *** Labour weights ***
% NB: These are the weights of labour in the utility function not in the production function
for k_ = 1 : nsectors
    param.omega_ns_(k_,1) = 1;
end

% Industry intermediate input shares
input_output_weights;
for k_ = 1 : nsectors
    for m_ = 1 : nsectors
    param.ii_weights_(m_,k_) = ii_weights(m_, k_);
    end
end

% Industry labour shares
param.lva_weights = lva_weights';
for k_ = 1 : nsectors
    param.omega_nd_(k_,1) = lva_weights(k_);
end

% Industry value added shares
param.va_weights = va_weights';
for k_ = 1 : nsectors
    param.omega_y_(k_,1) = va_weights(k_);
end

% Industry average tax rates
for k_ = 1 : nsectors
    param.tau_c_(k_,1) = io.taxes_ishare(k_);
end

% Industry productivity steady-states
for k_ = 1 : nsectors
    param.a_ss_(k_,1)  = 1;
    param.rho_a_(k_,1) = 0.5;
    param.sig_a_(k_,1) = 1;
end

% Markup steady-states
for k_ = 1 : nsectors
    param.mu_ss_(k_,1)  = 1;
    param.rho_mu_(k_,1) = 0.5;
    param.sig_mu_(k_,1) = 1;
end

% Markup correlations across sectors
for k_ = 1 : nsectors
    param.c_mu_mm_(k_,1) = 0;
end
% Mining markup correlated with manufacturing markup
% param.c_mu_mm_(find(strcmp(io.industry_names,'Mining')),1) = 1;
param.c_mu_mm_(find(strcmp(io.industry_names,'Manufacturing')),1) = -0.05;

% TFP-Markup correlations across sectors
for k_ = 1 : nsectors
    param.c_tfp_mm_(k_,1) = 0;
end

% Mining markup correlated with manufacturing markup
param.c_tfp_mm_(find(strcmp(io.industry_names,'Mining')),1) = 1;
param.c_tfp_mm_(find(strcmp(io.industry_names,'Manufacturing')),1) = 0;

% Industry Phillips curve stead-states
param.thet_  = 0.8 + zeros(nsectors,1); 
param.chi_p_ = 0.2 + zeros(nsectors,1); 

% Flexible price sectors
param.thet_(idx_flex,1) = 0.25;

% Semi-sticky price sectors
param.thet_(idx_semi_flex,1) = 0.65;

% Shares of GDP
param.consumption_share = io.consumption_share;
param.investment_share  = io.investment_share; 
param.government_share  = io.government_share;
param.employment_share  = empl.employment_share;

SET.param       = param;
SET.param_init  = param;       % Save initial set of parameter values