param       = SET.param;
ss          = SET.ss;


% An individual parameter (this is a scalar)
d_eta                    = param.eta;
d_eta_g                  = param.eta_g;
d_zeta                   = param.zeta;
d_vphi                   = param.vphi;
d_psi                    = param.psi;
d_xi                     = param.xi;
d_bet                    = param.bet;
d_h                      = param.h;
d_capup                  = param.capup;
d_eps_w                  = param.eps_w;
d_nu                     = param.nu;
d_delt                   = param.delt;
d_Spp                    = param.Spp;
d_chi_w                  = param.chi_w;
d_rho_r                  = param.rho_r;
d_phi_pie                = param.phi_pie;
d_phi_pie_core           = param.phi_pie_core;
d_phi_pie_NYFed          = param.phi_pie_NYFed;
d_phi_pie_bar            = param.phi_pie_bar;
d_phi_pie_forecast       = param.phi_pie_forecast;
d_phi_pie_headmin        = param.phi_pie_headmin;
d_phi_y                  = param.phi_y;
d_phi_gap                = param.phi_gap;
d_sig_r                  = param.sig_r;
d_rho_g                  = param.rho_g;
d_sig_g                  = param.sig_g;
d_ss_n                   = param.n_ss;
d_A_N                    = param.A_N;
d_thet_w                 = param.thet_w;
d_ss_G                   = param.G_ss;
d_rk_wedge               = param.rk_wedge;
d_ss_cshift              = param.cshift_ss;
d_rho_cshift             = param.rho_cshift;
d_sig_cshift             = param.sig_cshift;
d_omega_r                = param.omega_r;
d_rho_trans              = param.rho_trans;
d_sig_trans              = param.sig_trans;
d_lva_weights            = param.lva_weights;
d_va_weights             = param.va_weights;
d_consumption_share      = param.consumption_share;
d_investment_share       = param.investment_share;
d_government_share       = param.government_share;
d_employment_share       = param.employment_share;
d_A_N                    = param.A_N;
d_ss_omega_y_            = param.omega_y_;
d_ss_omega_ns_           = param.omega_ns_;
d_n_ait                  = param.n_ait;
d_rho_tfp                = param.rho_tfp;
d_sig_tfp                = param.sig_tfp;
d_sig_mm                = param.sig_mm;

% A sectoral parameter (this is a vector)
d_omega_c_               = param.omega_c_;
d_rho_om_c_              = param.rho_om_c_;
d_sig_om_c_              = param.sig_om_c_;
d_tau_c_                 = param.tau_c_;
d_omega_i_               = param.omega_i_;
d_omega_g_               = param.omega_g_;
d_omega_ns_              = param.omega_ns_;
d_omega_nd_              = param.omega_nd_;
d_ii_weights_            = param.ii_weights_;
d_omega_y_               = param.omega_y_;
d_rho_a_                 = param.rho_a_;
d_rho_mu_                = param.rho_mu_;
d_sig_a_                 = param.sig_a_;
d_sig_mu_                = param.sig_mu_;
d_thet_                  = param.thet_;
d_chi_p_                 = param.chi_p_;
d_ss_a_                  = param.a_ss_;
d_ss_mu_                 = param.mu_ss_;
d_c_mu_mm_               = param.c_mu_mm_;
d_c_tfp_mm_              = param.c_tfp_mm_;


% An individual steady state parameter

d_ss_c_r                 = ss.c_r;
d_ss_g                   = ss.g;
d_ss_rk                  = ss.rk;
d_ss_gam_I               = ss.gam_I;
d_ss_gam_G               = ss.gam_G;
d_ss_c_nr                = ss.c_nr;
d_ss_c                   = ss.c;
d_ss_lam_r               = ss.lam_r;
d_ss_lam_nr              = ss.lam_nr;
d_ss_lam                 = ss.lam;
d_ss_lab                 = ss.lab;
d_ss_trans               = ss.trans;
d_ss_invest              = ss.invest;
d_ss_r                   = ss.r;
d_ss_ngdp_p              = ss.ngdp_p;
d_ss_y_va                = ss.y_va;
d_ss_ngdp_e              = ss.ngdp_e;
d_ss_w                   = ss.w;
d_ss_hrs                 = ss.hrs;
d_ss_ngdp_i              = ss.ngdp_i;
d_ss_r                   = ss.r;
d_ss_rr                  = ss.rr;

% A sectoral steady state parameter
d_ss_w_                  = ss.w_;
d_ss_x_                  = ss.x_;
d_ss_gam_                = ss.gam_;
d_ss_gam_x_              = ss.gam_x_;
d_ss_c_                  = ss.c_;
d_ss_gam_f_              = ss.gam_f_;
d_ss_xx_                 = ss.xx_;
d_ss_f_                  = ss.f_;
d_ss_n_                  = ss.n_;
d_ss_k_                  = ss.k_;
d_ss_ks_                 = ss.ks_;
d_ss_i_                  = ss.i_;
d_ss_g_                  = ss.g_;
d_ss_y_                  = ss.y_;
d_ss_mc_                 = ss.mc_;
d_ss_rk_                 = ss.rk_;
d_ss_z_                  = ss.z_;
d_ss_nva_                = ss.nva_;
d_ss_y_va_               = ss.y_va_;
d_ss_gam_va_             = ss.gam_va_;
d_ss_omega_c_            = ss.omega_c_;




% Set parameters to be read by Dynare
eta                      = d_eta;
eta_g                    = d_eta_g;
zeta                     = d_zeta;
vphi                     = d_vphi;
psi                      = d_psi;
xi                       = d_xi;
bet                      = d_bet;
h                        = d_h;
capup                    = d_capup;
eps_w                    = d_eps_w;
nu                       = d_nu;
delt                     = d_delt;
Spp                      = d_Spp;
chi_w                    = d_chi_w;
rho_r                    = d_rho_r;
phi_pie                  = d_phi_pie;
phi_pie_core             = d_phi_pie_core;
phi_pie_NYFed            = d_phi_pie_NYFed;
phi_pie_bar              = d_phi_pie_bar;
phi_pie_forecast         = d_phi_pie_forecast;
phi_pie_headmin          = d_phi_pie_headmin;
phi_y                    = d_phi_y;
phi_gap                  = d_phi_gap;
sig_r                    = d_sig_r;
rho_g                    = d_rho_g;
sig_g                    = d_sig_g;
ss_n                     = d_ss_n;
A_N                      = d_A_N;
thet_w                   = d_thet_w;
ss_G                     = d_ss_G;
rk_wedge                 = d_rk_wedge;
ss_cshift                = d_ss_cshift;
rho_cshift               = d_rho_cshift;
sig_cshift               = d_sig_cshift;
omega_r                  = d_omega_r;
rho_trans                = d_rho_trans;
sig_trans                = d_sig_trans;
lva_weights              = d_lva_weights;
va_weights               = d_va_weights;
consumption_share        = d_consumption_share;
investment_share         = d_investment_share;
government_share         = d_government_share;
employment_share         = d_employment_share;
ss_c_r                   = d_ss_c_r;
ss_g                     = d_ss_g;
ss_rk                    = d_ss_rk;
ss_gam_I                 = d_ss_gam_I;
ss_gam_G                 = d_ss_gam_G;
ss_c_nr                  = d_ss_c_nr;
ss_c                     = d_ss_c;
ss_lam_r                 = d_ss_lam_r;
ss_lam_nr                = d_ss_lam_nr;
ss_lam                   = d_ss_lam;
ss_lab                   = d_ss_lab;
ss_trans                 = d_ss_trans;
ss_invest                = d_ss_invest;
ss_r                     = d_ss_r;
ss_ngdp_p                = d_ss_ngdp_p;
ss_y_va                  = d_ss_y_va;
ss_ngdp_e                = d_ss_ngdp_e;
ss_w                     = d_ss_w;
ss_hrs                   = d_ss_hrs;
ss_ngdp_i                = d_ss_ngdp_i;
ss_r                     = d_ss_r;
ss_rr                    = d_ss_rr;
n_ait                    = d_n_ait;
rho_tfp                  = d_rho_tfp;
sig_tfp                  = d_sig_tfp;
sig_mm                   = d_sig_mm;

paramlist = {'omega_c_', 'sig_om_c_', 'rho_om_c_', 'sig_om_c_', 'tau_c_', ...
    'omega_i_', 'omega_g_', 'omega_ns_', 'omega_nd_', 'omega_y_', 'rho_a_', ...
    'rho_mu_', 'sig_a_', 'sig_mu_', 'thet_', 'chi_p_', 'ss_a_', 'ss_mu_', ...
    'ss_w_', 'ss_x_', 'ss_omega_y_', 'ss_omega_ns_', 'ss_gam_', 'ss_gam_x_', ...
    'ss_c_', 'ss_gam_f_', 'ss_n_', 'ss_k_', 'ss_ks_', 'ss_i_', 'ss_g_', 'ss_y_', ...
    'ss_a_', 'ss_mc_', 'ss_rk_', 'ss_z_', 'ss_nva_', 'ss_y_va_', 'ss_gam_va_', ...
    'ss_omega_c_', 'c_mu_mm_', 'c_tfp_mm_'};

paramlist_ii_ = {'ii_weights_', 'ss_xx_'};

for se = 1 : nsectors

    for parm = 1:length(paramlist);


        eval([char(paramlist(parm)), num2str(se), ' = d_', char(paramlist(parm)), '(', num2str(se), ',1);' ])
    end


    for se2 = 1 : nsectors
        for parm2 = 1 : length(paramlist_ii_);
            eval([char(paramlist_ii_(parm2)), num2str(se), '_', num2str(se2) ' = d_', char(paramlist_ii_(parm2)), '(', num2str(se), ',' num2str(se2), ');' ])

        end
    end

end

% for se = 1 : nsectors
%
% end

%
%
%
% for se = 1 : nsectors
%
%     ss_nva_@{se}         = d_ss_nva_(@{se});
%     ss_y_va_@{se}        = d_ss_y_va_(@{se});
%     ss_gam_va_@{se}      = d_ss_gam_va_(@{se});
%     ss_omega_c_@{se}     = d_ss_omega_c_(@{se});
%     c_mu_mm_@{se}        = d_c_mu_mm_(@{se});
%     c_tfp_mm_@{se}       = d_c_tfp_mm_(@{se});
%
%     for se2 in 1 : nsector
%         ii_weights_@{se}_@{se2}     = d_ii_weights_(@{se},@{se2});
%         ss_xx_@{se}_@{se2}          = d_ss_xx_(@{se},@{se2});
%     end
%
% endfor

