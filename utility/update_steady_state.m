% update_steady_state
%
% Update the steady state values
% =================================================================
function out = update_steady_state(SET, ss_params) 

nsectors = SET.nsectors ;
param    = SET.param ;

% Update steady state values:

% === Update steady state values ===
ss_params(end+1)=1 ;
j_ = 1 ;
%ss.lab      = ss_params(j_) ; j_ = j_ + 1 ; 
ss.c_r      = ss_params(j_) ; j_ = j_ + 1 ; 
ss.g        = ss_params(j_) ; j_ = j_ + 1 ;

for m_ = 1 : nsectors 
    ss.w_(m_,1)   = ss_params(j_) ; j_ = j_ + 1 ; 
    ss.x_(m_,1)   = ss_params(j_) ; j_ = j_ + 1 ; 
%    ss.gam_(m_,1) = ss_params(j_) ; j_ = j_ + 1 ;    
end

% for m_ = 1 : nsectors
%     ss.gam_(m_,1) = ss_params(j_) ; j_ = j_ + 1 ;
% end

ss.gam_ = ss_params(j_:end) ;

ss.gam_(end) = ((1 - sum(param.omega_c_(1:end-1).*(ss.gam_(1:end-1).^(1-param.eta))))./...
                param.omega_c_(end))^(1/(1-param.eta)) ;

% === Additional steady state values ===

% Return on capital
ss.rk  = param.rk_wedge*(1 / param.bet - 1 + param.delt) ;

% Relative prices of investment and public demand and intermediates, given goods prices
ss.gam_I  = (sum(param.omega_i_ .* (ss.gam_.^(1-param.eta))))^(1/(1-param.eta));
ss.gam_G  = (sum(param.omega_g_ .* (ss.gam_.^(1-param.eta_g))))^(1/(1-param.eta));
ss.gam_x_ =  (param.ii_weights_' * (ss.gam_.^(1 - param.psi))).^(1/(1-param.psi)) ;   % Checked the matrix multiplication on this

% Assume that transfers equalise consumption
ss.c_nr   = ss.c_r ;

% Aggregate consumption
ss.c   = param.omega_r * ss.c_r + (1 - param.omega_r) * ss.c_nr ;

% MUC for Ricardian consumers, given c_r
ss.lam_r   = param.cshift_ss * (1 - param.h * param.bet) / (1 - param.h) / ss.c_r ;
ss.lam_nr  = ss.lam_r ;
ss.lam     = param.omega_r * ss.lam_r + (1 - param.omega_r) * ss.lam_nr ; 
    
% Consumption of individual varieties, given prices and aggregate consumption
ss.c_ = param.omega_c_ .* (ss.gam_.^(-param.eta)) .* ss.c ;

% Price of value-added given factor prices
ss.gam_f_= (param.omega_nd_ .* ss.w_.^(1-param.zeta) ...
            + ( 1- param.omega_nd_) .* ss.rk.^(1-param.zeta)).^(1/(1-param.zeta)) ; ...

% Demand for intermediates, given prices and intermediates demand
ss.xx_ = param.ii_weights_ .* (ss.gam_ * (1./ss.gam_x_)').^(-param.psi) * diag(ss.x_) ;

% Value added, given prices and intermediates demand
ss.f_ = param.omega_y_./( 1- param.omega_y_) .* ...
              ss.gam_x_.^param.vphi .* ss.x_ ./ ss.gam_f_.^param.vphi; ...

% Factor demands, given relative prices and value added          
ss.n_ = param.omega_nd_ .* (ss.w_ ./ ss.gam_f_).^(-param.zeta) .* ss.f_; 
ss.k_ = (1 - param.omega_nd_) .* (ss.rk ./ ss.gam_f_).^(-param.zeta) .* ss.f_ ;
ss.ks_ = ss.k_ ;
        
% Investment and government spending demand, given aggregate demand and relative prices 
ss.i_ = param.omega_i_.* (ss.gam_ ./ ss.gam_I).^(-param.eta) * sum(param.delt * ss.k_) ;
ss.g_ = param.omega_g_.* (ss.gam_ ./ ss.gam_G).^(-param.eta_g) .* ss.g ;

% Gross output
ss.y_ = ss.c_ + ss.i_ + ss.g_ + sum(ss.xx_,2) ; 

% Aggregate labour
ss.lab =  (sum(param.omega_ns_.^(-1/param.xi).*ss.n_.^((param.xi+1)/param.xi))) ^(param.xi/(param.xi+1));

% Add additional steady states
ss.trans    = ss.c_nr - sum(ss.w_ .* ss.n_) ;

ss.invest   = sum(param.delt * ss.k_) ;     % Aggregate investment
ss.r        = 1 / param.bet ;               % Policy rate
ss.rr       = 1 / param.bet ;               % Real rate (assuming steady state = 1 ;
ss.ngdp_p   = sum(ss.gam_'*ss.y_ - ss.gam_x_'*ss.x_) ;
ss.y_va     = ss.ngdp_p ;                   % Steady state of value added
ss.ngdp_e   = ss.c + ss.gam_I * ss.invest + ...
                ss.gam_G * ss.g;            % Expenditure measure of GDP
ss.w        = 0 ;                           % Initialise aggregate wage index
ss.y_va     = 0 ;
ss.hrs      = 0 ;
ss.mc_      = ones(nsectors,1) ;

for j_ = 1 : nsectors
    
    ss.rk_(j_,1)  = param.rk_wedge*(1 / param.bet - 1 + param.delt) ;
    ss.z_(j_,1)   = param.delt * ss.k_(j_) ;
    ss.w          = ss.w + param.omega_ns_(j_) * ss.w_(j_)^(1 + param.xi) ;
    ss.nva_(j_,1)    = ss.gam_(j_) * ss.y_(j_) - ss.gam_x_(j_) * ss.x_(j_) ;
    ss.y_va_(j_,1)   = ss.y_(j_) - ss.x_(j_) ; 
    ss.gam_va_(j_,1) = ss.nva_(j_,1)./ ss.y_va_(j_,1) ;

%    ss.gam_va_(j_,1) = (ss.gam_(j_) * ss.y_(j_)) / ss.nva_(j_,1) ;
%    ss.y_va_(j_,1)   = ss.nva_(j_,1)./ss.gam_va_(j_,1) ;
%    ss.y_va          = ss.y_va + ss.nva_(j_) ;
    ss.y_va          = ss.y_va + ss.y_va_(j_) ;
    ss.hrs           = ss.hrs + ss.n_(j_) ;
    
end

ss.ngdp_i    = ((1 + param.tau_c_).*ss.w_)' * ss.n_ + ((1 + param.tau_c_).*ss.rk_)' * ss.k_ ...
                + (param.tau_c_.*ss.gam_x_)'*ss.x_;    % Income measure of GDP

ss.w        = ss.w^(1/(1+param.xi)) ;               % Wage index
ss.a_       = SET.param.a_ss_ ;
ss.mu_      = SET.param.mu_ss_ ;

ss.omega_c_ = SET.param.omega_c_ ;

if abs(ss.ngdp_p - ss.ngdp_e)>0.001 || abs(ss.ngdp_p - ss.ngdp_i) > 0.001
    disp('Oh no - steady E, I and P are not the same!')
ss.ngdp_p
ss.ngdp_e
ss.ngdp_i
end

out = ss ;
