% calc_ss
%
% Function to calculate the steady state of the system, given sectoral productivity
% levels.
% 
% =========================================================================
function out = calc_ss(SET) 

param    = SET.param ;
nsectors = SET.nsectors ;

% Obvious steady state values
r   = 1 / param.bet ; 
rk  = param.rk_wedge*(1 / param.bet - 1 + param.delt) ;

% Load initial parameter guess
x0 = SET.x0 ;

% Many, many evaluations
SET.options.MaxFunctionEvaluations=1000000;
SET.options.Algorithm = "trust-region-dogleg";

% Call minimiser
[x_ss, fval, exitflag] = fsolve(@ssfunc, x0, SET.options) ;

% if exitflag ==0
%     disp('Warning: exitflag equal to 0, calibration might be wrong') ; 
%     fval
% end

if max(abs(fval))>0.1
    disp('Error too big: max error is')
    max(abs(fval))
    end

out = exp(x_ss) ;

% ========================================================================
% FUNCTIONS
% ========================================================================
function F = ssfunc(x0)
    
x0 = exp(x0) ;

x0(end+1)=1;

j_ = 1 ;
    
c_r     = x0(j_) ; j_ = j_ + 1 ;
g       = x0(j_) ; j_ = j_ + 1 ;

for m_ = 1 : nsectors

    w_(m_,1)      = x0(j_) ; j_ = j_ + 1 ;
    x_(m_,1)      = x0(j_) ; j_ = j_ + 1 ;

end

gam_ = x0(j_ : end) ;
    
gam_(end) = ((1 - sum(param.omega_c_(1:end-1).*(gam_(1:end-1).^(1-param.eta))))./param.omega_c_(end))^(1/(1-param.eta)) ;

% Relative prices of investment and public demand and intermediates, given goods prices
gam_I  = (sum(param.omega_i_ .* (gam_.^(1-param.eta))))^(1/(1-param.eta));
gam_G  = (sum(param.omega_g_ .* (gam_.^(1-param.eta_g))))^(1/(1-param.eta));
gam_x_ =  (param.ii_weights_' * (gam_.^(1 - param.psi))).^(1/(1-param.psi)) ;   % Checked the matrix multiplication on this

% Assume that transfers equalise consumption
c_nr   = c_r ;

% Aggregate consumption
c   = param.omega_r * c_r + (1 - param.omega_r) * c_nr ;

% MUC for Ricardian consumers, given c_r
lam_r   = param.cshift_ss * (1 - param.h * param.bet) / (1 - param.h) / c_r ;
lam_nr  = lam_r ;
lam     = param.omega_r * lam_r + (1 - param.omega_r) * lam_nr ; 

% Consumption of individual varieties, given prices and aggregate consumption
c_ = param.omega_c_ .* (gam_.^(-param.eta)) .* c ;

% Price of value-added given factor prices
gam_f_= (param.omega_nd_ .* w_.^(1-param.zeta) ...
            + ( 1- param.omega_nd_) .* rk.^(1-param.zeta)).^(1/(1-param.zeta)) ; ...

% Demand for intermediates, given prices and intermediates demand
xx_ = param.ii_weights_ .* (gam_ * (1./gam_x_)').^(-param.psi) * diag(x_) ;

% Value added, given prices and intermediates demand
f_ = param.omega_y_./( 1- param.omega_y_) .* ...
              gam_x_.^param.vphi .* x_ ./ gam_f_.^param.vphi; ...        
        
% Factor demands, given relative prices and value added          
n_ = param.omega_nd_ .* (w_ ./ gam_f_).^(-param.zeta) .* f_; 
k_ = (1 - param.omega_nd_) .* (rk ./ gam_f_).^(-param.zeta) .* f_ ;
        
% Investment and government spending demand, given aggregate demand and relative prices 
i_ = param.omega_i_.* (gam_ ./ gam_I).^(-param.eta) * sum(param.delt * k_) ;
g_ = param.omega_g_.* (gam_ ./ gam_G).^(-param.eta_g) .* g ;

% Gross output
y_ = c_ + i_ + g_ + sum(xx_,2) ; 

% Aggregate labour
lab =  (sum(param.omega_ns_.^(-1/param.xi).*n_.^((param.xi+1)/param.xi))) ^(param.xi/(param.xi+1));

% Productivity

% Value added
y_va = sum(y_) - sum(x_) ;


%%

% *** First all the steady state conditions including only aggregate variables ***
F = [ 

% Intermediates must satisfy production function
y_ - param.a_ss_ .* (param.omega_y_.^(1/param.vphi) .* f_ ...
     .^((param.vphi-1)/param.vphi) + (1 - param.omega_y_).^(1/param.vphi) .* x_ ...
     .^((param.vphi-1)/param.vphi)).^(param.vphi/(param.vphi-1)) ; ...

% Wages must satisfy labour supply condition
param.omega_ns_.^(-1/param.xi) .* param.A_N .* lab.^(param.nu - 1 / param.xi) ...
            .* n_.^(1/param.xi) - lam .* w_; ...

% Check that MC = 1 holds with new relative prices
param.a_ss_ - (param.omega_y_ .* (gam_f_ .*(1 + param.tau_c_)./ gam_).^(1-param.vphi) +...
     (1 - param.omega_y_) .* (gam_x_ .*(1 + param.tau_c_)./ gam_).^(1-param.vphi)).^(1/(1-param.vphi)) ; 
            
% Government spending must equal right share of GDP
gam_G * g / (c + gam_I * sum(param.delt * k_) + gam_G * g) - param.government_share

] ;

end

end


